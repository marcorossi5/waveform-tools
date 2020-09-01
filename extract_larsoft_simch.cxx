#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "boost/program_options.hpp"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "cnpy.h"

using namespace art;
using namespace std;
using namespace std::chrono;

namespace po = boost::program_options;

enum class Format { Text, Numpy };


template<class T>
void save_to_file(std::string const& outfile,
                  std::vector<std::vector<T> > v,
                  Format format,
                  bool append)
{
    switch(format){

    case Format::Text:
    {
        // Open in append mode because we 
        std::ofstream fout(outfile, append ? ios::app : ios_base::out);
        for(auto const& v1 : v){
            for(auto const& s : v1){
                fout << s << " ";
            }
            fout << std::endl;
        }
    }
    break;
    
    case Format::Numpy:
    {
        // Do nothing if the vector is empty
        if(v.empty() || v[0].empty()) break;
        // cnpy needs a single contiguous array of data, so do that conversion
        std::vector<T> tmp;
        tmp.reserve(v.size()*v[0].size());
        for(auto const& v1 : v){
            for(auto const& s : v1){
                tmp.push_back(s);
            }
        }
        cnpy::npy_save(outfile, &tmp[0], {v.size(), v[0].size()}, append ? "a" : "w");
    }
    break;
    }
}

// Write `nevents` events of data from `filename` to text files. The
// sim::SimChannels are written to `outfile`.
// 
// Each line in `outfile` has the format:
//
// event_no, channel_no, tdc, t_ID, charge, energy, x, y, z 

void
extract_larsoft_waveforms(std::string const& tag,
                          std::string const& filename,
                          std::string const& outfile,
                          Format format,
                          int nevents, int nskip,
                          int triggerType,
                          bool timestampInFilename)
{
    InputTag simch_tag{ tag };
    std::string ext = (format==Format::Text) ? ".txt" : ".npy";
    // Create a vector of length 1, containing the given filename.
    vector<string> filenames(1, filename);

    int iev=0;
    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
        vector<vector<int> > samples;
        vector<vector<float> > trueIDEs;

        std::set<int> channelsWithSignal;
        if(iev<nskip) continue;
        if(iev>=nevents+nskip) break;
        if(triggerType!=-1){
            auto& timestamp=*ev.getValidHandle<std::vector<raw::RDTimeStamp>>(
                            InputTag{"timingrawdecoder:daq:DecoderandReco"});
            assert(timestamp.size()==1);
            if(timestamp[0].GetFlags()!=triggerType){
                std::cout << "Skipping event " << ev.eventAuxiliary().event()
                          << " with trigger type " << timestamp[0].GetFlags()
                          << std::endl;
                
                continue;
            }
            else{
                std::cout << "Using event " << ev.eventAuxiliary().event()
                          << " with trigger type " << timestamp[0].GetFlags()
                          << std::endl;
            }
        }
        std::cout << "Event " << ev.eventAuxiliary().id() << std::endl;

        // Get the SimChannels so we can see
        //where the actual energy depositions were
        auto& simchs=*ev.getValidHandle<std::vector<sim::SimChannel>>(
                     InputTag{"tpcrawdecoder:simpleSC:Detsim"});
        
        for(auto&& simch: simchs){
            channelsWithSignal.insert(simch.Channel());

            double charge=0;
            int tID;
            double energy;
            for (const auto& TDCinfo: simch.TDCIDEMap()) {
                auto const tdc = TDCinfo.first;
                for (const sim::IDE& ide: TDCinfo.second) {
                    tID = ide.trackID;
                    charge = ide.numElectrons;
                    energy = ide.energy;
                    trueIDEs.push_back(std::vector<float>{(float)iev,
                                                          (float)simch.Channel(),
                                                          (float)tdc,
                                                          (float)tID,
                                                          (float)charge,
                                                          (float)energy,
                                                          (float)ide.x,
                                                          (float)ide.y,
                                                          (float)ide.z});
                } // for IDEs
            } // for TDCs
        } // loop over SimChannels

        std::string this_outfile(outfile);
        size_t dotpos=outfile.find_last_of(".");
        if(dotpos==std::string::npos){
            dotpos=outfile.length();
        }
        std::ostringstream iss, timestampStr;
        if(timestampInFilename){
            auto& rdtimestamps=*ev.getValidHandle<std::vector<raw::RDTimeStamp>>(
                                InputTag{"timing:daq:RunRawDecoder"});
            assert(rdtimestamps.size()==1);
            // std::cout << "eventAuxiliary value is " << ev.eventAuxiliary().time().value() << std::endl;
            timestampStr << "_t0x" << std::hex << rdtimestamps[0].GetTimeStamp();
        }
        iss << outfile.substr(0, dotpos) << "_evt" << ev.eventAuxiliary().event()
            << timestampStr.str() <<  outfile.substr(dotpos, outfile.length()-dotpos) << ext;
        std::cout << "Writing event " << ev.eventAuxiliary().event()
                  << " to file " << iss.str() << std::endl;
        save_to_file<float>(iss.str(), trueIDEs, format, iev!=0);
        ++iev;
    } // end loop over events
}

int main(int argc, char** argv)
    {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<string>(), "input file name")
        ("output,o", po::value<string>(), "base output file name. Individual output files will be created for each event, with \"_evtN\" inserted before the extension, or at the end if there is no extension")
        ("tag, g", po::value<string>()->default_value("tpcrawdecoder:simpleSC:DetsimStage1"),
                   "input tag (aka \"module label:product instance name: process name\") for sim::SimChannels")
        ("nevent,n", po::value<int>()->default_value(1), "number of events to save")
        ("nskip,k", po::value<int>()->default_value(0), "number of events to skip")
        ("numpy", "use numpy output format instead of text")
        ("trig", po::value<int>()->default_value(-1), "select events with given trigger type")
        ("ts", "add event timestamp to filename")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if(vm.count("help") || vm.empty()) {
        cout << desc << "\n";
        return 1;
    }

    if(!vm.count("input")){
        cout << "No input file specified" << endl;
        cout << desc << endl;
        return 1;
    }

    if(!vm.count("output")){
        cout << "No output file specified" << endl;
        cout << desc << endl;
        return 1;
    }

    extract_larsoft_waveforms(vm["tag"].as<string>(),
                              vm["input"].as<string>(),
                              vm["output"].as<string>(),
                              vm.count("numpy") ? Format::Numpy : Format::Text,
                              vm["nevent"].as<int>(),
                              vm["nskip"].as<int>(),
                              vm["trig"].as<int>(),
                              vm.count("ts"));
    return 0;
}
