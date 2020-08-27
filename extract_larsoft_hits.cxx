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
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RecoBase/Hit.h"

#include "cnpy.h"

using namespace art;
using namespace std;
using namespace std::chrono;

namespace po = boost::program_options;

enum class Format { Text, Numpy };

// DecoderandReco | timingrawdecoder | daq.................. | std::vector<raw::RDTimeStamp>....................................... | ....1

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
// raw waveforms are written to `outfile`, while the true energy
// depositions are written to `truth_outfile` (unless it is an empty
// string, in which case no truth file is produced). If `onlySignal`
// is true, then only channels with some true energy deposition are
// written out; otherwise all channels are written out.
// 
// Each line in `outfile` has the format:
//
// event_no channel_no sample_0 sample_1 ... sample_N
//
// Each line in `truth_outfile` has the format
//
// event_no channel_no tdc total_charge
void
extract_larsoft_hits(std::string const& tag,
                     std::string const& filename,
                     std::string const& outfile,
                     Format format,
                     int nevents, int nskip,
                     int triggerType)
{
    InputTag daq_tag{ tag };
    // Create a vector of length 1, containing the given filename.
    vector<string> filenames(1, filename);

    int iev=0;
    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
        vector<vector<int> > samples;

        if(iev<nskip) continue;
        if(iev>=nevents+nskip) break;
        if(triggerType!=-1){
            auto& timestamp=*ev.getValidHandle<std::vector<raw::RDTimeStamp>>(InputTag{"timingrawdecoder:daq:DecoderandReco"});
            assert(timestamp.size()==1);
            if(timestamp[0].GetFlags()!=triggerType){
                std::cout << "Skipping event " << ev.eventAuxiliary().event()  << " with trigger type " << timestamp[0].GetFlags() << std::endl;
                continue;
            }
            else{
                std::cout << "Using event " << ev.eventAuxiliary().event()  << " with trigger type " << timestamp[0].GetFlags() << std::endl;
            }
        }
        std::cout << "Event " << ev.eventAuxiliary().id() << std::endl;

        //------------------------------------------------------------------
        // Look at the hits
        auto& hits =
            *ev.getValidHandle<vector<recob::Hit>>(daq_tag);
        if(hits.empty()){
            std::cout << "Hits vector is empty" << std::endl;
        }
        for(auto&& hit: hits){
            samples.push_back({
                    (int)hit.Channel(), hit.StartTick(), hit.EndTick(), (int)hit.SummedADC(), (int)hit.RMS()
                        });
        } // end loop over digits (=?channels)
        std::string this_outfile(outfile);
        size_t dotpos=outfile.find_last_of(".");
        if(dotpos==std::string::npos){
            dotpos=outfile.length();
        }
        std::ostringstream iss;
        //auto& rdtimestamps=*ev.getValidHandle<std::vector<raw::RDTimeStamp>>(InputTag{"timing:daq:RunRawDecoder"});
        //assert(rdtimestamps.size()==1);

        iss << outfile.substr(0, dotpos) << "_evt" << ev.eventAuxiliary().event() << "_t0x" << outfile.substr(dotpos, outfile.length()-dotpos);
        std::cout << "Writing event " << ev.eventAuxiliary().event() << " to file " << iss.str() << std::endl;
        save_to_file<int>(iss.str(), samples, format, false);
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
        ("tag,g", po::value<string>()->default_value("daq"), "input tag (aka \"module label\") of input digits")
        ("nevent,n", po::value<int>()->default_value(1), "number of events to save")
        ("nskip,k", po::value<int>()->default_value(0), "number of events to skip")
        ("numpy", "use numpy output format instead of text")
        ("trig", po::value<int>()->default_value(-1), "select events with given trigger type")
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

    extract_larsoft_hits(vm["tag"].as<string>(),
                         vm["input"].as<string>(),
                         vm["output"].as<string>(),
                         vm.count("numpy") ? Format::Numpy : Format::Text,
                         vm["nevent"].as<int>(),
                         vm["nskip"].as<int>(),
                         vm["trig"].as<int>());
    return 0;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
