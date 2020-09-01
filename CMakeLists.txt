cmake_minimum_required(VERSION 3.6)
project(waveformtools CXX)
find_package( ZLIB REQUIRED )
include_directories($ENV{BOOST_INC})
include_directories($ENV{CANVAS_INC})
include_directories($ENV{CANVAS_ROOT_IO_INC})
include_directories($ENV{CETLIB_EXCEPT_INC})
include_directories($ENV{CETLIB_INC})
include_directories($ENV{FHICLCPP_INC})
include_directories($ENV{GALLERY_INC})
include_directories($ENV{LARCOREOBJ_INC})
include_directories($ENV{LARDATAOBJ_INC})
include_directories($ENV{NUSIMDATA_INC})
include_directories($ENV{ROOT_INC})
include_directories($ENV{HEP_CONCURRENCY_INC})
link_directories($ENV{ROOTSYS}/lib $ENV{CANVAS_LIB} $ENV{CETLIB_EXCEPT_LIB} $ENV{CETLIB_LIB} $ENV{GALLERY_LIB} $ENV{NUSIMDATA_LIB} $ENV{LARCOREOBJ_LIB} $ENV{LARDATAOBJ_LIB} $ENV{BOOST_LIB})

set(MY_LIBS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread MultiProc pthread canvas cetlib_except cetlib gallery nusimdata_SimulationBase larcoreobj_SummaryData lardataobj_RecoBase lardataobj_RawData boost_program_options ${ZLIB_LIBRARIES})

add_executable(extract_larsoft_waveforms extract_larsoft_waveforms.cxx cnpy.cpp)
set_property(TARGET extract_larsoft_waveforms PROPERTY CXX_STANDARD 17)
target_link_libraries(extract_larsoft_waveforms ${MY_LIBS})

add_executable(extract_larsoft_hits extract_larsoft_hits.cxx cnpy.cpp)
set_property(TARGET extract_larsoft_hits PROPERTY CXX_STANDARD 17)
target_link_libraries(extract_larsoft_hits  ${MY_LIBS})

add_executable(extract_photon_waveforms extract_photon_waveforms.cxx cnpy.cpp)
set_property(TARGET extract_photon_waveforms PROPERTY CXX_STANDARD 17)
target_link_libraries(extract_photon_waveforms ${MY_LIBS})

add_executable(extract_larsoft_wires extract_larsoft_wires.cxx cnpy.cpp)
set_property(TARGET extract_larsoft_wires PROPERTY CXX_STANDARD 17)
target_link_libraries(extract_larsoft_wires ${MY_LIBS})

add_executable(extract_larsoft_simch extract_larsoft_simch.cxx cnpy.cpp)
set_property(TARGET extract_larsoft_simch PROPERTY CXX_STANDARD 17)
target_link_libraries(extract_larsoft_simch ${MY_LIBS})


add_subdirectory(test)
