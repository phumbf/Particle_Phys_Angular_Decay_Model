add_library(RareDecayModelLib SHARED Amplitudes.cpp  EvtGenFormFactor.cpp  FFcalculation.cpp	FFparser.cpp	Kinematic.cpp   RandomSampling.cpp   StatTools.cpp EvtGenAmplitudes.cpp  EvtGenWilson.cpp  HistogramInterval.cpp Parameters.cpp SphericalHarmonics.cpp  Wilson.cpp LineShapes.cpp ResonanceBase.cpp ResonanceParser.cpp AmplitudeSet.cpp ScalarResonance.cpp VectorResonance.cpp HHLineShapeBase.cpp generator.cpp)

target_link_libraries(RareDecayModelLib PUBLIC root_interface)
target_include_directories(RareDecayModelLib PUBLIC ${CMAKE_SOURCE_DIR}/include)

set_target_properties(RareDecayModelLib PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
