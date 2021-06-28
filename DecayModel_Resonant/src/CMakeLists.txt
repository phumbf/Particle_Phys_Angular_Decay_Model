add_library(ResonantDecayModelLib Event.cc FitParameter.cc HelperFunctions.cc KinematicFunctions.cc KPiResonanceBase.cc KPiResonanceBW.cc KPiResonanceBWNorm.cc KPiResonanceFlatte.cc KPiResonanceLASS.cc KPiResonanceNR.cc OptionsParser.cc Parser.cc PsiKPiAmplitudeModel.cc PsiKPiFitter.cc PsiPiResonance.cc SplineModel.cc Variables.cc Wignerd.cc )

target_link_libraries(ResonantDecayModelLib PUBLIC root_interface)
target_include_directories(ResonantDecayModelLib PUBLIC ${CMAKE_SOURCE_DIR}/include)

set_target_properties(ResonantDecayModelLib PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})


