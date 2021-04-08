add_executable(
  ActsExampleFatrasDD4hep
  DD4hepFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasDD4hep
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorDD4hep)

install(
  TARGETS ActsExampleFatrasDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
add_executable(
  ActsExampleFatrasDD4hep2
  DD4hepFatrasExample2.cpp)
target_link_libraries(
  ActsExampleFatrasDD4hep2
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorDD4hep ActsExamplesEventRecording ActsExamplesGeant4 ActsExamplesGeant4DD4hep)

install(
  TARGETS ActsExampleFatrasDD4hep2
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
