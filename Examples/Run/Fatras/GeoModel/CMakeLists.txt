add_executable(
  ActsExampleFatrasGeoModel
  GeoModelFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasGeoModel
  PRIVATE ActsExamplesFatrasCommon ActsExamplesEventRecording ActsExamplesGeant4 ActsExamplesGeant4GeoModel ActsExamplesDetectorTGeo)

install(
  TARGETS ActsExampleFatrasGeoModel
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
