get_directory_property(incdirs INCLUDE_DIRECTORIES)
message(STATUS "Include Paths: ${incdirs}")

foreach(idir ${AINCDIRS})
  message(STATUS "Adding include dir to Eunice's project.. ${idir}")
  include_directories(${idir})
endforeach()

add_library(Smearception SHARED Smearception.cxx Smearception.hxx)

add_executable(TestSmear TestSmear.cxx)

target_link_libraries(TestSmear reweight;${AEVENT_LIBS};vali;Smearception)
set(FOMP "-fopenmp")
set_target_properties(TestSmear PROPERTIES LINK_FLAGS "${CMAKE_LINK_FLAGS} ${FOMP}")

add_executable(ValidRejSampling ValidRejSampling.cxx Smearception.cxx Smearception.hxx)

target_link_libraries(ValidRejSampling reweight;${AEVENT_LIBS};vali)
set(FOMP "-fopenmp")
set_target_properties(ValidRejSampling PROPERTIES LINK_FLAGS "${CMAKE_LINK_FLAGS} ${FOMP}")

install (TARGETS TestSmear ValidRejSampling DESTINATION bin)
install (TARGETS Smearception DESTINATION lib)
