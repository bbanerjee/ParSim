# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics/ParticleCreator

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_executable(utTest0 utTest0.cpp)
target_link_libraries(utTest0 ${Boost_LIBRARIES} ${MPI_LIBRARIES})
add_test(utTest0 ./utTest0)
