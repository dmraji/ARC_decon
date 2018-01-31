cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# setting a variable equal to the executable name
set(name decon_main)

# name the project
project(decon_main)

# get the implementation files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

# get the header files
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# create the executable
add_executable(${name} ${name}.cc ${sources} ${headers})

# get the external libraries
link_directories(/usr/lib/x86_64-linux-gnu)
target_link_libraries(${name} pqxx-4.0.so armadillo.so)
