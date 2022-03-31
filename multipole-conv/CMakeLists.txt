cmake_minimum_required(VERSION 3.18.0)

project( multipole-conv VERSION 0.5.0 DESCRIPTION "A tool to convert between spherical mutlipole moments and Cartesian multipole moments." LANGUAGES CXX )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server (e.g.
				      # to recognise the include libraries)
add_subdirectory(src)

add_executable(multipole-conv src/multipole-conv.cpp)
 # inherits the include/ directory from multipole_conv_funcs
 target_link_libraries(multipole-conv PRIVATE multipole_conv_funcs
   matrix_operations command_line_parser)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
