add_library(matrix_operations STATIC matrix_operations.cpp
  ../include/matrix_operations.h)
target_include_directories(matrix_operations PUBLIC ../include/)

add_library(math STATIC math.cpp ../include/math.h)
target_include_directories(math PUBLIC ../include/)

add_library(multipole_conv_funcs STATIC multipole_conv_funcs.cpp
  ../include/multipole_conv_funcs.h)
# inherits include directory from math 
target_link_libraries(multipole_conv_funcs PUBLIC math) 

find_package(Boost REQUIRED COMPONENTS
             program_options)
add_library(command_line_parser STATIC command_line_parser.cpp
  ../include/command_line_parser.h)
target_include_directories(command_line_parser PUBLIC ../include/
  ${PROJECT_BINARY_DIR})
target_link_libraries(command_line_parser Boost::program_options)
