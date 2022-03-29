add_library(matrix_operations STATIC matrix_operations.cpp ../include/matrix_operations.h)
target_include_directories(matrix_operations PUBLIC ../include/)

add_library(math STATIC math.cpp ../include/math.h)
target_include_directories(math PUBLIC ../include/)

add_library(multipole_conv_funcs STATIC multipole_conv_funcs.cpp
  ../include/multipole_conv_funcs.h)
# inherits include directory from math 
target_link_libraries(multipole_conv_funcs PUBLIC math) 
