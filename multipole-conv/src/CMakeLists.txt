add_library(math STATIC math.cpp ../include/math.h)
target_include_directories(math PUBLIC ../include)

# header only libraries
add_library(options INTERFACE ../include/options.h)
target_include_directories(options INTERFACE ../include)

add_library(square_matrix INTERFACE ../include/square_matrix.h)
target_include_directories(square_matrix INTERFACE ../include)
