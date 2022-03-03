add_library(math STATIC math.cpp ../include/math.h)
target_include_directories(math PUBLIC ../include)

# header only file
add_library(options INTERFACE ../include/options.h)
target_include_directories(options INTERFACE ../include)

add_executable(multipole-conv multipole-conv.cpp)
target_link_libraries(multipole-conv PRIVATE math options)
