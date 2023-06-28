cmake_minimum_required(VERSION 3.17)
project(wav)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc")
set(SFML_STATIC_LIBRARIES TRUE)





add_executable(wav main.cpp)


set(SFML_DIR "D:\\cLion\\programming\\lib\\sfml\\GCC 7.3.0 MinGW (SEH) - 64-bit\\SFML-2.5.1\\lib\\cmake\\SFML")
target_link_libraries(wav sfml-graphics)
target_link_libraries(wav sfml-audio)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_modules")
find_package(SFML 2.5 REQUIRED graphics network audio window system)
include_directories("D:\\cLion\\programming\\lib\\sfml\\GCC 7.3.0 MinGW (SEH) - 64-bit\\SFML-2.5.1\\include\\SFML")

include_directories("D:\\cLion\\programming\\lib\\lame\\lame-3.100\\include")
