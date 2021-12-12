cmake_minimum_required(VERSION 3.20)
project(new)

set(CMAKE_CXX_STANDARD 11)
include_directories(./include/json)

#link_directories(D:/new)
#add_library(lib libjsoncpp-25.dll)
#set_target_properties(lib PROPERTIES LINKER_LANGUAGE CPP)

#add_library(bar SHARED IMPORTED)
#set_property(TARGET bar PROPERTY
#        IMPORTED_LOCATION "D:/new/libjsoncpp-25.dll")
#set_property(TARGET bar PROPERTY
#        IMPORTED_IMPLIB "D:/new/libjsoncpp.dll.a")

add_executable(new main.cpp)

#set_target_properties(new PROPERTIES LINKER_LANGUAGE C)
#target_include_directories(new PUBLIC C:/Users/'Confidence'F/Desktop/new/include/json)
target_link_libraries(new PUBLIC "C:/Users/'Confidence'F/Documents/dev netlist/new/libjsoncpp-25.dll")


#target_link_libraries(new lib)

#target_link_libraries(new PRIVATE bar)

