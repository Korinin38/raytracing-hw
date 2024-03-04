cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 11)

add_library(pcg
            pcg/include/pcg_random.hpp
            pcg/include/pcg_extras.hpp
            pcg/include/pcg_uint128.hpp)

set_target_properties(pcg PROPERTIES LINKER_LANGUAGE CXX)
include_directories(${CMAKE_SOURCE_DIR}/pcg/include)