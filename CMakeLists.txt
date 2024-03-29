cmake_minimum_required(VERSION 3.16)
project(solution)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/src )

find_package(OpenMP)

option(USE_CALLBACK "display progress and time" OFF)

add_executable(solution
        src/main.cpp
        src/core/scene.h
        src/core/scene.cpp
        src/geometry/primitive.h
        src/geometry/primitive.cpp
        src/geometry/light_source.h
        src/geometry/light_source.cpp
        src/core/camera.cpp
        src/core/camera.h
        src/render/canvas.cpp
        src/render/canvas.h
        src/utils/base.h
        src/utils/base.cpp
        src/utils/timer.h
        src/utils/vector.h
        src/utils/random.h
        src/utils/random.cpp)

if (USE_CALLBACK)
    target_compile_definitions(solution PRIVATE "USE_CALLBACK")
endif()

target_link_libraries(solution PRIVATE OpenMP::OpenMP_CXX)