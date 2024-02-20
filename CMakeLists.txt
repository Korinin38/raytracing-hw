cmake_minimum_required(VERSION 3.16)
project(solution)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(solution
        src/main.cpp
        src/core/scene.h
        src/core/scene.cpp
        src/utils/base.h
        src/geometry/primitive.h
        src/utils/base.cpp
        src/geometry/plane.h
        src/geometry/objects.h
        src/geometry/ellipsoid.h
        src/geometry/box.h
        src/geometry/primitive.cpp
        src/core/camera.cpp
        src/core/camera.h
        src/render/canvas.cpp
        src/render/canvas.h
        src/geometry/light_source.cpp
        src/geometry/light_source.h
        src/utils/timer.h)
