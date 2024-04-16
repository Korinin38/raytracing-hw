#pragma once

#include <core/scene.h>

#include <string>

Scene parse_scene_gltf(const std::string &filename, int width, int height, int samples);
