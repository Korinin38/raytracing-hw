#pragma once

#include <string>
#include <core/scene.h>

Scene parse_scene_naive(const std::string &filename);
Scene parse_scene_gltf(const std::string &filename);
