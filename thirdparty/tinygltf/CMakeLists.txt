set(TINYGLTF_USE_RAPIDJSON ON CACHE INTERNAL "" FORCE)
add_library(tinygltf INTERFACE)
target_include_directories(tinygltf INTERFACE .)

