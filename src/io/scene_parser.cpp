#include "scene_parser.h"

#include <utils/matrix.h>

#include <cstdlib>
#include <iostream>

// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <tinygltf/tiny_gltf.h>

struct ScenePartial {
    camera_uniq_ptr camera = nullptr;
    vector3f bg_color{};
    int ray_depth = 6;
    int samples = 256;
    std::vector<primitive_sh_ptr> objects;
    float max_distance = 1e9;
};

Scene parse_scene_gltf(const std::string &filename, int width, int height, int samples) {
    ScenePartial scene;

    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, filename);

    if (!warn.empty()) {
        std::cout << "[tinygltf] Warning: " << warn << std::endl;
    }

    if (!err.empty()) {
        std::cout << "[tinygltf] Error: " << err << std::endl;
    }

    if (!ret)
        throw std::runtime_error("[tinygltf] Failed to parse \"" + filename + "\"");

    int node_with_camera = -1;

    // calculate transformation for every node
    for (int i = 0; i < model.nodes.size(); ++i) {
        tinygltf::Node &node = model.nodes[i];
        // for every
        if (node.camera != -1 && node_with_camera == -1) {
            node_with_camera = i;
        }
        matrix4d transform;
        if (node.translation.empty() && node.rotation.empty() && node.scale.empty()) {
            if (!node.matrix.empty()) {
                transform = matrix4d(node.matrix);
            } else {
                transform = matrix4d::eye();
            }
        } else {
            vector3f translation{0.f, 0.f, 0.f};
            vector4f rotation{0.f, 0.f, 0.f, 1.f};
            vector3f scale{1.f, 1.f, 1.f};
            if (!node.translation.empty()) {
                for (int j = 0; j < 3; ++j) {
                    translation[j] = (float) node.translation[j];
                }
            }
            if (!node.rotation.empty()) {
                for (int j = 0; j < 4; ++j) {
                    rotation[j] = (float) node.rotation[j];
                }
            }
            if (!node.scale.empty()) {
                for (int j = 0; j < 3; ++j) {
                    scale[j] = (float) node.scale[j];
                }
            }
            transform = matrix4d::TRS(translation, rotation, scale);
            if (!node.matrix.empty()) {
                transform = multiply(matrix4d(node.matrix), transform);
            }
        }

        node.matrix = transform;

        for (int child : node.children) {
            std::vector<double> &sub_transform = model.nodes[child].matrix;
            matrix4d new_sub_transform = transform;
            if (!sub_transform.empty()) {
                new_sub_transform = multiply(transform, matrix4d(sub_transform));
            }
            sub_transform = new_sub_transform;
        }
    }

    if (node_with_camera == -1)
        throw std::runtime_error("[gltf check] No camera found");

    {
        vector2i cam_canvas{width, height};
        vector3f cam_position = multiply(matrix4d(model.nodes[node_with_camera].matrix), vector3f{0.f, 0.f, 0.f});
        vector3f cam_axes[3];
        for (int i = 0; i < 3; ++i) {
            vector4f direction{};
            direction[i] = 1.f;
            direction.z = -direction.z;
            vector4f res = multiply(matrix4d(model.nodes[node_with_camera].matrix), direction);
            if (!model.nodes[node_with_camera].scale.empty()) {
                for (int j = 0; j < 3; ++j)
                    res[j] /= (float)model.nodes[node_with_camera].scale[j];
            }
            cam_axes[i] = normal(vector3f{res.x, res.y, res.z});
        }
        auto perspective = model.cameras[model.nodes[node_with_camera].camera].perspective;
        if (perspective.zfar > 0)
            scene.max_distance = (float) perspective.zfar;
        auto cam_fov_y = (float) perspective.yfov;
        float aspect_ratio = (perspective.aspectRatio != 0) ? (float) perspective.aspectRatio : (float) width /
                                                                                                (float) height;
        float cam_fov_x = 2.f * std::atan(std::tan(cam_fov_y * 0.5f) * aspect_ratio);

        scene.camera = std::make_unique<Camera>(cam_canvas, cam_position, cam_axes, cam_fov_x, cam_fov_y);
    }
    scene.samples = samples;

    for (auto &node : model.nodes) {
        if (node.mesh == -1)
            continue;

        tinygltf::Mesh &mesh = model.meshes[node.mesh];
        matrix4d transform = node.matrix;

        for (auto &p : mesh.primitives) {
            Material material;

            if (p.material != -1) {
                tinygltf::Material mesh_material = model.materials[p.material];
                for (int j = 0; j < 3; ++j)
                    material.color[j] = (float)mesh_material.pbrMetallicRoughness.baseColorFactor[j];
                if (mesh_material.pbrMetallicRoughness.baseColorFactor[3] < 1.) {
                    material.type = Material::Type::Dielectric;
                    if (mesh_material.extensions.count("KHR_materials_ior")) {
                        material.ior = (float)mesh_material.extensions["KHR_materials_ior"].Get("ior").GetNumberAsDouble();
                    } else {
                        material.ior = 1.5;
                    }
                }
                if (mesh_material.pbrMetallicRoughness.metallicFactor > 0.) {
                    material.type = Material::Type::Metallic;
                }
                for (int j = 0; j < 3; ++j)
                    material.emission[j] = (float)mesh_material.emissiveFactor[j];

                if (mesh_material.extensions.count("KHR_materials_emissive_strength")) {
                    float emission_strength = (float) mesh_material.extensions["KHR_materials_emissive_strength"].Get(
                            "emissiveStrength").GetNumberAsDouble();
                    material.emission *= emission_strength;
                }
            } else {
                material.type = Material::Type::Diffuse;
                material.color = {1.f, 1.f, 1.f};
            }

            auto indices_accessor = model.accessors[p.indices];
            auto position_accessor = model.accessors[p.attributes["POSITION"]];

            if (indices_accessor.type != TINYGLTF_TYPE_SCALAR)
                throw std::runtime_error("Index type not supported");
            if (position_accessor.type != TINYGLTF_TYPE_VEC3)
                throw std::runtime_error("Position type not supported");
            if (position_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
                throw std::runtime_error("Position component type not supported");

            tinygltf::Buffer &indices_buffer = model.buffers[model.bufferViews[indices_accessor.bufferView].buffer];
            tinygltf::Buffer &position_buffer = model.buffers[model.bufferViews[position_accessor.bufferView].buffer];

            size_t index_start = model.bufferViews[indices_accessor.bufferView].byteOffset + indices_accessor.byteOffset;
            size_t position_start = model.bufferViews[position_accessor.bufferView].byteOffset + position_accessor.byteOffset;

            unsigned char *indices_view = &indices_buffer.data[index_start];
            unsigned char *position_view = &position_buffer.data[position_start];

            for (size_t i = 0; i < indices_accessor.count / 3; ++i) {
                scene.objects.emplace_back(new Primitive());
                Primitive &primitive = *scene.objects.back();
                primitive.material = material;

                size_t index[3];
                switch (indices_accessor.componentType) {
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
                        auto *indices_view_16u = (uint16_t *)indices_view;
                        index[0] = indices_view_16u[i * 3];
                        index[1] = indices_view_16u[i * 3 + 1];
                        index[2] = indices_view_16u[i * 3 + 2];
                        break;
                    }
                    case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: {
                        auto *indices_view_32u = (uint32_t *)indices_view;
                        index[0] = indices_view_32u[i * 3];
                        index[1] = indices_view_32u[i * 3 + 1];
                        index[2] = indices_view_32u[i * 3 + 2];
                        break;
                    }
                    default:
                        throw std::runtime_error("Index component type not supported");
                }
                auto *position_view_32f = (float_t *)position_view;
                for (int v = 0; v < 3; ++v) {
                    for (int j = 0; j < 3; ++j) {
                        primitive.position[v][j] = position_view_32f[index[v] * 3 + j];
                    }
                    primitive.position[v] = multiply(matrix4d(node.matrix), primitive.position[v]);
                }
                primitive.position[1] = primitive.position[1] - primitive.position[0];
                primitive.position[2] = primitive.position[2] - primitive.position[0];
            }
        }
    }

    return {scene.camera, scene.objects, scene.bg_color, scene.ray_depth, scene.samples, scene.max_distance};
}
