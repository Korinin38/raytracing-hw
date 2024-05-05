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
    int ray_depth = 6;
    int samples = 256;
    std::vector<Mesh> meshes;
    std::vector<Texture> textures;
    std::vector<Primitive> objects;
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

    // camera
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
            cam_axes[i] = normal(res.reduce());
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

    // textures
    for (auto &gltf_image : model.images) {
        scene.textures.emplace_back();
        Texture &t = scene.textures.back();

        t.data = std::move(gltf_image.image);
        t.width = gltf_image.width;
        t.height = gltf_image.height;
        t.channels = gltf_image.component;
        if (gltf_image.bits != 8)
            throw std::runtime_error("Texture format not supported: only 8 bit channels are supported");
        t.bytes_per_channel = gltf_image.bits / 8;
    }

    // meshes (gltf_mesh != mesh; gltf_primitive == mesh. confusing, but simplifies material handling)
    for (auto &node : model.nodes) {
        if (node.mesh == -1)
            continue;


        tinygltf::Mesh &gltf_mesh = model.meshes[node.mesh];

        for (auto &p : gltf_mesh.primitives) {
            scene.meshes.emplace_back();
            Mesh &mesh = scene.meshes.back();

            matrix4d transform = node.matrix;
            mesh.normal_transform = inverse(transpose(transform));
            matrix4d normal_transform = mesh.normal_transform;

            // material
            {
                Material &material = mesh.material;

                if (p.material != -1) {
                    tinygltf::Material gltf_material = model.materials[p.material];
                    for (int j = 0; j < 3; ++j)
                        material.base_color[j] = (float) gltf_material.pbrMetallicRoughness.baseColorFactor[j];
                    if (gltf_material.pbrMetallicRoughness.baseColorFactor[3] < 1.) {
                        material.alpha = (float) gltf_material.pbrMetallicRoughness.baseColorFactor[3];
                        if (gltf_material.extensions.count("KHR_materials_ior")) {
                            material.ior = (float) gltf_material.extensions["KHR_materials_ior"].Get(
                                    "ior").GetNumberAsDouble();
                        } else {
                            material.ior = 1.5;
                        }
                    }
                    material.metallic = (float) gltf_material.pbrMetallicRoughness.metallicFactor;
                    material.roughness2 = (float) gltf_material.pbrMetallicRoughness.roughnessFactor;
                    // make it squared
                    material.roughness2 *= material.roughness2;
                    for (int j = 0; j < 3; ++j)
                        material.emission[j] = (float) gltf_material.emissiveFactor[j];

                    if (gltf_material.extensions.count("KHR_materials_emissive_strength")) {
                        float emission_strength = (float) gltf_material.extensions["KHR_materials_emissive_strength"].Get(
                                "emissiveStrength").GetNumberAsDouble();
                        material.emission *= emission_strength;
                    }
                    auto get_image = [&](int i) {return (i == -1) ? -1 : model.textures[i].source; };

                    material.base_color_i = get_image(gltf_material.pbrMetallicRoughness.baseColorTexture.index);
                    material.normal_i = get_image(gltf_material.normalTexture.index);
                    material.metallic_roughness_i = get_image(gltf_material.pbrMetallicRoughness.metallicRoughnessTexture.index);
                    material.emission_i = get_image(gltf_material.emissiveTexture.index);
                }
            }

            unsigned char *indices_view;
            float_t *position_view;
            float_t *normal_view;
            float_t *texcoord_view;
            float_t *tangent_view;

            size_t indices_count;
            int indices_component_type;

            {
                auto indices_accessor = model.accessors[p.indices];
                if (indices_accessor.type != TINYGLTF_TYPE_SCALAR)
                    throw std::runtime_error("Index type not supported");
                tinygltf::Buffer &indices_buffer = model.buffers[model.bufferViews[indices_accessor.bufferView].buffer];
                size_t index_start =
                        model.bufferViews[indices_accessor.bufferView].byteOffset + indices_accessor.byteOffset;
                indices_view = &indices_buffer.data[index_start];


                auto position_accessor = model.accessors[p.attributes["POSITION"]];
                if (position_accessor.type != TINYGLTF_TYPE_VEC3)
                    throw std::runtime_error("Position type not supported");
                if (position_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
                    throw std::runtime_error("Position component type not supported");
                tinygltf::Buffer &position_buffer = model.buffers[model.bufferViews[position_accessor.bufferView].buffer];
                size_t position_start =
                        model.bufferViews[position_accessor.bufferView].byteOffset + position_accessor.byteOffset;
                position_view = (float_t *) (&position_buffer.data[position_start]);


                if (p.attributes.count("NORMAL") == 0) {
                    normal_view = nullptr;
                } else {
                    auto normal_accessor = model.accessors[p.attributes["NORMAL"]];

                    indices_count = indices_accessor.count;
                    indices_component_type = indices_accessor.componentType;

                    if (normal_accessor.type != TINYGLTF_TYPE_VEC3)
                        throw std::runtime_error("Normal type not supported");
                    if (normal_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
                        throw std::runtime_error("Normal component type not supported");

                    tinygltf::Buffer &normal_buffer = model.buffers[model.bufferViews[normal_accessor.bufferView].buffer];

                    size_t normal_start =
                            model.bufferViews[normal_accessor.bufferView].byteOffset + normal_accessor.byteOffset;

                    normal_view = (float_t *) (&normal_buffer.data[normal_start]);
                }

                if (p.attributes.count("TEXCOORD_0") == 0) {
                    texcoord_view = nullptr;
                } else {
                    auto texcoord_accessor = model.accessors[p.attributes["TEXCOORD_0"]];
                    if (texcoord_accessor.type != TINYGLTF_TYPE_VEC2)
                        throw std::runtime_error("Texcoord type not supported");
                    if (texcoord_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
                        throw std::runtime_error("Texcoord component type not supported");
                    tinygltf::Buffer &texcoord_buffer = model.buffers[model.bufferViews[texcoord_accessor.bufferView].buffer];
                    size_t texcoord_start =
                            model.bufferViews[texcoord_accessor.bufferView].byteOffset + texcoord_accessor.byteOffset;

                    texcoord_view = (float_t *) (&texcoord_buffer.data[texcoord_start]);
                }

                if (p.attributes.count("TANGENT") == 0) {
                    tangent_view = nullptr;
                } else {
                    auto tangent_accessor = model.accessors[p.attributes["TANGENT"]];
                    if (tangent_accessor.type != TINYGLTF_TYPE_VEC4)
                        throw std::runtime_error("Tangent type not supported");
                    if (tangent_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
                        throw std::runtime_error("Tangent component type not supported");
                    tinygltf::Buffer &tangent_buffer = model.buffers[model.bufferViews[tangent_accessor.bufferView].buffer];
                    size_t tangent_start =
                            model.bufferViews[tangent_accessor.bufferView].byteOffset + tangent_accessor.byteOffset;

                    tangent_view = (float_t *) (&tangent_buffer.data[tangent_start]);
                }
            }

            for (size_t i = 0; i < indices_count / 3; ++i) {
                scene.objects.emplace_back();
                Primitive &primitive = scene.objects.back();
                primitive.mesh_id = (int)scene.meshes.size() - 1;
                primitive.mesh = &scene.meshes.back();

                size_t index[3];
                switch (indices_component_type) {
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

                for (int v = 0; v < 3; ++v) {
                    for (int j = 0; j < 3; ++j) {
                        primitive.position[v][j] = position_view[index[v] * 3 + j];
                    }
                    primitive.position[v] = multiply(transform, primitive.position[v]);
                }
                primitive.position[1] = primitive.position[1] - primitive.position[0];
                primitive.position[2] = primitive.position[2] - primitive.position[0];

                for (int v = 0; v < 3; ++v) {
                    if (normal_view) {
                        for (int j = 0; j < 3; ++j) {
                            primitive.normal[v][j] = normal_view[index[v] * 3 + j];
                        }
                    } else {
                        primitive.normal[v] = {0, 0, 1};
                    }
                    primitive.normal[v] = ::normal(multiplyVector(normal_transform, primitive.normal[v]));
                }

                if (texcoord_view) {
                    for (int v = 0; v < 3; ++v) {
                        for (int j = 0; j < 2; ++j) {
                            primitive.texcoord[v][j] = texcoord_view[index[v] * 2 + j];
                        }
                    }
                }

                if (tangent_view) {
                    for (int v = 0; v < 3; ++v) {
                        for (int j = 0; j < 4; ++j) {
                            primitive.tangent[v][j] = tangent_view[index[v] * 4 + j];
                        }
                    }
                    float side = primitive.tangent[0].w;
                    if (primitive.tangent[1].w != side || primitive.tangent[2].w != side)
                        std::cout << "Warning: tangent.w values differ within triangle" << std::endl;
                }
            }
        }
    }

    return {scene.camera, scene.meshes, scene.objects, scene.textures, scene.ray_depth, scene.samples,
            scene.max_distance};
}
