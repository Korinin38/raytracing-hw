#include "base.h"
#include <string>

int int_from_string(const std::string &s, size_t pos, size_t *idx) {
    size_t end_pos = s.find(' ', pos);
    std::string raw_str = s.substr(pos, end_pos - pos);

    return std::stoi(raw_str, idx);
}

float float_from_string(const std::string &s, size_t pos, size_t *idx) {
    size_t end_pos = s.find(' ', pos);
    std::string raw_str = s.substr(pos, end_pos - pos);

    return std::stof(raw_str, idx);
}

vector2i vec2i_from_string(const std::string &s, size_t pos) {
    vector2i vec{};
    for (int i = 0; i < 2; ++i) {
        size_t new_pos;
        vec[i] = int_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

vector3f vec3f_from_string(const std::string &s, size_t pos) {
    vector3f vec{};
    for (int i = 0; i < 3; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

vector4f vec4f_from_string(const std::string &s, size_t pos) {
    vector4f vec{};
    for (int i = 0; i < 4; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
        pos += new_pos + 1;

    }
    return vec;
}

vector3f normal(vector3f v) {
    float mod = 0.f;
    for (int i = 0; i < 3; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod};
}

vector4f normal(vector4f v) {
    float mod = 0.f;
    for (int i = 0; i < 4; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod, v.w * mod};
}

void normal(vector3f &v) {
    float mod = 0.f;
    for (int i = 0; i < 3; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    for (int i = 0; i < 3; ++i) {
        v[i] *= mod;
    }
}

void normal(vector4f &v) {
    float mod = 0.f;
    for (int i = 0; i < 4; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    for (int i = 0; i < 4; ++i) {
        v[i] *= mod;
    }
}
