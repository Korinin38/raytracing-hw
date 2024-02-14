#include "base.h"

#include <string>
#include <cmath>

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

vector3si vec3si_from_string(const std::string &s, size_t pos) {
    vector3si vec{};
    for (int i = 0; i < 3; ++i) {
        size_t new_pos;
        vec[i] = int_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

vector2f vec2f_from_string(const std::string &s, size_t pos) {
    vector2f vec{};
    for (int i = 0; i < 2; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
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

vector2f normal(vector2f v) {
    float mod = 0.f;
    for (int i = 0; i < 2; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    return {v.x * mod, v.y * mod};
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

void normalize(vector2f &v) {
    float mod = 0.f;
    for (int i = 0; i < 2; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    for (int i = 0; i < 2; ++i) {
        v[i] *= mod;
    }
}

void normalize(vector3f &v) {
    float mod = 0.f;
    for (int i = 0; i < 3; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    for (int i = 0; i < 3; ++i) {
        v[i] *= mod;
    }
}

void normalize(vector4f &v) {
    float mod = 0.f;
    for (int i = 0; i < 4; ++i) {
        mod += v[i] * v[i];
    }

    mod = 1.f / mod;
    for (int i = 0; i < 4; ++i) {
        v[i] *= mod;
    }
}

uint8_t normal_to_ch8bit(float val) {
    return std::floor(val * 255);
}

vector3si normal_to_ch8bit(vector3f val) {
    vector3si result{};
    for (int i = 0; i < 3; ++i) {
        result[i] = normal_to_ch8bit(val[i]);
    }
    return result;
}

float ch8bit_to_normal(int val) {
    return (float)val / 255.f;
}

vector3f ch8bit_to_normal(vector3si val) {
    vector3f result{};
    for (int i = 0; i < 3; ++i) {
        result[i] = ch8bit_to_normal(val[i]);
    }
    return result;
}

vector3f operator+(vector3f a, vector3f b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

vector3f operator*(vector3f v, float t) {
    return {v.x * t, v.y * t, v.z * t};
}

vector3f operator*(float t, vector3f v) {
    return v * t;
}

vector3f operator/(vector3f a, vector3f b) {
    vector3f result{};
    for (int i = 0; i < 3; ++i)
        result[i] = a[i] / b[i];
    return result;
}

vector4f operator*(vector4f q) {
    return {-q.x, -q.y, -q.z, q.w};
}

float dot(vector3f a, vector3f b) {
    float s = 0.f;
    for (int i = 0; i < 3; ++i)
        s += a[i] * b[i];
    return s;
}

vector3f cross(vector3f a, vector3f b) {
    return {
        a[1] * b[2] - b[1] * a[2],
        a[2] * b[0] - b[2] * a[0],
        a[0] * b[1] - b[0] * a[1]
    };
}

vector3f rotate(vector3f v, vector4f q) {
    vector3f u{q.x, q.y, q.z};
    float s = q.w;

    return 2.f * dot(u, v) * u
         + (s * s - dot(u, u)) * v
         + 2.f * s * cross(u, v);
}
