#include "vector.h"

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

float length(vector2f v) {
    float mod = 0.f;
    for (int i = 0; i < 2; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

float length(vector3f v) {
    float mod = 0.f;
    for (int i = 0; i < 3; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

float length(vector4f v) {
    float mod = 0.f;
    for (int i = 0; i < 4; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

vector2f normal(vector2f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod};
}

vector3f normal(vector3f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod};
}

vector4f normal(vector4f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod, v.w * mod};
}

void normalize(vector2f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 2; ++i) {
        v[i] *= mod;
    }
}

void normalize(vector3f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 3; ++i) {
        v[i] *= mod;
    }
}

void normalize(vector4f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 4; ++i) {
        v[i] *= mod;
    }
}

uint8_t normal_to_ch8bit(float val) {
    return (uint8_t)std::round(clamp(val * 255, 0.f, 255.f));
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

vector3f &vector3f::operator+=(vector3f b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] += b[i];
    return *this;
}

vector3f operator+(vector3f a, vector3f b) {
    vector3f r = a;
    r += b;
    return r;
}

vector3f &vector3f::operator*=(vector3f b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] *= b[i];
    return *this;
}

vector3f &vector3f::operator*=(float b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] *= b;
    return *this;
}

vector3f operator-(vector3f a, vector3f b) {
    return a + (-b);
}

vector3f operator*(vector3f v, float t) {
    vector3f r = v;
    r *= t;
    return r;
}

vector3f operator*(float t, vector3f v) {
    return v * t;
}

vector3f operator*(vector3f a, vector3f b) {
    vector3f r = a;
    r *= b;
    return r;
}

vector3f operator/(vector3f a, vector3f b) {
    vector3f result{};
    for (int i = 0; i < 3; ++i)
        result[i] = a[i] / b[i];
    return result;
}

vector3f pow(vector3f v, float power) {
    return {std::pow(v.x, power), std::pow(v.y, power), std::pow(v.z, power)};
}

vector3f operator-(vector3f v) {
    return {-v.x, -v.y, -v.z};
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

vector3f rotate(vector3f v, vector3f axis, float angle) {
    float sin = std::sin(angle) / 2;
    axis = axis * sin;
    return rotate(v, {axis.x, axis.y, axis.z, std::cos(angle / 2)});
}

vector3f operator+(vector3f a, float t) {
    return {a.x + t, a.y + t, a.z + t};
}

template<class T>
T clamp(T v, T min, T max) {
    return std::max(std::min(v, max), min);
}

template<>
vector3f clamp<vector3f> (vector3f v, vector3f min, vector3f max) {
    vector3f res{};
    for (int i = 0; i < 3; ++i) {
        res[i] = clamp(v[i], min[i], max[i]);
    }
    return res;
}

vector3f saturate(const vector3f &color) {
    return clamp(color, {0.f, 0.f, 0.f}, {1.f, 1.f, 1.f});
}

vector3f aces_tonemap(const vector3f &x) {
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    return saturate((x*(a*x+b))/(x*(c*x+d)+e));
}
