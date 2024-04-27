#pragma once

#include <cmath>
#include <limits>
#include <string>

#pragma pack(push, 1)
struct vector2i {
    int x;
    int y;

    int &operator[](int i) { return *(&x + i); }
    const int &operator[](int i) const { return *(&x + i); }
};

struct vector3si {
    uint8_t x;
    uint8_t y;
    uint8_t z;

    uint8_t &operator[](int i) { return *(&x + i); }
    const uint8_t &operator[](int i) const { return *(&x + i); }
};

struct vector2f {
    float x;
    float y;

    float &operator[](int i) { return *(&x + i); }
    const float &operator[](int i) const { return *(&x + i); }
};

struct vector3f {
    float x;
    float y;
    float z;

    float &operator[](int i) { return *(&x + i); }
    const float &operator[](int i) const { return *(&x + i); }

    vector3f add(float t) const;
    vector3f &operator+=(vector3f b);
    vector3f &operator*=(vector3f b);
    vector3f &operator*=(float b);

    [[nodiscard]] bool is_zero() const {return x == 0 && y == 0 && z == 0; };
};

struct vector4f {
    float x;
    float y;
    float z;
    float w;

    float &operator[](int i) { return *(&x + i); }
    const float &operator[](int i) const { return *(&x + i); }
};
#pragma pack(pop)

inline int int_from_string(const std::string &s, size_t pos, size_t *idx = nullptr) {
    size_t end_pos = s.find(' ', pos);
    std::string raw_str = s.substr(pos, end_pos - pos);

    return std::stoi(raw_str, idx);
}

inline float float_from_string(const std::string &s, size_t pos, size_t *idx = nullptr) {
    size_t end_pos = s.find(' ', pos);
    std::string raw_str = s.substr(pos, end_pos - pos);

    double res = std::stod(raw_str, idx);
    return (float)res;
}

// extracts ints from string of format "X Y"
inline vector2i vec2i_from_string(const std::string &s, size_t pos = 0) {
    vector2i vec{};
    for (int i = 0; i < 2; ++i) {
        size_t new_pos;
        vec[i] = int_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

// extracts ints from string of format "X Y Z"
inline vector3si vec3si_from_string(const std::string &s, size_t pos = 0) {
    vector3si vec{};
    for (int i = 0; i < 3; ++i) {
        size_t new_pos;
        vec[i] = int_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

// extracts floats from string of format "X.XXX Y.YYY"
inline vector2f vec2f_from_string(const std::string &s, size_t pos = 0) {
    vector2f vec{};
    for (int i = 0; i < 2; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ"
inline vector3f vec3f_from_string(const std::string &s, size_t pos = 0) {
    vector3f vec{};
    for (int i = 0; i < 3; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
        pos += new_pos + 1;
    }
    return vec;
}

// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ W.WWW"
inline vector4f vec4f_from_string(const std::string &s, size_t pos = 0) {
    vector4f vec{};
    for (int i = 0; i < 4; ++i) {
        size_t new_pos;
        vec[i] = float_from_string(s, pos, &new_pos);
        pos += new_pos + 1;

    }
    return vec;
}

inline float length(vector2f v) {
    float mod = 0.f;
    for (int i = 0; i < 2; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

inline float length(vector3f v) {
    float mod = 0.f;
    for (int i = 0; i < 3; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

inline float length(vector4f v) {
    float mod = 0.f;
    for (int i = 0; i < 4; ++i) {
        mod += v[i] * v[i];
    }
    return std::sqrt(mod);
}

[[nodiscard]]
inline vector2f normal(vector2f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod};
}

[[nodiscard]]
inline vector3f normal(vector3f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod};
}

[[nodiscard]]
inline vector4f normal(vector4f v) {
    float mod = length(v);
    mod = 1.f / mod;
    return {v.x * mod, v.y * mod, v.z * mod, v.w * mod};
}

inline void normalize(vector2f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 2; ++i) {
        v[i] *= mod;
    }
}

inline void normalize(vector3f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 3; ++i) {
        v[i] *= mod;
    }
}

inline void normalize(vector4f &v) {
    float mod = length(v);
    mod = 1.f / mod;
    for (int i = 0; i < 4; ++i) {
        v[i] *= mod;
    }
}

template<class T>
inline T clamp(T v, T min, T max) {
    return std::max(std::min(v, max), min);
}

template<>
inline vector3f clamp<vector3f>(vector3f v, vector3f min, vector3f max) {
    vector3f res{};
    for (int i = 0; i < 3; ++i) {
        res[i] = clamp(v[i], min[i], max[i]);
    }
    return res;
}

// translation from [0.0, 1.0] to [0, 255]
inline uint8_t normal_to_ch8bit(float val) {
    return (uint8_t) std::round(clamp(val * 255, 0.f, 255.f));
}

// translation from [0.0, 1.0] to [0, 255]
inline vector3si normal_to_ch8bit(vector3f val) {
    vector3si result{};
    for (int i = 0; i < 3; ++i) {
        result[i] = normal_to_ch8bit(val[i]);
    }
    return result;
}

// translation from [0, 255] to [0.0, 1.0]
inline float ch8bit_to_normal(int val) {
    return (float) val / 255.f;
}

// translation from [0, 255] to [0.0, 1.0]
inline vector3f ch8bit_to_normal(vector3si val) {
    vector3f result{};
    for (int i = 0; i < 3; ++i) {
        result[i] = ch8bit_to_normal(val[i]);
    }
    return result;
}

inline vector3f operator-(vector3f v) {
    return {-v.x, -v.y, -v.z};
}

// conjugate
inline vector4f operator*(vector4f q) {
    return {-q.x, -q.y, -q.z, q.w};
}

inline vector3f &vector3f::operator+=(vector3f b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] += b[i];
    return *this;
}

inline vector3f operator+(vector3f a, vector3f b) {
    vector3f r = a;
    r += b;
    return r;
}

inline vector3f &vector3f::operator*=(vector3f b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] *= b[i];
    return *this;
}

inline vector3f &vector3f::operator*=(float b) {
    for (int i = 0; i < 3; ++i)
        (*this)[i] *= b;
    return *this;
}

inline vector3f vector3f::add(float t) const {
    return {x + t, y + t, z + t};
}

inline vector3f operator-(vector3f a, vector3f b) {
    return a + (-b);
}

inline vector3f operator*(vector3f v, float t) {
    vector3f r = v;
    r *= t;
    return r;
}

inline vector3f operator*(float t, vector3f v) {
    return v * t;
}

inline vector3f operator*(vector3f a, vector3f b) {
    vector3f r = a;
    r *= b;
    return r;
}

inline vector3f operator/(vector3f a, vector3f b) {
    vector3f result{};
    for (int i = 0; i < 3; ++i)
        result[i] = a[i] / b[i];
    return result;
}

inline bool operator==(vector3f a, vector3f b) {
    return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}

inline bool operator!=(vector3f a, vector3f b) {
    return !(a == b);
}

inline vector3f pow(const vector3f v, float power) {
    return {std::pow(v.x, power), std::pow(v.y, power), std::pow(v.z, power)};
}

inline float dot(const vector3f a, const vector3f b) {
    float s = 0.f;
    for (int i = 0; i < 3; ++i)
        s += a[i] * b[i];
    return s;
}

inline vector3f cross(const vector3f a, const vector3f b) {
    return {
            a[1] * b[2] - b[1] * a[2],
            a[2] * b[0] - b[2] * a[0],
            a[0] * b[1] - b[0] * a[1]
    };
}

inline vector3f rotate(const vector3f v, const vector4f q) {
    vector3f u{q.x, q.y, q.z};
    float s = q.w;

    return 2.f * dot(u, v) * u
           + (s * s - dot(u, u)) * v
           + 2.f * s * cross(u, v);
}

inline vector3f rotate(vector3f v, vector3f axis, float angle) {
    float sin = std::sin(angle) / 2;
    axis = axis * sin;
    return rotate(v, {axis.x, axis.y, axis.z, std::cos(angle / 2)});
}


inline vector3f saturate(const vector3f &color) {
    return clamp(color, {0.f, 0.f, 0.f}, {1.f, 1.f, 1.f});
}

inline vector3f aces_tonemap(const vector3f &x) {
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    return saturate((x * ((a * x).add(b))) / (x * ((c * x).add(d))).add(e));
}


inline vector3f get_max_vec3f() {
    const float m = std::numeric_limits<float>::max();
    return {m, m, m};
}

inline vector3f get_min_vec3f() {
    const float m = -std::numeric_limits<float>::max();
    return {m, m, m};
}

inline vector3f min(vector3f a, vector3f b) {
    return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}

inline vector3f max(vector3f a, vector3f b) {
    return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}

inline vector3f lerp(vector3f val1, vector3f val2, float coeff) {
    return val1 * (1 - coeff) + val2 * coeff;
}

// see https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#specular-brdf
inline float ggx_microfacet_distribution(float alpha, vector3f N, vector3f H) {
    float alpha2 = alpha * alpha;
    float NdotH = dot(N, H);
    if (NdotH <= 0)
        return 0;
    return alpha2 * M_1_PIf32 / ((NdotH * NdotH * (alpha2 - 1) + 1) * (NdotH * NdotH * (alpha2 - 1) + 1));
}

// see https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#specular-brdf
inline float smith_joint_masking_shadowing_function(float alpha, vector3f N, vector3f H, vector3f V, vector3f L) {
    float alpha2 = alpha * alpha;
    if (dot(N, H) <= 0 || dot(V, H) <= 0)
        return 0;
    float Ndot[2] = {std::abs(dot(N, L)), std::abs(dot(N, V))};
    float res = 1.f;
    for (auto ndot : Ndot) {
        res *= 2 * ndot / (ndot + std::sqrt(alpha2 + (1 - alpha2) * ndot * ndot));
    }
    return res;
}

inline vector4f quat_from_two_vectors(vector3f u, vector3f v) {
//    vector3f w = normal(cross(u, v));
    vector3f w = cross(u, v);
    return normal(vector4f{w.x, w.y, w.z, 1.f + dot(u, v)});
}