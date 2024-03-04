#pragma once

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

    vector3f &operator+=(vector3f b);
    vector3f &operator*=(vector3f b);
    vector3f &operator*=(float b);
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

extern int int_from_string(const std::string &s, size_t pos, size_t *idx = nullptr);
extern float float_from_string(const std::string &s, size_t pos, size_t *idx = nullptr);
// extracts ints from string of format "X Y"
extern vector2i vec2i_from_string(const std::string &s, size_t pos = 0);
// extracts ints from string of format "X Y Z"
extern vector3si vec3si_from_string(const std::string &s, size_t pos = 0);
// extracts floats from string of format "X.XXX Y.YYY"
extern vector2f vec2f_from_string(const std::string &s, size_t pos = 0);
// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ"
extern vector3f vec3f_from_string(const std::string &s, size_t pos = 0);
// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ W.WWW"
extern vector4f vec4f_from_string(const std::string &s, size_t pos = 0);

float length(vector2f v);
float length(vector3f v);
float length(vector4f v);

[[nodiscard]] vector2f normal(vector2f v);
[[nodiscard]] vector3f normal(vector3f v);
[[nodiscard]] vector4f normal(vector4f v);

void normalize(vector2f &v);
void normalize(vector3f &v);
void normalize(vector4f &v);

// translation from [0.0, 1.0] to [0, 255]
uint8_t normal_to_ch8bit(float val);
vector3si normal_to_ch8bit(vector3f val);

// translation from [0, 255] to [0.0, 1.0]
float ch8bit_to_normal(int val);
vector3f ch8bit_to_normal(vector3si val);

float dot(vector3f a, vector3f b);
vector3f cross(vector3f a, vector3f b);

vector3f operator+(vector3f a, float t);
vector3f operator+(vector3f a, vector3f b);
vector3f operator-(vector3f a, vector3f b);
vector3f operator*(vector3f v, float t);
vector3f operator*(float t, vector3f v);
vector3f operator*(vector3f a, vector3f b);
vector3f operator/(vector3f a, vector3f b);

vector3f pow(vector3f v, float power);

vector3f operator-(vector3f v);
// conjugate
vector4f operator*(vector4f q);

vector3f rotate(vector3f v, vector4f q);

vector3f rotate(vector3f v, vector3f axis, float angle);

template<class T>
T clamp(T v, T min, T max);

vector3f saturate(vector3f const & color);
vector3f aces_tonemap(vector3f const & x);
