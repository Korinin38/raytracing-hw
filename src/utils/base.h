#pragma once

#include <string>

struct vector2i {
    int x;
    int y;

    int &operator[](int i) { return *(&x + i); }
};

struct vector3f {
    float x;
    float y;
    float z;

    float &operator[](int i) { return *(&x + i); }
};

struct vector4f {
    float x;
    float y;
    float z;
    float w;

    float &operator[](int i) { return *(&x + i); }
};

extern int int_from_string(const std::string &s, size_t pos, size_t *idx = nullptr);
extern float float_from_string(const std::string &s, size_t pos, size_t *idx = nullptr);
// extracts ints from string of format "X Y"
extern vector2i vec2i_from_string(const std::string &s, size_t pos = 0);
// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ"
extern vector3f vec3f_from_string(const std::string &s, size_t pos = 0);
// extracts floats from string of format "X.XXX Y.YYY Z.ZZZ W.WWW"
extern vector4f vec4f_from_string(const std::string &s, size_t pos = 0);

[[nodiscard]] vector3f normal(vector3f v);
[[nodiscard]] vector4f normal(vector4f v);

void normal(vector3f &v);
void normal(vector4f &v);