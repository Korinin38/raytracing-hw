#pragma once

#include <utils/vector.h>

#include <vector>

template <typename T>
struct matrix4 {
    T data[16] = {};

    matrix4() = default;
    matrix4(const std::vector<T> &a) { // NOLINT(*-explicit-constructor)
        std::copy(a.begin(), a.end(), data);
    }

    T* operator[](int i) { return &data[4 * i]; } // get column
    const T* operator[](int i) const { return &data[4 * i]; } // get column

    static matrix4 eye() {
        return {{1.f, 0.f, 0.f, 0.f,
                0.f, 1.f, 0.f, 0.f,
                0.f, 0.f, 1.f, 0.f,
                0.f, 0.f, 0.f, 1.f}};
    }

    static matrix4 TRS(vector3f t, vector4f r, vector3f s) {
        return {{
                (1.0f-2.0f*(r.y*r.y+r.z*r.z))*s.x,  (r.x*r.y+r.z*r.w)*s.x*2.0f,         (r.x*r.z-r.y*r.w)*s.x*2.0f,         0.f,
                (r.x*r.y-r.z*r.w)*s.y*2.0f,         (1.0f-2.0f*(r.x*r.x+r.z*r.z))*s.y,  (r.y*r.z+r.x*r.w)*s.y*2.0f,         0.f,
                (r.x*r.z+r.y*r.w)*s.z*2.0f,         (r.y*r.z-r.x*r.w)*s.z*2.0f,         (1.0f-2.0f*(r.x*r.x+r.y*r.y))*s.z,  0.f,
                t.x,                                t.y,                                t.z,                                1.f
                }};
    }

    operator std::vector<T> () const {
        std::vector<T> res;
        res.insert(res.end(), &data[0], &data[16]);
        return res;
    }
};

template <typename T>
matrix4<T> multiply(const matrix4<T> &a, const matrix4<T> &b) {
    matrix4<T> dest;
    dest[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0] + a[0][3] * b[3][0];
    dest[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1] + a[0][3] * b[3][1];
    dest[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2] + a[0][3] * b[3][2];
    dest[0][3] = a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3] * b[3][3];
    dest[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0] + a[1][3] * b[3][0];
    dest[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1] + a[1][3] * b[3][1];
    dest[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2] + a[1][3] * b[3][2];
    dest[1][3] = a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3] * b[3][3];
    dest[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0] + a[2][3] * b[3][0];
    dest[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1] + a[2][3] * b[3][1];
    dest[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2] + a[2][3] * b[3][2];
    dest[2][3] = a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3] * b[3][3];
    dest[3][0] = a[3][0] * b[0][0] + a[3][1] * b[1][0] + a[3][2] * b[2][0] + a[3][3] * b[3][0];
    dest[3][1] = a[3][0] * b[0][1] + a[3][1] * b[1][1] + a[3][2] * b[2][1] + a[3][3] * b[3][1];
    dest[3][2] = a[3][0] * b[0][2] + a[3][1] * b[1][2] + a[3][2] * b[2][2] + a[3][3] * b[3][2];
    dest[3][3] = a[3][0] * b[0][3] + a[3][1] * b[1][3] + a[3][2] * b[2][3] + a[3][3] * b[3][3];
    return dest;
}

template <typename T>
vector4f multiply(const matrix4<T> &mat, const vector4f &t) {
    vector4f res{};
    for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
        res[i] += mat[j][i] * t[j];
    }
    return res;
}

template <typename T>
vector3f multiply(const matrix4<T> &mat, const vector3f &t) {
    vector4f res = multiply(mat, {t.x, t.y, t.z, 1.f});
    return {res.x, res.y, res.z};
}

typedef matrix4<float> matrix4f;
typedef matrix4<double> matrix4d;

