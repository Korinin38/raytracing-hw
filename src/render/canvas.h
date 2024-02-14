#pragma once

#include "utils/base.h"

#include "vector"

class Canvas {
public:
    explicit Canvas(vector2i size);

    vector2i size() const;
    int width() const;
    int height() const;

    void set(vector2i pos, vector3si color);
    vector3si get(vector2i pos);

    const uint8_t max_color_value_ = 255;

    void write_to(const std::string &filename) const;
private:
    vector2i size_;
    std::vector<vector3si> data_;
};
