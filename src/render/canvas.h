#pragma once

#include <utils/vector.h>

#include <stdexcept>
#include <fstream>

class Canvas {
public:
    explicit Canvas(vector2i size) : size_(size) {
        data_.resize(size_.x * size_.y);
    }

    vector2i size() const { return size_; }
    int width() const { return size_[0]; }
    int height() const {return size_[1]; }

    void set(vector2i pos, vector3si color);
    vector3si get(vector2i pos) const;

    const uint8_t max_color_value_ = 255;

    void write_to(const std::string &filename) const;

    void smooth();
private:
    vector2i size_;
    std::vector<vector3si> data_;
};

inline void Canvas::set(vector2i pos, vector3si color) {
    if (pos.x < 0 || pos.x >= size_.x
        || pos.y < 0 || pos.y >= size_.y) {
        throw std::runtime_error("Pixel out of canvas");
    }
    ptrdiff_t offset = size_.x * pos.y + pos.x;
    data_[offset] = color;
}

inline vector3si Canvas::get(vector2i pos) const {
    if (pos.x < 0 || pos.x >= size_.x
        || pos.y < 0 || pos.y >= size_.y) {
        throw std::runtime_error("Pixel out of canvas");
    }
    ptrdiff_t offset = size_.x * pos.y + pos.x;
    return data_[offset];
}

inline void Canvas::smooth() {
    std::vector<vector3si> res(size_.x * size_.y);
    for (int i = 0; i < width(); ++i)
    for (int j = 0; j < height(); ++j) {
        int weight = 0;
        vector3f value{};

        for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj) {
            if (i + di < 0 || i + di >= width()
             || j + dj < 0 || j + dj >= height())
                continue;
            int w = 1 * (1 + std::abs(di)) * (1 + std::abs(dj));
            weight += w;

            ptrdiff_t offset = size_.x * (j + dj) + i + di;

            value += (ch8bit_to_normal(data_[offset])) * (float)w;
        }

        ptrdiff_t offset = size_.x * j + i;
        res[offset] = normal_to_ch8bit(value * (1 / (float)weight));
    }

    std::swap(data_, res);
}

inline void Canvas::write_to(const std::string &filename) const {
    std::ofstream out(filename, std::ios_base::binary);

    if (!out)
        throw std::runtime_error("File open error");

    out << "P6\n"
        << width() << " " << height() << std::endl
        << "255" << std::endl;

    char *data_char = (char *)(void *)data_.data();
    out.write(data_char, data_.size() * sizeof(vector3si));
    out.close();
}

