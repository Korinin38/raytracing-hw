#include "canvas.h"

#include <stdexcept>
#include <fstream>

Canvas::Canvas(vector2i size) : size_(size) {
    data_.resize(size_.x * size_.y);
}

void Canvas::set(vector2i pos, vector3si color) {
    if (pos.x < 0 || pos.x >= size_.x
     || pos.y < 0 || pos.y >= size_.y) {
        throw std::runtime_error("Pixel out of canvas");
    }
    ptrdiff_t offset = size_.x * pos.y + pos.x;
    data_[offset] = color;
}

vector3si Canvas::get(vector2i pos) {
    if (pos.x < 0 || pos.x >= size_.x
        || pos.y < 0 || pos.y >= size_.y) {
        throw std::runtime_error("Pixel out of canvas");
    }
    ptrdiff_t offset = size_.x * pos.y + pos.x;
    return data_[offset];
}

vector2i Canvas::size() const {
    return size_;
}

int Canvas::width() const {
    return size_[0];
}

int Canvas::height() const {
    return size_[1];
}

void Canvas::write_to(const std::string &filename) const {
    std::ofstream out(filename, std::ios_base::binary);

    if (!out)
        throw std::runtime_error("File open error");

    out << "P6\n"
        << width() << " " << height() << std::endl
        << "255" << std::endl;

    char *data_char = (char *)(void *)data_.data();
    out.write(data_char, data_.size() * sizeof(vector3si));
//    for (auto i : data_) {
//        const unsigned char a[3] = {(unsigned char)i.x, (unsigned char)i.y, (unsigned char)i.z};
//        out.write(a, 1);
//    }
    out.close();
}
