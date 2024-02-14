#include "camera.h"

#include <cmath>

Camera::Camera(vector2i size, vector3f position, vector3f *axes, float fov_x)
        : canvas_(size),
          position_(position),
          axes_{axes[0], axes[1], axes[2]} {
    fov_.x = fov_x;
    float ratio_inv =  (float)size.y / (float)size.x;
    fov_.y = 2.f * std::atan(std::tan(fov_.x * 0.5f) * ratio_inv);
}

vector2i Camera::get_canvas_size() const {
    return canvas_.size();
}

vector3f Camera::get_position() const {
    return position_;
}

vector3f Camera::set_position(vector3f position) {
    position_ = position;
}

vector3f Camera::get_axis(int dim) const {
    return axes_[dim];
}

vector3f Camera::get_axis_right() const {
    get_axis(0);
}

vector3f Camera::get_axis_up() const {
    get_axis(1);
}

vector3f Camera::get_axis_forward() const {
    get_axis(2);
}

vector2f Camera::get_fov() const {
    return fov_;
}

Ray Camera::cast_in_pixel(vector2i p) {
    vector3f t{};
    t.x = (2.f * ((float)p.x + 0.5f) / (float)canvas_.width() - 1) * std::tan(fov_.x / 2);
    t.y = -(2.f * ((float)p.y + 0.5f) / (float)canvas_.height() - 1) * std::tan(fov_.y / 2);
    t.z = 1;

    vector3f d{};
    for (int i = 0; i < 3; ++i) {
        d = d + t[i] * axes_[i];
    }

    return {position_, d};
}
