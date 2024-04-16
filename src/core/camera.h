#pragma once

#include "utils/base.h"
#include "render/canvas.h"
#include "geometry/primitive.h"

#include <memory>

class Camera;

typedef std::unique_ptr<Camera> camera_uniq_ptr;

class Camera {
public:
    Camera(vector2i size, vector3f position, const vector3f (&axes)[3], float fov_x, float fov_y);
    Camera(vector2i size, vector3f position, const vector3f (&axes)[3], float fov_x);

    Canvas canvas;

    vector2i get_canvas_size() const;

    vector3f get_position() const;
    void set_position(vector3f position);

    vector3f get_axis(int dim) const;
    vector3f get_axis_right() const;
    vector3f get_axis_up() const;
    vector3f get_axis_forward() const;

    vector2f get_fov() const;

    // NB: casts into a center of a pixel
    Ray cast_in_pixel(vector2i p, vector2f rand_offset = {});

private:
    vector3f position_{};

    // 0 - right
    // 1 - up
    // 2 - forward
    vector3f axes_[3]{};

    // in radians
    vector2f fov_{};
};
