#include "camera.h"

//Camera::Camera()
Camera::Camera(vector2i size, vector3f position, vector3f *axes, float fov_x)
        : canvas_size_(size),
          position_(position),
          axes_{axes[0], axes[1], axes[2]},
          fov_{fov_x, fov_x / (float)size.x * (float)size.y}
          {}
