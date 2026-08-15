#include "boundary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

double clamp(double value, double minimum, double maximum) {
    return std::max(minimum, std::min(value, maximum));
}

}

RectangularPatch makeRectangularPatch(const Vec3& center,
                                      const Vec3& normal,
                                      const Vec3& width_direction,
                                      double width, double height) {
    if (width <= 0.0 || height <= 0.0) {
        throw std::invalid_argument("Boundary dimensions must be positive");
    }

    const double normal_length = normal.length();
    if (normal_length == 0.0) {
        throw std::invalid_argument("Boundary normal must be non-zero");
    }
    const Vec3 unit_normal = normal * (1.0 / normal_length);

    const Vec3 tangent = width_direction
                       - unit_normal * width_direction.dot(unit_normal);
    const double tangent_length = tangent.length();
    if (tangent_length == 0.0) {
        throw std::invalid_argument("Width direction must not be parallel to normal");
    }

    return {{center, unit_normal}, tangent * (1.0 / tangent_length), width, height};
}

SurfaceQuery closestPoint(const RectangularPatch& patch, const Vec3& point) {
    const Vec3 height_axis = patch.heightAxis();
    const Vec3 relative = point - patch.plane.point;

    const double width_coordinate = clamp(relative.dot(patch.width_axis),
                                          -0.5 * patch.width, 0.5 * patch.width);
    const double height_coordinate = clamp(relative.dot(height_axis),
                                           -0.5 * patch.height, 0.5 * patch.height);

    const Vec3 closest = patch.plane.point
                       + patch.width_axis * width_coordinate
                       + height_axis * height_coordinate;
    const Vec3 separation = point - closest;
    const double distance_squared = separation.dot(separation);

    const Vec3 normal = distance_squared > 0.0
                      ? separation * (1.0 / std::sqrt(distance_squared))
                      : patch.plane.normal;

    return {closest, normal, distance_squared};
}

void applyBoundaryContact(Particle& particle, const Boundary& boundary,
                          const SurfaceQuery& query) {
    const double radius_squared = particle.radius * particle.radius;
    if (query.distance_squared >= radius_squared) {
        return;
    }

    const double overlap = particle.radius - std::sqrt(query.distance_squared);
    const double normal_velocity = particle.vel.dot(query.normal);
    const double force = boundary.material.normal_stiffness * overlap
                       - boundary.material.normal_damping * normal_velocity;

    if (force > 0.0) {
        particle.acc = particle.acc + query.normal * (force / particle.mass);
    }
}
