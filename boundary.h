#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "particles.h"

struct Plane {
	Vec3 point;
	Vec3 normal;

	double sdf(const Vec3& p) const {
		return (p - point).dot(normal);
	}
};

struct RectangularPatch {
    Plane plane;  // plane.point is the patch center
    Vec3 width_axis;
    double width, height;

    Vec3 heightAxis() const { return plane.normal.cross(width_axis); }
};

struct BoundaryMaterial {
    double normal_stiffness;
    double normal_damping;
    double friction_coefficient = 0.0;
};

struct Boundary {
    RectangularPatch wall;
    BoundaryMaterial material;
};

struct SurfaceQuery {
    Vec3 closest_point;
    Vec3 normal;
    double distance_squared;
};

RectangularPatch makeRectangularPatch(const Vec3& center,
                                      const Vec3& normal,
                                      const Vec3& width_direction,
                                      double width, double height);

SurfaceQuery closestPoint(const RectangularPatch&, const Vec3& point);

void applyBoundaryContact(
    Particle& particle,
    const Boundary& boundary,
    const SurfaceQuery& query);

#endif
