#include <vector>
#include <cmath>

#ifndef MOTION_H
#define MOTION_H

struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& v) const { return {x+v.x, y+v.y, z+v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x-v.x, y-v.y, z-v.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    double length() const { return sqrt(x*x + y*y + z*z); }
    double dot(const Vec3& v)  const { return x*v.x + y*v.y + z*v.z; }
    Vec3 cross(const Vec3& v) const {
        return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }
    Vec3 normalized() const { double l = length(); return {x/l, y/l, z/l}; }
};

struct Particle {
    Vec3 pos, vel, acc;
    double radius, mass;
};

struct Boundary;

void velocityVerlet(std::vector<Particle>& particles, double dt,
                    double stiffness, const Vec3& gravity,
                    const Boundary& boundary);

#endif
