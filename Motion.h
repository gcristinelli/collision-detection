//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include <vector>

#ifndef ASSIGNMENT_MOTION_H
#define ASSIGNMENT_MOTION_H

struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& v) const { return {x+v.x, y+v.y, z+v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x-v.x, y-v.y, z-v.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    double length() const { return sqrt(x*x + y*y + z*z); }
    Vec3 normalized() const { double l = length(); return {x/l, y/l, z/l}; }
};

struct Particle {
    Vec3 pos, vel, acc;
    double radius, mass;
};

struct Contact {
    int i, j;
    double overlap;
};

void velocityVerlet(std::vector<Particle>& particles, double dt,
                    double stiffness, const Vec3& gravity,
                    double planeAngle, double planeY);

#endif //ASSIGNMENT_MOTION_H