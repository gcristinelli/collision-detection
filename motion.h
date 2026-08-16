#ifndef MOTION_H
#define MOTION_H

#include <vector>

#include "boundary.h"
#include "particles.h"
#include "utils.h"

void velocityVerlet(std::vector<Particle>& particles, double dt,
                    double stiffness, const Vec3& gravity,
                    const std::vector<Boundary>& boundaries);

#endif
