#include <iostream>
#include <cmath>
#include "motion.h"
#include "boundary.h"
#include "output.h"
#include <vector>

int main() {
    // Parameters
    const double dt = 1e-4;
    const double stiffness = 1e7;
    const double damping = 0.3;
    const Vec3 gravity = {0.0, 0.0, -9.81};
    const double planeAngle = 15.0 * std::acos(-1.0) / 180.0;
    const double planeWidth = 50.0;
    const double planeDepth = 50.0;
    const double planeAffz = 0.5;
    const int numSteps = 50000;
    const int outputInterval = 100;

    const Boundary plane = {
        makeRectangularPatch(
            {0.0, 0.0, planeAffz},
            {-std::sin(planeAngle), 0.0, std::cos(planeAngle)},
            {std::cos(planeAngle), 0.0, std::sin(planeAngle)},
            planeWidth / std::cos(planeAngle), planeDepth),
        {stiffness, damping}
    };

    // Initialize particles
    std::vector<Particle> particles;
    particles.push_back({{-0.9, 0.5, 3.5}, {0.0, 0.0,0.0}, gravity, 0.4, 150.0});
    particles.push_back({{0.0, 0.5, 6}, {0.0, 0.0,-3.0}, gravity, 0.3, 50.0});
    particles.push_back({{0.0, 0.5, 4}, {0.0, 0.0,0.0}, gravity, 0.4, 50.0});
    particles.push_back({{0.0, 0.6, 1.5}, {0.0, 0.0, 0}, gravity, 0.3, 2.0});
    particles.push_back({{0.0, 0.5, 3}, {0.0, 0.0, 2}, gravity, 0.2, 10.0});

    writeBoundaryVTK(plane, "plane");

    // Time loop
    for (int step = 0; step < numSteps; step++) {
        velocityVerlet(particles, dt, stiffness, gravity, plane);

        if (step % outputInterval == 0) {
            writeVTK(particles, step/outputInterval, "particles_");
            std::cout << "Step " << step << std::endl;
        }
    }

    return 0;
}
