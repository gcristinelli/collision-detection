#include <iostream>
#include <cmath>
#include "motion.h"
#include "boundary.h"
#include "output.h"
#include "particles.h"
#include <vector>

int main() {
    // Parameters
    const double dt = 1e-4;
    const double stiffness = 1e7;
    const double damping = 0.6;
    const Vec3 gravity = {0.0, 0.0, -9.81};
    const double planeAngle = 5.0 * std::acos(-1.0) / 180.0;
    const double planeWidth = 40.0;
    const double planeDepth = 40.0;
    const double planeAffz = 0.5;
    const double wallHeight = 10.0;
    const int numSteps = 50000;
    const int outputInterval = 200;

    const BoundaryMaterial wallMaterial = {stiffness, damping};
    const std::vector<Boundary> boundaries = {
        // Inclined floor: z = tan(planeAngle) * x + planeAffz.
        {makeRectangularPatch(
             {0.0, 0.0, planeAffz},
             {-std::sin(planeAngle), 0.0, std::cos(planeAngle)},
             {std::cos(planeAngle), 0.0, std::sin(planeAngle)},
             planeWidth / std::cos(planeAngle), planeDepth),
         wallMaterial},

        // Four vertical walls enclosing the floor footprint.
        {makeRectangularPatch(
             {-0.5 * planeWidth, 0.0, planeAffz},
             {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, planeDepth, wallHeight),
         wallMaterial},
        {makeRectangularPatch(
             {0.5 * planeWidth, 0.0, planeAffz},
             {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, planeDepth, wallHeight),
         wallMaterial},
        {makeRectangularPatch(
             {0.0, -0.5 * planeDepth, planeAffz},
             {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, planeWidth, wallHeight),
         wallMaterial},
        {makeRectangularPatch(
             {0.0, 0.5 * planeDepth, planeAffz},
             {0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, planeWidth, wallHeight),
         wallMaterial}
    };

    std::vector<Particle> particles = loadParticlesCsv("particles.csv", gravity);

    writeBoundariesVTK(boundaries, "boundaries");

    // Time loop
    for (int step = 0; step < numSteps; step++) {
        velocityVerlet(particles, dt, stiffness, gravity, boundaries);

        if (step % outputInterval == 0) {
            writeVTK(particles, step/outputInterval, "particles_");
            std::cout << "Step " << step << std::endl;
        }
    }

    return 0;
}
