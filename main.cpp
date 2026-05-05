#include <iostream>
#include "Motion.h"
#include "Output.h"
#include <vector>

// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    // Parameters
    const double dt = 1e-4;
    const double stiffness = 1e7;
    const double dumping = 0.3;
    const Vec3 gravity = {0.0, 0.0, -9.81};
    const double planeAngle = 30.0 * M_PI / 180.0;
    const double planeWidth = 50.0;
    const double planeDepth = 50.0;
    const double planeAffz = 0.5;
    const int numSteps = 50000;
    const int outputInterval = 100;

    // Initialize particles
    std::vector<Particle> particles;
    particles.push_back({{0.0, 0.5, 4}, {0.0, 5.0,1.0}, gravity, 0.4, 50.0});
    particles.push_back({{0.0, 0.6, 1.5}, {0.0, 0.0, 0}, gravity, 0.3, 2.0});
    particles.push_back({{0.0, 0.5, 3}, {0.0, 0.0, 2}, gravity, 0.2, 10.0});

    writePlaneVTK(planeAngle, planeWidth, planeDepth, planeAffz, "plane");

    // Time loop
    for (int step = 0; step < numSteps; step++) {
        velocityVerlet(particles, dt, stiffness, dumping, gravity, planeAngle, planeAffz);

        if (step % outputInterval == 0) {
            writeVTK(particles, step/outputInterval, "particles_");
            std::cout << "Step " << step << std::endl;
        }
    }

    return 0;
}