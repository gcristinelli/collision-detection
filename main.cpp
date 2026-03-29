#include <iostream>
#include <apache2/ap_compat.h>
#include "Motion.h"
#include "Output.h"
#include <array>
#include <vector>

// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    // Parameters
    const double dt = 1e-4;
    const double stiffness = 1e5;
    const double dumping = 0.3;
    const Vec3 gravity = {0.0, 0.0, -9.81};
    const double planeAngle = 15.0 * M_PI / 180.0;
    const double planeWidth = 30.0;
    const double planeDepth = 30.0;
    const int numSteps = 10000;
    const int outputInterval = 30;

    // Initialize particles
    std::vector<Particle> particles;
    particles.push_back({{0.0, 0.5, 0}, {0.0, 1.0, 0.0}, gravity, 0.4, 50.0});
    particles.push_back({{0.0, 0.5, 1.5}, {0.0, 0.0, -3.0}, gravity, 0.3, 2.0});
    particles.push_back({{0.0, 0.5, 3}, {0.0, 0.0, -5.0}, gravity, 0.2, 10.0});

    writePlaneVTK(planeAngle, planeWidth, planeDepth, "plane");

    // Time loop
    for (int step = 0; step < numSteps; step++) {
        velocityVerlet(particles, dt, stiffness, dumping, gravity, planeAngle, 0.0);

        if (step % outputInterval == 0) {
            writeVTK(particles, step/outputInterval, "particles_");
            std::cout << "Step " << step << std::endl;
        }
    }

    return 0;
}