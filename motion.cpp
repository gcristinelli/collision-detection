//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include "motion.h"
#include "contact.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

void computeForces(std::vector<Particle>& particles,
                   const std::vector<Contact>& contacts,
                   double stiffness, const Vec3& gravity) {
    // Reset accelerations
    for (auto& p : particles) {
        p.acc = gravity;
    }

    // Contact forces
    for (const auto& c : contacts) {
        Vec3 n = (particles[c.j].pos - particles[c.i].pos).normalized();
        double force = stiffness * c.overlap;

        particles[c.i].acc = particles[c.i].acc - n * (force / particles[c.i].mass);
        particles[c.j].acc = particles[c.j].acc + n * (force / particles[c.j].mass);
    }
}

void applyInclinedPlane(Particle& p, double angle, double affz,
                        double stiffness, double damping) {
    double tanA = std::tan(angle);

    Vec3 normal = {-tanA, 0.0, 1.0};
    normal = normal.normalized();

    double dist = (p.pos.z - tanA * p.pos.x - affz)
                  / std::sqrt(1.0 + tanA * tanA);

    if (dist < p.radius) {
        double overlap = p.radius - dist;
        double vn = p.vel.dot(normal);

        double forceMag = stiffness * overlap - damping * vn;
        if (forceMag < 0.0) forceMag = 0.0;

        p.acc = p.acc + normal * (forceMag / p.mass);
    }
}

void velocityVerlet(std::vector<Particle>& particles, double dt,
                    double stiffness, double dumping, const Vec3& gravity,
                    double planeAngle, double planeZ) {

    // Step 1: Update positions and save accelerations
    std::vector<Vec3> accOld(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        accOld[i]          = particles[i].acc;
        particles[i].pos   = particles[i].pos + particles[i].vel * dt + particles[i].acc * (0.5 * dt * dt);
    }

    // Step 2: Detect contacts and update accelerations
    auto contacts = detectContact_SP(particles);
    computeForces(particles, contacts, stiffness, gravity);

    // Apply boundary conditions
    for (auto& p : particles) {
        applyInclinedPlane(p, planeAngle, planeZ, stiffness, dumping);
    }

    // Step 3: Complete velocity Verlet algorithm by averaging new and old acceleration for the velocity
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].vel = particles[i].vel + (accOld[i] + particles[i].acc) * (0.5 * dt);
    }
}