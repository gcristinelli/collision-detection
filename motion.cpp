#include "motion.h"
#include "contact.h"
#include "boundary.h"
#include <vector>

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

void velocityVerlet(std::vector<Particle>& particles, double dt,
                    double stiffness, const Vec3& gravity,
                    const std::vector<Boundary>& boundaries) {

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
        for (const Boundary& boundary : boundaries) {
            const SurfaceQuery query = closestPoint(boundary.wall, p.pos);
            if (query.distance_squared < p.radius * p.radius) {
                applyBoundaryContact(p, boundary, query);
            }
        }
    }

    // Step 3: Complete velocity Verlet algorithm by averaging new and old acceleration for the velocity
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].vel = particles[i].vel + (accOld[i] + particles[i].acc) * (0.5 * dt);
    }
}
