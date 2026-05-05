//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include "Motion.h"
#include <cmath>
#include <vector>
#include <numeric>

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

std::vector<Contact> detectContact_SP(std::vector<Particle>& particles) {
    std::vector<Contact> contacts;

    std::vector<size_t> ids(particles.size());
    std::iota(ids.begin(), ids.end(), 0);

    std::sort(ids.begin(), ids.end(), [&](size_t a, size_t b) {
        return particles[a].pos.x - particles[a].radius <
               particles[b].pos.x - particles[b].radius;
    });

    for (size_t a = 0; a < ids.size(); ++a) {
        size_t i = ids[a];
        double iMaxX = particles[i].pos.x + particles[i].radius;

        for (size_t b = a + 1; b < ids.size(); ++b) {
            size_t j = ids[b];
            double jMinX = particles[j].pos.x - particles[j].radius;

            if (jMinX > iMaxX) {
                break;
            }

            Vec3 delta = particles[j].pos - particles[i].pos;
            double dist = delta.length();
            double sumRadii = particles[i].radius + particles[j].radius;

            if (dist < sumRadii) {
                Contact c;
                c.i = i;
                c.j = j;
                c.overlap = sumRadii - dist;
                contacts.push_back(c);
            }
        }
    }
    return contacts;
}


std::vector<Contact> detectContact(std::vector<Particle>& particles) {
    std::vector<Contact> contacts;
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Vec3 delta = particles[j].pos - particles[i].pos;
            double dist = delta.length();
            double sumRadii = particles[i].radius + particles[j].radius;
            if (dist < sumRadii) {
                Contact c;
                c.i = i;
                c.j = j;
                c.overlap = sumRadii - dist;
                contacts.push_back(c);
            }
        }

    }
    return contacts;
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