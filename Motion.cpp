//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include "Motion.h"
#include <cmath>
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

void applyInclinedPlane(Particle& p, double angle, double planeY, double stiffness) {
    // Plane equation: y = tan(angle) * x + planeY
    double normalY = cos(angle);
    double normalX = -sin(angle);
    double dist = (p.pos.y - planeY - tan(angle) * p.pos.x) * normalY;

    if (dist < p.radius) {
        double overlap = p.radius - dist;
        Vec3 normal = {normalX, normalY, 0};
        double force = stiffness * overlap / p.mass;
        p.acc = p.acc + normal * force;
    }
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
                    double stiffness, const Vec3& gravity,
                    double planeAngle, double planeY) {

    // Step 1: Update positions and half-step velocities
    for (auto& p : particles) {
        p.pos = p.pos + p.vel * dt + p.acc * (0.5 * dt * dt);
        //p.vel = p.vel + p.acc * (0.5 * dt);
    }

    // Step 2: Detect contacts and compute new forces
    auto contacts = detectContact(particles);
    computeForces(particles, contacts, stiffness, gravity);

    // Apply boundary conditions
    for (auto& p : particles) {
        applyInclinedPlane(p, planeAngle, planeY, stiffness);
    }

    // Step 3: Complete velocity update
    for (auto& p : particles) {
        p.vel = p.vel + p.acc * (0.5 * dt);
    }
}