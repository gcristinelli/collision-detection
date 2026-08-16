#include "contact.h"

#include <algorithm>
#include <numeric>


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
