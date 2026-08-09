#include "motion.h"
#include "contact.h"
#include "boundary.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>


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