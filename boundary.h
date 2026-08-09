#import "motion.h"

#ifndef BOUNDARY_H
#define BOUNDARY_H

void applyInclinedPlane(Particle& p, double angle, double affz,
                        double stiffness, double damping);

#endif