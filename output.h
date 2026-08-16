#include <vector>
#include <string>
#include "motion.h"
#include "boundary.h"

#ifndef OUTPUT_H
#define OUTPUT_H

void writeVTK(const std::vector<Particle>& particles, int step,
              const std::string& prefix);

void writeBoundariesVTK(const std::vector<Boundary>& boundaries,
                        const std::string& prefix);

#endif
