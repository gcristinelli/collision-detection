#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>

#include "boundary.h"
#include "particles.h"

void writeVTK(const std::vector<Particle>& particles, int step,
              const std::string& prefix);

void writeBoundariesVTK(const std::vector<Boundary>& boundaries,
                        const std::string& prefix);

#endif
