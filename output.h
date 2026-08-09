//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include <vector>
#include <string>
#include "motion.h"

#ifndef OUTPUT_H
#define OUTPUT_H

void writeVTK(const std::vector<Particle>& particles, int step,
              const std::string& prefix);

void writePlaneVTK(double angle, double width, double depth,
                   double affz,
                   const std::string& prefix);

#endif