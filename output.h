//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include <vector>
#include <string>
#include "Motion.h"

#ifndef DEM_OUTPUT_H
#define DEM_OUTPUT_H

void writeVTK(const std::vector<Particle>& particles, int step,
              const std::string& prefix);

void writePlaneVTK(double angle, double width, double depth,
                   double affz,
                   const std::string& prefix);

#endif //ASSIGNMENT_OUTPUT_H