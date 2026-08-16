#include "particles.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<double> parseRow(const std::string& line, std::size_t line_number) {
    std::stringstream stream(line);
    std::vector<double> values;
    std::string field;

    while (std::getline(stream, field, ',')) {
        values.push_back(std::stod(field));
    }

    if (values.size() != 9) {
        throw std::runtime_error("Expected 9 columns in particles CSV at line "
                                 + std::to_string(line_number));
    }

    for (double value : values) {
        if (!std::isfinite(value)) {
            throw std::runtime_error("Non-finite value in particles CSV at line "
                                     + std::to_string(line_number));
        }
    }

    return values;
}

} // namespace

std::vector<Particle> loadParticlesCsv(const std::string& filename,
                                       const Vec3& initial_acceleration) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open particles CSV: " + filename);
    }

    std::vector<Particle> particles;
    std::string line;
    std::size_t line_number = 0;
    bool header_read = false;

    while (std::getline(file, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') {
            continue;
        }
        if (!header_read) {
            if (line != "id,x,y,z,vx,vy,vz,radius,mass") {
                throw std::runtime_error("Invalid particles CSV header");
            }
            header_read = true;
            continue;
        }

        const std::vector<double> values = parseRow(line, line_number);
        const double radius = values[7];
        const double mass = values[8];
        if (radius <= 0.0 || mass <= 0.0) {
            throw std::runtime_error("Radius and mass must be positive at line "
                                     + std::to_string(line_number));
        }

        particles.push_back({
            {values[1], values[2], values[3]},
            {values[4], values[5], values[6]},
            initial_acceleration,
            radius,
            mass
        });
    }

    if (!header_read) {
        throw std::runtime_error("Particles CSV is missing its header");
    }
    if (particles.empty()) {
        throw std::runtime_error("Particles CSV contains no particles");
    }

    return particles;
}
