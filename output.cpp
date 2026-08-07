//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include "output.h"
#include "motion.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::filesystem::path create_folder() {
    namespace fs = std::filesystem;

    // Static: executed only once per program run
    static fs::path outputDir;
    static bool initialized = false;

    if (!initialized) {
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&t);

        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");

        outputDir = fs::path("results") / oss.str();
        fs::create_directories(outputDir);

        initialized = true;
    }

    return outputDir;
}

void writePlaneVTK(double angle, double width, double depth,
                   double affz,
                   const std::string& prefix) {
    std::filesystem::path dir = create_folder();
    std::ofstream file(dir / (prefix + ".vtk"));
    if (!file) {
        std::cerr << "Failed to open plane VTK file\n";
        return;
    }

    double tanA  = std::tan(angle);
    double halfW = width / 2.0;
    double halfD = depth / 2.0;

    // Plane: z = tanA*x + affz
    double x0 = -halfW;
    double x1 =  halfW;
    double y0 = -halfD;
    double y1 =  halfD;

    double z_x0 = tanA * x0 + affz;
    double z_x1 = tanA * x1 + affz;

    Vec3 normal = {-tanA, 0.0, 1.0};
    normal = normal.normalized();

    file << "# vtk DataFile Version 3.0\n";
    file << "Inclined Plane\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS 4 float\n";
    file << x0 << " " << y0 << " " << z_x0 << "\n";
    file << x1 << " " << y0 << " " << z_x1 << "\n";
    file << x0 << " " << y1 << " " << z_x0 << "\n";
    file << x1 << " " << y1 << " " << z_x1 << "\n";

    file << "\nPOLYGONS 2 8\n";
    file << "3 0 1 2\n";
    file << "3 1 3 2\n";

    file << "\nPOINT_DATA 4\n";
    file << "NORMALS normals float\n";
    for (int i = 0; i < 4; ++i)
        file << normal.x << " " << normal.y << " " << normal.z << "\n";
}


void writeVTK(const std::vector<Particle>& particles, int step,
              const std::string& prefix) {
    std::filesystem::path dir = create_folder();
    std::ofstream file(dir / (prefix + std::to_string(step) + ".vtk"));
    if (!file) {
        std::cerr << "Failed to open VTK file\n";
        return;
    }

    const std::size_t n = particles.size();
    file << "# vtk DataFile Version 3.0\n";
    file << "DEM particle centres\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS " << n << " float\n";
    for (const auto& p : particles)
        file << p.pos.x << " " << p.pos.y << " " << p.pos.z << "\n";

    // Vertex cells make every particle centre an explicit renderable point.
    file << "\nVERTICES " << n << " " << 2 * n << "\n";
    for (std::size_t i = 0; i < n; ++i)
        file << "1 " << i << "\n";

    file << "\nPOINT_DATA " << n << "\n";
    file << "SCALARS radius float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles)
        file << p.radius << "\n";

    file << "VECTORS velocity float\n";
    for (const auto& p : particles)
        file << p.vel.x << " " << p.vel.y << " " << p.vel.z << "\n";

    file << "SCALARS mass float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles)
        file << p.mass << "\n";
}