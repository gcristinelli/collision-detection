//
// Created by Giacomo Cristinelli on 18/12/2025.
//

#include "Output.h"
#include "Motion.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>

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

    // Calculate total points and cells for all spheres
    int resolution = 64;
    int pointsPerSphere = resolution * (resolution - 1) + 2;
    int cellsPerSphere = resolution * (resolution - 1) * 2;
    int totalPoints = particles.size() * pointsPerSphere;
    int totalCells = particles.size() * cellsPerSphere;

    file << "# vtk DataFile Version 3.0\n";
    file << "DEM Particles as Spheres\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << totalPoints << " float\n";

    // Generate sphere vertices for each particle
    for (const auto& p : particles) {
        // Top pole
        file << p.pos.x << " " << (p.pos.y + p.radius) << " " << p.pos.z << "\n";

        // Latitude rings
        for (int i = 1; i < resolution; i++) {
            double theta = M_PI * i / resolution;
            for (int j = 0; j < resolution; j++) {
                double phi = 2.0 * M_PI * j / resolution;
                double x = p.pos.x + p.radius * sin(theta) * cos(phi);
                double y = p.pos.y + p.radius * cos(theta);
                double z = p.pos.z + p.radius * sin(theta) * sin(phi);
                file << x << " " << y << " " << z << "\n";
            }
        }

        // Bottom pole
        file << p.pos.x << " " << (p.pos.y - p.radius) << " " << p.pos.z << "\n";
    }

    // Write triangles
    file << "\nPOLYGONS " << totalCells << " " << (totalCells * 4) << "\n";

    for (size_t p_idx = 0; p_idx < particles.size(); p_idx++) {
        int offset = p_idx * pointsPerSphere;

        // Top cap
        for (int j = 0; j < resolution; j++) {
            int next = (j + 1) % resolution;
            file << "3 " << offset << " " << (offset + 1 + j) << " "
                 << (offset + 1 + next) << "\n";
        }

        // Middle rings
        for (int i = 0; i < resolution - 2; i++) {
            for (int j = 0; j < resolution; j++) {
                int next = (j + 1) % resolution;
                int curr = offset + 1 + i * resolution + j;
                int currNext = offset + 1 + i * resolution + next;
                int below = offset + 1 + (i + 1) * resolution + j;
                int belowNext = offset + 1 + (i + 1) * resolution + next;

                file << "3 " << curr << " " << below << " " << currNext << "\n";
                file << "3 " << currNext << " " << below << " " << belowNext << "\n";
            }
        }

        // Bottom cap
        int bottomPole = offset + pointsPerSphere - 1;
        int lastRingStart = offset + 1 + (resolution - 2) * resolution;
        for (int j = 0; j < resolution; j++) {
            int next = (j + 1) % resolution;
            file << "3 " << bottomPole << " " << (lastRingStart + next) << " "
                 << (lastRingStart + j) << "\n";
        }
    }
}