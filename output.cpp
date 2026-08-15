#include "output.h"
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

void writeBoundaryVTK(const Boundary& boundary, const std::string& prefix) {
    std::filesystem::path dir = create_folder();
    std::ofstream file(dir / (prefix + ".vtk"));
    if (!file) {
        std::cerr << "Failed to open plane VTK file\n";
        return;
    }

    const RectangularPatch& patch = boundary.wall;
    const Vec3 height_axis = patch.heightAxis();
    const Vec3 width_offset = patch.width_axis * (0.5 * patch.width);
    const Vec3 height_offset = height_axis * (0.5 * patch.height);

    const Vec3 p0 = patch.plane.point - width_offset - height_offset;
    const Vec3 p1 = patch.plane.point + width_offset - height_offset;
    const Vec3 p2 = patch.plane.point - width_offset + height_offset;
    const Vec3 p3 = patch.plane.point + width_offset + height_offset;

    file << "# vtk DataFile Version 3.0\n";
    file << "Inclined Plane\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS 4 float\n";
    file << p0.x << " " << p0.y << " " << p0.z << "\n";
    file << p1.x << " " << p1.y << " " << p1.z << "\n";
    file << p2.x << " " << p2.y << " " << p2.z << "\n";
    file << p3.x << " " << p3.y << " " << p3.z << "\n";

    file << "\nPOLYGONS 2 8\n";
    file << "3 0 1 2\n";
    file << "3 1 3 2\n";

    file << "\nPOINT_DATA 4\n";
    file << "NORMALS normals float\n";
    for (int i = 0; i < 4; ++i)
        file << patch.plane.normal.x << " " << patch.plane.normal.y
             << " " << patch.plane.normal.z << "\n";
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
