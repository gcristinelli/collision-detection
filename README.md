# 🪨 3D Discrete Element Method (DEM) Particle Simulator

A high-performance, header-light **C++17** implementation of a three-dimensional Discrete Element Method solver for granular particle dynamics. The simulation models rigid spherical particles interacting through penalty-spring contact forces and integrates their equations of motion with a **velocity Verlet** scheme. Output is written in the **VTK legacy format** for direct visualisation in ParaView.

---

## Features

- **Velocity Verlet integration** — second-order accurate, symplectic time-stepping
- **Linear spring–dashpot contact model** — normal stiffness and viscous damping for both particle–particle and particle–wall interactions
- **Sweep-and-prune (SAP) broad-phase contact detection** — O(N log N) average complexity, ready to scale
- **Inclined plane boundary condition** — configurable tilt angle and position
- **VTK output** — high-resolution tessellated sphere geometry + plane mesh, compatible with ParaView and VisIt
- **Automatic timestamped output directories** — results are never overwritten between runs
- **Zero external dependencies** — pure C++ standard library (C++17)

---

## Quick Start

### Prerequisites

| Tool | Minimum version |
|---|---|
| C++ compiler (GCC / Clang / MSVC) | C++17 support |
| CMake *(optional)* | 3.15+ |
| ParaView *(visualisation)* | 5.x+ |

### Build with CMake

```bash
git clone https://github.com/<your-username>/collision-detection.git
cd collision-detection
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Build manually

```bash
g++ -std=c++17 -O2 -o collision-detection main.cpp Motion.cpp Output.cpp
```

### Run

```bash
./collision_detection
```

Output VTK files are written to `results/<YYYY-MM-DD_HH-MM-SS>/`.

---

## Visualisation

1. Open **ParaView**.
2. *File → Open* → select the `particles_*.vtk` file series.
3. Click **Apply**, then optionally add the `plane.vtk` file for the boundary.
4. Use the **Play** button to animate the particle trajectories.

---

## Project Structure

```
.
├── main.cpp        # Entry point — parameters, initialisation, time loop
├── Motion.h        # Data structures: Vec3, Particle, Contact
├── Motion.cpp      # Physics engine: contact detection, forces, Verlet integrator
├── Output.h        # VTK writer declarations
└── Output.cpp      # VTK file generation (spheres + plane)
```

---

## Physics Overview

The solver implements:

**Equations of motion** (Newton's 2nd law per particle):

$$m_i \ddot{\mathbf{x}}_i = \mathbf{F}_i^{\text{grav}} + \mathbf{F}_i^{\text{contact}} + \mathbf{F}_i^{\text{wall}}$$

**Contact force** (linear spring, normal direction):

$$\mathbf{F}_{ij} = k \,\delta_{ij}\,\hat{\mathbf{n}}_{ij}, \qquad \delta_{ij} = r_i + r_j - \|\mathbf{x}_j - \mathbf{x}_i\|$$

**Wall force** (spring–dashpot):

$$F^{\text{wall}} = k\,\delta^{\text{wall}} - \eta\,(\mathbf{v} \cdot \hat{\mathbf{n}}_{\text{wall}})$$

**Velocity Verlet integration**:

$$\mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{v}^n\Delta t + \tfrac{1}{2}\mathbf{a}^n\Delta t^2$$
$$\mathbf{v}^{n+1} = \mathbf{v}^n + \tfrac{1}{2}(\mathbf{a}^n + \mathbf{a}^{n+1})\Delta t$$

---

## Configuration

All physical and numerical parameters are set at the top of `main.cpp`:

| Parameter | Variable | Default |
|---|---|---|
| Time step | `dt` | `1e-4` s |
| Contact stiffness | `stiffness` | `1e7` N/m |
| Wall damping | `dumping` | `0.3` |
| Plane angle | `planeAngle` | `0°` |
| Plane z-intercept | `planeAffz` | `0.5` m |
| Total steps | `numSteps` | `50 000` |
| Output frequency | `outputInterval` | every 100 steps |

---

## Potential Extensions

- [ ] Tangential friction (Coulomb model)
- [ ] Hertzian (non-linear) contact mechanics
- [ ] Uniform grid / octree for O(N) contact detection
- [ ] OpenMP / MPI parallelisation for large-scale simulations
- [ ] Angular momentum and particle rotation

---

## References

- Cundall & Strack, *A discrete numerical model for granular assemblies*, Géotechnique, 1979.
- Verlet, *Computer "Experiments" on Classical Fluids*, Physical Review, 1967.
- Schroeder et al., *The Visualization Toolkit*, Kitware, 2006.

---

## License

This project is released under the [MIT License](LICENSE).
