Simple discrete element method for spheres on an inclined plane, using the Verlet algorithm and sweep-and-prune collision detection.

To build/run:

g++ -std=c++17 main.cpp motion.cpp output.cpp contact.cpp boundary.cpp particles.cpp -o DEM; ./DEM

Particle input is read from `particles.csv` with the required columns:

```text
id,x,y,z,vx,vy,vz,radius,mass
```


To visualise on Paraview:

Select Point Gaussian
Add Filters/Glyph
	Glyph type: Sphere
	Orientation array: radius
	Scale array: radius
	Scale factor: 2
	Glyph Mode: all points
