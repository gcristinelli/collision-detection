Simple discrete element method for spheres on an inclined plane, using the Verlet algorithm and sweep-and-prune collision detection.

To build/run:

g++ -std=c++17 main.cpp motion.cpp output.cpp contact.cpp boundary.cpp -o DEM; ./DEM


To visualise on Paraview:

Select Point Gaussian
Add Filters/Glyph
	Glyph type: Sphere
	Orientation array: radius
	Scale array: radius
	Scale factor: 2
	Glyph Mode: all points