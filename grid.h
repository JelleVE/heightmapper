
#ifndef GRID_H
#define GRID_H

// #include <array>

using std::vector;

// struct Grid {
// 	//indices into vectors
// 	std::array<int,3> vertices;

// 	//int material;
// 	Triangle(std::array<int,3> v) : vertices(v) {};
// };


struct Cell {
	//Cell contains a list (vector) of overlapping (2D) triangles
	vector<Triangle> triangles;
	float height=0;
};

#endif