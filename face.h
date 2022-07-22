
#ifndef FACE_H
#define FACE_H

#include <array>

struct Triangle {
	//indices into vectors of coordinates, vertex normals, and uv texture coordinates for each vertex
	std::array<int,3> vertices;

	//int material;
	Triangle(std::array<int,3> v) : vertices(v) {};
};

#endif