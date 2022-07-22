#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc

#include "tinyobj_loader.h"
#include "face.h"
#include "grid.h"

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>

using glm::vec2;
using glm::vec3;
using glm::uvec3;
using std::vector;

typedef vector< vector<Cell> > Grid;

const int MAX_H = 1000;
const int MAX_W = 1000;
const float RESOLUTION = 0.05;

////////////////////////////////////////
// https://github.com/canmom/rasteriser/blob/master/fileloader.cpp
////////////////////////////////////////


void components_to_vec2s(const vector<float> components, vector<vec2>& vecs) {
    //convert a vector of back-to-back vertex components to a vector of vec2 objects
    for(size_t vec_start = 0; vec_start < components.size(); vec_start+=2) {
        vecs.push_back(
            vec2(components[vec_start],
                components[vec_start+1]
            ));
    }
}


void components_to_vec3s(const vector<float> components, vector<vec3>& vecs) {
    //convert a vector of back-to-back vertex components to a vector of vec3 objects
    for(size_t vec_start = 0; vec_start < components.size(); vec_start+=3) {
        vecs.push_back(
            vec3(components[vec_start],
                components[vec_start+1],
                components[vec_start+2]
            ));
    }
}


void load_triangles(const tinyobj::shape_t & shape, vector<Triangle> & triangles) {
    //convert a tinyobjloader shape_t object containing indices into vertex properties and textures
    //into a vector of Triangle objects grouping these indices
    const vector<tinyobj::index_t> & indices = shape.mesh.indices;
    const vector<int> & mat_ids = shape.mesh.material_ids;

    for(size_t face_ind = 0; face_ind < mat_ids.size(); face_ind++) {
        triangles.push_back(
            Triangle(
                {indices[3*face_ind].vertex_index, indices[3*face_ind+1].vertex_index, indices[3*face_ind+2].vertex_index}
                ));
    }
}


void get_bounding_box(const vector<vec3>& vertices, float& min_x, float& max_x, float& min_y, float& max_y, const float& resolution) {
    min_x = 1e10, min_y = 1e10;
    max_x = -1e10, max_y = -1e10;

    for (size_t v_ind=0; v_ind < vertices.size(); v_ind++) {
        if (vertices[v_ind][0] < min_x) {
            min_x = vertices[v_ind][0];
        }
        if (vertices[v_ind][0] > max_x) {
            max_x = vertices[v_ind][0];
        }
        if (vertices[v_ind][1] < min_y) {
            min_y = vertices[v_ind][1];
        }
        if (vertices[v_ind][1] > max_y) {
            max_y = vertices[v_ind][1];
        }
    }
    min_x = std::floor(min_x*10)/10;
    min_y = std::floor(min_y*10)/10;

    max_x = std::ceil(max_x*10)/10;
    max_y = std::ceil(max_y*10)/10;
}


void initialize_grid(Grid & grid, vector<Triangle> & triangles, vector<vec3> & vertices, float grid_min_x, float grid_max_x, float grid_min_y, float grid_max_y, float resolution) {
    int nb_cols = grid.size();
    int nb_rows = grid[0].size();

    for (int i=0; i < triangles.size(); i++) {
        // Get triangle vertices
        vec3 v1 = vertices[triangles[i].vertices[0]];
        vec3 v2 = vertices[triangles[i].vertices[1]];
        vec3 v3 = vertices[triangles[i].vertices[2]];

        // Get bounding box triangle
        float min_x = std::min({v1[0], v2[0], v3[0]});
        float max_x = std::max({v1[0], v2[0], v3[0]});
        float min_y = std::min({v1[1], v2[1], v3[1]});
        float max_y = std::max({v1[1], v2[1], v3[1]});
        
        // Get overlapping cell indices
        int min_col_ind = std::floor((min_x - grid_min_x)/resolution);
        int max_col_ind = std::floor((max_x - grid_min_x)/resolution);
        int min_row_ind = std::floor((min_y - grid_min_y)/resolution);
        int max_row_ind = std::floor((max_y - grid_min_y)/resolution);

        for (int ind_x=min_col_ind; ind_x <= max_col_ind; ind_x++) {
            for (int ind_y=min_row_ind; ind_y <= max_row_ind; ind_y++) {
                grid[ind_x][ind_y].triangles.push_back(triangles[i]);
            }
        }
    }
}


float sign (float p1x, float p1y, float p2x, float p2y, float p3x, float p3y) {
    return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
}


bool point_in_triangle (float px, float py, float v1x, float v1y, float v2x, float v2y, float v3x, float v3y) {
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(px, py, v1x, v1y, v2x, v2y);
    d2 = sign(px, py, v2x, v2y, v3x, v3y);
    d3 = sign(px, py, v3x, v3y, v1x, v1y);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}


float calc_z(float x, float y, float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
    float z = (z3*(x-x1)*(y-y2) + z1*(x-x2)*(y-y3) + z2*(x-x3)*(y-y1) - z2*(x-x1)*(y-y3) - z3*(x-x2)*(y-y1) - z1*(x-x3)*(y-y2)) / ((x-x1)*(y-y2) + (x-x2)*(y-y3) + (x-x3)*(y-y1) - (x-x1)*(y-y3) - (x-x2)*(y-y1) - (x-x3)*(y-y2));
    return z;
}


void writeImage(int image[MAX_H][MAX_W], int height, int width) {
    std::ofstream ostr;
    ostr.open("outImage.pgm");
    if (ostr.fail()) {
        std::cout << "Unable to write file\n";
        exit(1);
    };

    // print the header
    ostr << "P2" << std::endl;
    // width, height
    ostr << width << ' ';
    ostr << height << std::endl;
    ostr << 255 << std::endl;

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            assert(image[row][col] < 256);
            assert(image[row][col] >= 0);
            ostr << image[row][col] << ' ';
            ostr << std::endl;
        }
    }
    ostr.close();
    return;
}


int main()
{
    const std::string MODEL_PATH = "finc.obj";
    vector<vec3> vertices;
    vector<Triangle> triangles;

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    bool success = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, MODEL_PATH.c_str());

    if (!err.empty()) {
        std::cerr << err << std::endl;
    }
    if (!success) {
        exit(1);
    }

    components_to_vec3s(attrib.vertices, vertices);

    for(auto shape = shapes.begin(); shape < shapes.end(); ++shape) {
        load_triangles(*shape, triangles);
    }

    //Calculate bounding box
    float min_x, max_x, min_y, max_y;
    get_bounding_box(vertices, min_x, max_x, min_y, max_y, RESOLUTION);

    //Create grid with reference to potential triangles in each cell
    //Uses bounding box of each 2D triangle
    int nb_rows = int((max_y-min_y)*(1/RESOLUTION))+1;
    int nb_cols = int((max_x-min_x)*(1/RESOLUTION))+1;

    Grid grid(nb_cols,vector<Cell>(nb_rows)); // access with x,y coordinates
    initialize_grid(grid, triangles, vertices, min_x, max_x, min_y, max_y, RESOLUTION);

    //Perform intersection between each center cell point and relevant triangles
    for (int col_ind=0; col_ind < nb_cols; col_ind++) {
        for (int row_ind=0; row_ind < nb_rows; row_ind++) {
            float x = min_x + col_ind*RESOLUTION + 0.5*RESOLUTION; //middle x of cell
            float y = min_y + row_ind*RESOLUTION + 0.5*RESOLUTION; //midle y of cell
            vector<Triangle> current_triangles = grid[col_ind][row_ind].triangles;

            // For each triangle, check if it indeed contains the cell's centerpoint
            // If so, calculate z-value
            for (auto triangle=current_triangles.begin(); triangle < current_triangles.end(); ++triangle) {
                float v1x = vertices[(*triangle).vertices[0]][0];
                float v1y = vertices[(*triangle).vertices[0]][1];
                float v1z = vertices[(*triangle).vertices[0]][2];

                float v2x = vertices[(*triangle).vertices[1]][0];
                float v2y = vertices[(*triangle).vertices[1]][1];
                float v2z = vertices[(*triangle).vertices[1]][2];

                float v3x = vertices[(*triangle).vertices[2]][0];
                float v3y = vertices[(*triangle).vertices[2]][1];
                float v3z = vertices[(*triangle).vertices[2]][2];

                bool in_triangle = point_in_triangle(x, y, v1x, v1y, v2x, v2y, v3x, v3y);
                if (in_triangle) {
                    float z = calc_z(x, y, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);
                    if (z > grid[col_ind][row_ind].height) {
                        grid[col_ind][row_ind].height = z;
                    }
                }
            }
        }
    }

    //Grid heights set at this point
    //Write image
    float max_z = 0;
    for (int col_ind=0; col_ind < nb_cols; col_ind++) {
        for (int row_ind=0; row_ind < nb_rows; row_ind++) {
            if (grid[col_ind][row_ind].height > max_z) {
                max_z = grid[col_ind][row_ind].height;
            }
        }
    }

    //Create image array
    int image[MAX_H][MAX_W];
    for (int i=- 0; i < MAX_H; i++) {
        for (int j=0; j < MAX_W; j++) {
            image[i][j] = 0;
        }
    }

    for (int col_ind=0; col_ind < nb_cols; col_ind++) {
        for (int row_ind=0; row_ind < nb_rows; row_ind++) {
                image[col_ind][row_ind] = int(255*grid[col_ind][row_ind].height/max_z);
        }
    }

    writeImage(image, nb_cols, nb_rows);

    return 0;
}