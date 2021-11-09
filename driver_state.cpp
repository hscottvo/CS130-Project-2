#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    std::cout<<"TODO: allocate and initialize state.image_depth."<<std::endl;
    state.image_color = new pixel[width*height];
    for (unsigned int i = 0; i < width*height; ++i){
        state.image_color[i] = make_pixel(0, 0, 0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry* geometry_arr = new data_geometry[state.num_vertices];
    for (unsigned int i = 0; i < state.num_vertices; ++i) {
        geometry_arr[i].data = new float[MAX_FLOATS_PER_VERTEX];
        data_vertex vertex;
        // std::cout << "TODO: Check if vertex.data is correct" << std::endl;
        vertex.data = &(state.vertex_data[i * state.floats_per_vertex]);
        state.vertex_shader(vertex, geometry_arr[i], state.uniform_data);
        // geometry_arr[i].gl_Position = 
    }

    switch(type){ 
        case render_type::triangle:

            for(unsigned int i = 0; i < state.num_vertices-2; i += 3) {
            // for(unsigned int i = 0; i < state.num_vertices*state.floats_per_vertex; i += state.floats_per_vertex) {
                // rasterize_triangle(state, 
                //                    geometry_arr[i],
                //                    geometry_arr[i+1],
                //                    geometry_arr[i+2]);
                clip_triangle(state, 
                              geometry_arr[i],
                              geometry_arr[i+1],
                              geometry_arr[i+2],
                              6);

            }
        break;

        case render_type::indexed:

        break;

        case render_type::fan:

        break;

        case render_type::strip:

        break;
    }

    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    vec3 point_a = vec3(v0.gl_Position[0]/v0.gl_Position[3], 
                        v0.gl_Position[1]/v0.gl_Position[3], 
                        v0.gl_Position[2]/v0.gl_Position[3]);
    point_a[0] = ((state.image_width-1)*((point_a[0]+1)/2)); 
    point_a[1] = ((state.image_height-1)*((point_a[1]+1)/2));

    vec3 point_b = vec3(v1.gl_Position[0]/v1.gl_Position[3], 
                        v1.gl_Position[1]/v1.gl_Position[3], 
                        v1.gl_Position[2]/v1.gl_Position[3]);
    point_b[0] = ((state.image_width-1)*((point_b[0]+1)/2)); 
    point_b[1] = ((state.image_height-1)*((point_b[1]+1)/2));

    vec3 point_c = vec3(v2.gl_Position[0]/v2.gl_Position[3], 
                        v2.gl_Position[1]/v2.gl_Position[3], 
                        v2.gl_Position[2]/v2.gl_Position[3]);
    point_c[0] = ((state.image_width-1)*((point_c[0]+1)/2));
    point_c[1] = ((state.image_height-1)*((point_c[1]+1)/2));

    // can get bounding box of triangle and iterate over that instead of the whole image
    float total_area = 0.5 * ((point_b[0]*point_c[1] - point_c[0]*point_b[1]) + 
                              (point_c[0]*point_a[1]-point_a[0]*point_c[1]) + 
                              (point_a[0]*point_b[1]-point_b[0]*point_a[1]));
    float alpha;
    float beta;
    float gamma;
    int x_min = std::min(std::min(point_a[0], point_b[0]), point_c[0]);
    int y_min = std::min(std::min(point_a[1], point_b[1]), point_c[1]);

    int x_max = std::max(std::max(point_a[0], point_b[0]), point_c[0]);
    int y_max = std::max(std::max(point_a[1], point_b[1]), point_c[1]);

    // for (unsigned int i = 0; i < state.image_width; ++i) {
    //     for (unsigned int j = 0; j < state.image_height; ++j) {
    for (unsigned int i = x_min; i < x_max; ++i) {
        for (unsigned int j = y_min; j < y_max; ++j) {
            vec3 point_p = vec3(i, j, 0);
            alpha = (0.5 * ((point_b[0]*point_c[1] - point_c[0]*point_b[1]) + 
                         (point_c[0]*point_p[1]-point_p[0]*point_c[1]) + 
                         (point_p[0]*point_b[1]-point_b[0]*point_p[1])) / total_area);
            beta = (0.5 * ((point_p[0]*point_c[1] - point_c[0]*point_p[1]) + 
                        (point_c[0]*point_a[1]-point_a[0]*point_c[1]) +
                        (point_a[0]*point_p[1]-point_p[0]*point_a[1])) / total_area);
            gamma = (0.5 * ((point_b[0]*point_p[1]-point_p[0]*point_b[1]) + 
                        (point_p[0]*point_a[1]-point_a[0]*point_p[1]) +
                        (point_a[0]*point_b[1]-point_b[0]*point_a[1])) / total_area);

            

            if (alpha > 0 && beta > 0 && gamma > 0) {
                state.image_color[j*state.image_width+i] = make_pixel(255, 255, 255);
            }
        }
    }
    // std::cout<<"TODO: implement rasterization"<<std::endl;

}

