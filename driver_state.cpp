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
    for (int i = 0; i < width*height; ++i){
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
    for (int i = 0; i < state.num_vertices; ++i) {
        geometry_arr[i].data = new float[MAX_FLOATS_PER_VERTEX];
        data_vertex vertex;
        vertex.data = &(state.vertex_data[i * state.floats_per_vertex]);
        state.vertex_shader(vertex, geometry_arr[i], state.uniform_data);
    }

    switch(type){ 
        case render_type::triangle:

            for(int i = 0; i < state.num_vertices-2; i += 3) {
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

        default:
            exit(0);
    }
    for(int i = 0; i < state.num_vertices; ++i) {
        delete[] geometry_arr[i].data;
    }
    delete[] geometry_arr;
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
    point_a[0] = ((state.image_width-1)*((point_a[0]+1)/2.0)); 
    point_a[1] = ((state.image_height-1)*((point_a[1]+1)/2.0));

    vec3 point_b = vec3(v1.gl_Position[0]/v1.gl_Position[3], 
                        v1.gl_Position[1]/v1.gl_Position[3], 
                        v1.gl_Position[2]/v1.gl_Position[3]);
    point_b[0] = ((state.image_width-1)*((point_b[0]+1)/2.0)); 
    point_b[1] = ((state.image_height-1)*((point_b[1]+1)/2.0));

    vec3 point_c = vec3(v2.gl_Position[0]/v2.gl_Position[3], 
                        v2.gl_Position[1]/v2.gl_Position[3], 
                        v2.gl_Position[2]/v2.gl_Position[3]);
    point_c[0] = ((state.image_width-1)*((point_c[0]+1)/2.0));
    point_c[1] = ((state.image_height-1)*((point_c[1]+1)/2.0));

    // can get bounding box of triangle and iterate over that instead of the whole image
    float total_area = 0.5 * ((point_b[0]*point_c[1] - point_c[0]*point_b[1]) + 
                              (point_c[0]*point_a[1]-point_a[0]*point_c[1]) + 
                              (point_a[0]*point_b[1]-point_b[0]*point_a[1]));
    float alpha;
    float beta;
    float gamma;

    int x_min = std::max(std::min(std::min(point_a[0], point_b[0]), point_c[0]) - 1, float(0));
    int y_min = std::max(std::min(std::min(point_a[1], point_b[1]), point_c[1]) - 1, float(0));

    int x_max = std::min(std::max(std::max(point_a[0], point_b[0]), point_c[0]) + 1, float(state.image_width));
    int y_max = std::min(std::max(std::max(point_a[1], point_b[1]), point_c[1]) + 1, float(state.image_height));



    // for (int i = 0; i < state.image_width; ++i) {
    //     for (int j = 0; j < state.image_height; ++j) {
    for (int i = x_min; i < x_max; ++i) {
        for (int j = y_min; j < y_max; ++j) {
            // figure out what to put for z 
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

            

            if (alpha > -.001 && beta > -0.001 && gamma > -0.001) {
            
                data_fragment fragment_color;
                fragment_color.data = new float[MAX_FLOATS_PER_VERTEX];
                data_output pixel_color;
                state.fragment_shader(fragment_color, pixel_color, state.uniform_data);

                int pixel_r = pixel_color.output_color[0]*255;
                int pixel_g = pixel_color.output_color[1]*255;
                int pixel_b = pixel_color.output_color[2]*255;


                if (state.floats_per_vertex > 3) {
                    switch(state.interp_rules[3]){
                        case interp_type::flat:
                            pixel_r = v0.data[3]*255;
                            pixel_g = v0.data[4]*255;
                            pixel_b = v0.data[5]*255;
                        break;

                        case interp_type::smooth:

                        break;

                        case interp_type::noperspective:
                            pixel_r = v0.data[3]*255*alpha + v1.data[3]*255*beta + v2.data[3]*255*gamma;
                            pixel_g = v0.data[4]*255*alpha + v1.data[4]*255*beta + v2.data[4]*255*gamma;
                            pixel_b = v0.data[5]*255*alpha + v1.data[5]*255*beta + v2.data[5]*255*gamma;
                        break;

                        default:
                            std::cout << "Incorrect interp_rules" << std::endl;
                            exit(0);
                        break;
                    }
                } 

                state.image_color[j*state.image_width+i] = make_pixel(pixel_r, 
                                                                      pixel_g, 
                                                                      pixel_b);
                for (int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                    delete fragment_color.data[i];
                }
                delete[] fragment_color.data;
                // std::cout << "r: " << pixel_color.output_color[0] << " g: " << pixel_color.output_color[1]
                //         << " b: " << pixel_color.output_color[2] << std::endl;

                
                // state.image_color[j*state.image_width+i] = make_pixel(255, 255, 255);
            }
        }
    }
    // std::cout<<"TODO: implement rasterization"<<std::endl;

}

