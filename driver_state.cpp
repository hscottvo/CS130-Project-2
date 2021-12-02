#include "driver_state.h"
#include <cstring>
#include <limits>

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
    state.image_color = new pixel[width*height];
    state.image_depth = new float[width*height];
    for (int i = 0; i < width*height; ++i){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = std::numeric_limits<float>::max();
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
                // std::cout << "v0: " << geometry_arr[i].gl_Position << std::endl;  
                // std::cout << "v1: " << geometry_arr[i+1].gl_Position << std::endl;  
                // std::cout << "v2: " << geometry_arr[i+2].gl_Position << std::endl;
                clip_triangle(state, 
                              geometry_arr[i],
                              geometry_arr[i+1],
                              geometry_arr[i+2],
                              0);
                
            }
        break;

        case render_type::indexed:
            // std::cout<<"TODO: implement rendering indexed triangles."<<std::endl;
            for(int i = 0; i < state.num_triangles; i++) {
                clip_triangle(state, 
                              geometry_arr[state.index_data[i*state.floats_per_vertex]],
                              geometry_arr[state.index_data[i*state.floats_per_vertex+1]],
                              geometry_arr[state.index_data[i*state.floats_per_vertex+2]],
                              0);
            }
        break;

        case render_type::fan:
            // std::cout<<"TODO: implement rendering fanned triangles."<<std::endl;
            for(int i = 0; i < state.num_vertices-2; i++) {
                clip_triangle(state,
                              geometry_arr[0],
                              geometry_arr[i+1],
                              geometry_arr[i+2],
                              0);
            }

        break;

        case render_type::strip:
            for (int i = 0; i < state.num_vertices - 2; i++) {
                switch(i % 2) {
                    case 0:
                        clip_triangle(state, 
                              geometry_arr[i],
                              geometry_arr[i+2],
                              geometry_arr[i+1],
                              0);
                    break;

                    case 1:
                        clip_triangle(state, 
                              geometry_arr[i],
                              geometry_arr[i+1],
                              geometry_arr[i+2],
                              0);
                    break;

                    default:
                        exit(0);
                    break;
                }
            }
        break;

        default:
            exit(0);
    }
    for(int i = 0; i < state.num_vertices; ++i) {
        delete[] geometry_arr[i].data;
    }
    delete[] geometry_arr;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face){

    float lambda_a_b;
    float lambda_b_c;
    float lambda_c_a;
    short checking_plane;
    float plane_val;

    short a_out;
    short b_out;
    short c_out;
    short out_in_string = 0x00; // 0000 0000 0000 0abc - 1 out 0 in
    
    bool interp_persp = false;
    // // sets up new vertices to pass in if needed
    // data_geometry a;
    // a.gl_Position = v0.gl_Position / v0.gl_Position[3];
    // data_geometry b;
    // b.gl_Position = v1.gl_Position / v1.gl_Position[3];
    // data_geometry c;
    // c.gl_Position = v2.gl_Position / v2.gl_Position[3];
    data_geometry p;
    data_geometry q;
    // a.data = new float[MAX_FLOATS_PER_VERTEX];
    // b.data = new float[MAX_FLOATS_PER_VERTEX];
    // c.data = new float[MAX_FLOATS_PER_VERTEX];
    p.data = new float[MAX_FLOATS_PER_VERTEX];
    q.data = new float[MAX_FLOATS_PER_VERTEX];

    // data_vertex vertex_0;
    // data_vertex vertex_1;
    // data_vertex vertex_2;
    // vertex.data = &(state.vertex_data[i * state.floats_per_vertex]);
    // state.vertex_shader(vertex, geometry_arr[i], state.uniform_data);
    
    if (state.floats_per_vertex > 3) {
        if (state.interp_rules[3] == interp_type::smooth) interp_persp = true;
    }

    switch(face) {
        case 0: //-x
            checking_plane = 0;
            // plane_val = (interp_persp)? v0.gl_Position[3]:-1;
            plane_val = -1;
            // std::cout << "v0: " << v0.gl_Position[0] << ' ' << v0.gl_Position[1] << ' ' << v0.gl_Position[2] << std::endl;  
            // std::cout << "v1: " << v1.gl_Position[0] << ' ' << v1.gl_Position[1] << ' ' << v1.gl_Position[2] << std::endl;  
            // std::cout << "v2: " << v2.gl_Position[0] << ' ' << v2.gl_Position[1] << ' ' << v2.gl_Position[2] << std::endl;  
        break;

        case 1: //+x
            checking_plane = 0;
            // plane_val = (interp_persp)? v0.gl_Position[3]:1;
            plane_val = 1;
        break; 
 
        case 2: //-y
            checking_plane = 1;
            // plane_val = (interp_persp)? v0.gl_Position[3]:-1;
            plane_val = -1;
        break; 
 
        case 3: //+y
            checking_plane = 1;
            // plane_val = (interp_persp)? v0.gl_Position[3]:1;
            plane_val = 1;
        break; 
 
        case 4: //-z
            checking_plane = 2;
            // plane_val = (interp_persp)? v0.gl_Position[3]:-1;
            plane_val = -1;
        break; 
 
        case 5: //+z
            checking_plane = 2;
            // plane_val = (interp_persp)? v0.gl_Position[3]:1;
            plane_val = 1;
        break; 
 
        case 6: //finish
            // data_geometry a;
            // data_geometry b;
            // data_geometry c;
            // a.data = new float[MAX_FLOATS_PER_VERTEX];
            // b.data = new float[MAX_FLOATS_PER_VERTEX];
            // c.data = new float[MAX_FLOATS_PER_VERTEX];
            // a.gl_Position = v0.gl_Position/v0.gl_Position[3];
            // b.gl_Position = v1.gl_Position/v1.gl_Position[3];
            // c.gl_Position = v2.gl_Position/v2.gl_Position[3];
            
            rasterize_triangle(state, v0, v1, v2);
            return;
        break;

        default:
            exit(0);
        break;
    }
    if (plane_val == 1) {
        
        a_out = (v0.gl_Position[checking_plane] > v0.gl_Position[3])? 0x04: 0x00;
        b_out = (v1.gl_Position[checking_plane] > v1.gl_Position[3])? 0x02: 0x00;
        c_out = (v2.gl_Position[checking_plane] > v2.gl_Position[3])? 0x01: 0x00;
    } else {
        a_out = (v0.gl_Position[checking_plane] < -v0.gl_Position[3])? 0x04: 0x00;
        b_out = (v1.gl_Position[checking_plane] < -v1.gl_Position[3])? 0x02: 0x00;
        c_out = (v2.gl_Position[checking_plane] < -v2.gl_Position[3])? 0x01: 0x00;
    }
    out_in_string = a_out | b_out | c_out;
    // std::cout << "face " << face << " out_in: " << out_in_string << std::endl;
    switch(out_in_string){
        case 0x00: // all in
            clip_triangle(state, v0, v1, v2, face+1);
        break;

        case 0x01: // a b in, c out

            lambda_b_c = (plane_val - v1.gl_Position[checking_plane]) / (v2.gl_Position[checking_plane] - v1.gl_Position[checking_plane]);
            lambda_c_a = (plane_val - v0.gl_Position[checking_plane]) / (v2.gl_Position[checking_plane] - v0.gl_Position[checking_plane]);

            p.gl_Position = v0.gl_Position + lambda_c_a * (v2.gl_Position - v0.gl_Position);
            q.gl_Position = v1.gl_Position + lambda_b_c * (v2.gl_Position - v1.gl_Position);

            if (interp_persp){
                float k = 1 / (lambda_c_a * v0.gl_Position[3] + (1 - lambda_c_a) * v2.gl_Position[3]);
                lambda_c_a /= (k * v0.gl_Position[3]);
                k = 1 / (lambda_b_c * v1.gl_Position[3] + (1 - lambda_b_c) * v2.gl_Position[3]);
                lambda_b_c /= (k * v1.gl_Position[3]);
            }
            for (int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_c_a * v0.data[i] + (1-lambda_c_a)*v2.data[i];
                q.data[i] = lambda_b_c * v1.data[i] + (1-lambda_b_c)*v2.data[i];
            }
            clip_triangle(state, v0, v1, p, face+1);
            clip_triangle(state, p, v1, q, face+1);
        break;

        case 0x02: // a c in, b out
            lambda_b_c = (plane_val - v2.gl_Position[checking_plane]) / (v1.gl_Position[checking_plane] - v2.gl_Position[checking_plane]);
            lambda_c_a = (plane_val - v0.gl_Position[checking_plane]) / (v1.gl_Position[checking_plane] - v0.gl_Position[checking_plane]);

            p.gl_Position = v2.gl_Position + lambda_b_c * (v1.gl_Position - v2.gl_Position);
            q.gl_Position = v0.gl_Position + lambda_c_a * (v1.gl_Position - v0.gl_Position);

            if (interp_persp) {
                float k = 1 / (lambda_b_c * v2.gl_Position[3] + (1 - lambda_b_c) * v1.gl_Position[3]);
                lambda_b_c /= (k * v2.gl_Position[3]);
                k = 1 / (lambda_c_a * v0.gl_Position[3] + (1 - lambda_c_a) * v1.gl_Position[3]);
            }
            for (int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_b_c * v2.data[i] + (1-lambda_b_c*v1.data[i]);
                q.data[i] = lambda_c_a * v0.data[i] + (1-lambda_c_a*v1.data[i]);
            }
            clip_triangle(state, v2, v0, p, face+1);
            clip_triangle(state, p, v0, q, face+1);
        break;

        case 0x03: // a in, b c out
            lambda_a_b = (plane_val - v0.gl_Position[checking_plane]) / (v1.gl_Position[checking_plane] - v0.gl_Position[checking_plane]);
            lambda_c_a = (plane_val - v0.gl_Position[checking_plane]) / (v2.gl_Position[checking_plane] - v0.gl_Position[checking_plane]);

            p.gl_Position = v0.gl_Position + lambda_c_a * (v2.gl_Position - v0.gl_Position);
            q.gl_Position = v0.gl_Position + lambda_a_b * (v1.gl_Position - v0.gl_Position);

            if (interp_persp) {
                float k = 1 / (lambda_c_a * v0.gl_Position[3] + (1 - lambda_c_a) * v2.gl_Position[3]);
                lambda_c_a /= (k * v0.gl_Position[3]);
                k = 1 / (lambda_a_b * v0.gl_Position[3] + (1-lambda_a_b) * v1.gl_Position[3]);
                lambda_a_b /= k * v0.gl_Position[3];
            }

            for (int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_c_a * v0.data[i] + (1-lambda_c_a) * v2.data[i];
                q.data[i] = lambda_a_b * v0.data[i] + (1-lambda_a_b) * v1.data[i];                
            }
            //std::cout << "case 3" << std::endl;


            // std::cout << "a_b: " << lambda_a_b << std::endl;
            // std::cout << "c_a: " << lambda_c_a << std::endl;

            // std::cout << "p: " << p.gl_Position[0] << ' ' << p.gl_Position[1] << ' ' << p.gl_Position[2] << std::endl;  
            // std::cout << "q" << q.gl_Position[0] << ' ' << q.gl_Position[1] << ' ' << q.gl_Position[2] << std::endl;
            clip_triangle(state, v0, q, p, face+1);

        break;

        case 0x04: // b c in, \a out
            lambda_a_b = (plane_val - v1.gl_Position[checking_plane]) / (v0.gl_Position[checking_plane] - v1.gl_Position[checking_plane]);
            lambda_c_a = (plane_val - v2.gl_Position[checking_plane]) / (v0.gl_Position[checking_plane] - v2.gl_Position[checking_plane]);

            p.gl_Position = v1.gl_Position + lambda_a_b * (v0.gl_Position - v1.gl_Position);
            q.gl_Position = v2.gl_Position + lambda_c_a * (v0.gl_Position - v2.gl_Position);

            if (interp_persp) {
                float k = 1 / (lambda_a_b * v1.gl_Position[3] + (1-lambda_a_b) * v0.gl_Position[3]);
                lambda_a_b /= k * v1.gl_Position[3];
                k = 1 / (lambda_c_a * v2.gl_Position[3] + (1-lambda_c_a) * v0.gl_Position[3]);
                lambda_c_a /= k * v2.gl_Position[3];
            }
            for(int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_a_b * v1.data[i] + (1-lambda_a_b) * v0.data[i];
                q.data[i] = lambda_c_a * v2.data[i] + (1-lambda_c_a) * v0.data[i];
            }
            // std::cout << "case 4" << std::endl;
            // std::cout << "a_b: " << lambda_a_b << std::endl;
            // std::cout << "c_a: " << lambda_c_a << std::endl;
            // std::cout << "v0: " << v0.gl_Position[0] << ' ' << v0.gl_Position[1] << ' ' << v0.gl_Position[2] << std::endl;  
            // std::cout << "v1: " << v1.gl_Position[0] << ' ' << v1.gl_Position[1] << ' ' << v1.gl_Position[2] << std::endl;  
            // std::cout << "v2: " << v2.gl_Position[0] << ' ' << v2.gl_Position[1] << ' ' << v2.gl_Position[2] << std::endl;  
            // std::cout << "p: " << p.gl_Position[0] << ' ' << p.gl_Position[1] << ' ' << p.gl_Position[2] << std::endl;
            // std::cout << "q" << q.gl_Position[0] << ' ' << q.gl_Position[1] << ' ' << q.gl_Position[2] << std::endl;
            clip_triangle(state, v1, v2, p, face+1);
            clip_triangle(state, p, v2, q, face+1);

        break;

        case 0x05: // b in, a c out
            lambda_a_b = (plane_val - v1.gl_Position[checking_plane]) / (v0.gl_Position[checking_plane] - v1.gl_Position[checking_plane]);
            lambda_b_c = (plane_val - v2.gl_Position[checking_plane]) / (v2.gl_Position[checking_plane] - v1.gl_Position[checking_plane]);

            p.gl_Position = v1.gl_Position + lambda_a_b * (v0.gl_Position - v1.gl_Position);
            q.gl_Position = v1.gl_Position + lambda_b_c * (v2.gl_Position - v1.gl_Position);

            if (interp_persp) {
                float k = 1 / (lambda_a_b * v1.gl_Position[3] + (1-lambda_a_b) * v0.gl_Position[3]);
                lambda_a_b /= k * v1.gl_Position[3];
                k = 1 / (lambda_b_c * v1.gl_Position[3] + (1-lambda_b_c) * v2.gl_Position[3]);
                lambda_c_a /= k * v1.gl_Position[3];
            }
            for(int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_a_b * v1.data[i] + (1-lambda_a_b) * v0.data[i];
                q.data[i] = lambda_b_c * v1.data[i] + (1-lambda_b_c) * v2.data[i];
            }
            // std::cout << "case 5" << std::endl;
            // std::cout << "p: " << p.gl_Position[0] << ' ' << p.gl_Position[1] << ' ' << p.gl_Position[2] << std::endl;
            // std::cout << "q" << q.gl_Position[0] << ' ' << q.gl_Position[1] << ' ' << q.gl_Position[2] << std::endl;
            clip_triangle(state, v1, q, p, face+1);

        break;

        case 0x06: // c in, a b out
            lambda_b_c = (plane_val - v2.gl_Position[checking_plane]) / (v1.gl_Position[checking_plane] - v2.gl_Position[checking_plane]);
            lambda_c_a = (plane_val - v2.gl_Position[checking_plane]) / (v0.gl_Position[checking_plane] - v2.gl_Position[checking_plane]);

            p.gl_Position = v2.gl_Position + lambda_b_c * (v1.gl_Position - v2.gl_Position);
            q.gl_Position = v2.gl_Position + lambda_c_a * (v0.gl_Position - v2.gl_Position);

            if (interp_persp) {
                float k = 1 / (lambda_b_c * v2.gl_Position[3] + (1-lambda_b_c) * v1.gl_Position[3]);
                lambda_b_c /= k * v2.gl_Position[3];
                k = 1 / (lambda_c_a * v2.gl_Position[3] + (1-lambda_c_a) * v0.gl_Position[3]);
                lambda_c_a /= k * v2.gl_Position[3];
            }
            for(int i = 0; i < MAX_FLOATS_PER_VERTEX; ++i) {
                p.data[i] = lambda_b_c * v2.data[i] + (1-lambda_b_c) * v1.data[i];
                q.data[i] = lambda_c_a * v2.data[i] + (1-lambda_c_a) * v0.data[i];
            }

            clip_triangle(state, v2, q, p, face+1);

        break;

        case 0x07: // all out
        break;

        default:
            std::cout << "Inside-outisde for clipping failed\n";
            exit(0);
        break;
    }
    if (p.data != nullptr) {
        delete[] p.data;
        p.data = nullptr;
    }
    if (q.data != nullptr) {
        delete[] q.data;
        q.data = nullptr;        
    }
    // delete[] a.data;
    // delete[] b.data;
    // delete[] c.data;

    // clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    float w_a = v0.gl_Position[3];
    float w_b = v1.gl_Position[3];
    float w_c = v2.gl_Position[3];

    vec3 point_a = vec3(v0.gl_Position[0]/w_a, 
                        v0.gl_Position[1]/w_a, 
                        v0.gl_Position[2]/w_a);
    point_a[0] = ((state.image_width-1)*((point_a[0]+1)/2.0)); 
    point_a[1] = ((state.image_height-1)*((point_a[1]+1)/2.0));

    vec3 point_b = vec3(v1.gl_Position[0]/w_b, 
                        v1.gl_Position[1]/w_b, 
                        v1.gl_Position[2]/w_b);
    point_b[0] = ((state.image_width-1)*((point_b[0]+1)/2.0)); 
    point_b[1] = ((state.image_height-1)*((point_b[1]+1)/2.0));

    vec3 point_c = vec3(v2.gl_Position[0]/w_c, 
                        v2.gl_Position[1]/w_c, 
                        v2.gl_Position[2]/w_c);
    point_c[0] = ((state.image_width-1)*((point_c[0]+1)/2.0));
    point_c[1] = ((state.image_height-1)*((point_c[1]+1)/2.0));

    // can get bounding box of triangle and iterate over that instead of the whole image
    float total_area = 0.5 * ((point_b[0]*point_c[1] - point_c[0]*point_b[1]) + 
                              (point_c[0]*point_a[1]-point_a[0]*point_c[1]) + 
                              (point_a[0]*point_b[1]-point_b[0]*point_a[1]));
    float alpha;
    float beta;
    float gamma;

    // int x_min = std::max(std::min(std::min(point_a[0], point_b[0]), point_c[0]) - 1, float(0));
    // int y_min = std::max(std::min(std::min(point_a[1], point_b[1]), point_c[1]) - 1, float(0));

    // int x_max = std::min(std::max(std::max(point_a[0], point_b[0]), point_c[0]) + 1, float(state.image_width));
    // int y_max = std::min(std::max(std::max(point_a[1], point_b[1]), point_c[1]) + 1, float(state.image_height));

    // for (int i = x_min; i < x_max; ++i) {
    //     for (int j = y_min; j < y_max; ++j) {
    for (int i = 0; i < state.image_width; ++i) {
        for (int j = 0; j < state.image_height; ++j) {
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
                point_p[2] = alpha * point_a[2] + beta * point_b[2] + gamma * point_c[2] + 1;
                data_fragment fragment_color;
                fragment_color.data = new float[MAX_FLOATS_PER_VERTEX];
                data_output pixel_color;
                state.fragment_shader(fragment_color, pixel_color, state.uniform_data);

                int pixel_r = pixel_color.output_color[0]*255;
                int pixel_g = pixel_color.output_color[1]*255;
                int pixel_b = pixel_color.output_color[2]*255;


                if (state.floats_per_vertex > 3) {
                    
                    float k = (alpha / w_a) + (beta / w_b) + (gamma / w_c);

                    switch(state.interp_rules[3]){
                        case interp_type::flat:
                            pixel_r = v0.data[3]*255;
                            pixel_g = v0.data[4]*255;
                            pixel_b = v0.data[5]*255;
                        break;

                        case interp_type::smooth:
                            alpha = alpha / (k * w_a);
                            beta = beta / (k * w_b);
                            gamma = gamma / (k * w_c);

                            pixel_r = v0.data[3]*255*alpha + v1.data[3]*255*beta + v2.data[3]*255*gamma;
                            pixel_g = v0.data[4]*255*alpha + v1.data[4]*255*beta + v2.data[4]*255*gamma;
                            pixel_b = v0.data[5]*255*alpha + v1.data[5]*255*beta + v2.data[5]*255*gamma;
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

                if (point_p[2] < state.image_depth[j*state.image_width+i]) {
                    state.image_color[j*state.image_width+i] = make_pixel(pixel_r, 
                                                                          pixel_g, 
                                                                          pixel_b);
                    state.image_depth[j*state.image_width+i] = point_p[2];
                }
                // std::cout << point_p[2] << std::endl;
                // state.image_color[j*state.image_width+i] = make_pixel(pixel_r, 
                //                                                       pixel_g, 
                //                                                       pixel_b);
                delete[] fragment_color.data;
                // std::cout << "r: " << pixel_color.output_color[0] << " g: " << pixel_color.output_color[1]
                //         << " b: " << pixel_color.output_color[2] << std::endl;

                
                // state.image_color[j*state.image_width+i] = make_pixel(255, 255, 255);
            }
        }
    }
    // std::cout<<"TODO: implement rasterization"<<std::endl;

}

