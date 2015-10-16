// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//
// Geometry shader to display velocity vectors
// 
// we draw velocity vector at the corresponding position
//
// ============================================================================
#version 330 core
#extension GL_ARB_explicit_uniform_location : enable

//uniform float alpha;           // alpha color factor
uniform mat4 modelviewMatrix;
uniform mat4 projMatrix;

in vec4 g_velocity[]; // velocity vector for each particles
uniform float vel_factor;

layout (points) in;
layout (line_strip, max_vertices = 2) out;

in VS_OUT {
    vec4 color;
} gs_in[];

out vec4 fColor;
// ============================================================================
void buildVelocityVectors(vec4 position)
{           
    fColor = gs_in[0].color; // gs_in[0] since there's only one input vertex

    // draw first point
    mat4 mvp=projMatrix*modelviewMatrix;  // the model view projection matrix
    gl_Position = mvp * position; // at the vertex position
    EmitVertex();

    // draw second point
    // which represent veolicty vector at the selected position
    // pos+vel

    vec4 vel =  g_velocity[0]*vel_factor;

    // It's important to keep 4th vertex coordinates to 1.0, otherwise
    // result is totally wrong !!!
    gl_Position = mvp * vec4(position.xyz+vel.xyz,1.0f) ; // at vertex position + vel vector

    EmitVertex();

    EndPrimitive();
}


void main()
{
  buildVelocityVectors(gl_in[0].gl_Position);
}

// ============================================================================
