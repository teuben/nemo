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
// Pixel shader to display velocity vectors
// 
//
// with ATI hardware, uniform variable MUST be used by output          
// variables. That's why factor_size is used by gl_FrontColor    
//
// !!!!!Attribute variable CAN'T be modified (ex: gl_Color)!!!!!!
//
// ============================================================================
#version 330 core


uniform float alpha;           // alpha color factor
uniform vec3 color;
//uniform float vel_factor;

in vec3 position;
in vec3 velocity; // velocity vector for each particles

out vec4 g_velocity;

out VS_OUT {
    vec4 color;
} vs_out;

// ============================================================================
void main()                                                            
{
  //color=vec3(0.0f,0.5f,1.0f);
  //alpha=1.0f;
  gl_Position = vec4(position.x, position.y, position.z,  1.0f);
  g_velocity = vec4(velocity.x,velocity.y,velocity.z, 1.0f);//*vel_factor;
  vs_out.color =  vec4( color.x ,color.y ,color.z , alpha);
}

// ============================================================================
