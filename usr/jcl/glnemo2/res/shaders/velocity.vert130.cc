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
#version 130

uniform mat4 modelviewMatrix;
uniform mat4 projMatrix;

uniform float vel_factor;

in vec3 position;
in vec3 velocity;

// ============================================================================
void main()                                                            
{
  //vec4 o_position = vec4(position.x* vel_factor,position.y* vel_factor,position.z* vel_factor, 1.0f);
  vec4 o_position = vec4(position+velocity*vel_factor, 1.0);
  gl_Position = projMatrix*modelviewMatrix * o_position;

}

// ============================================================================
