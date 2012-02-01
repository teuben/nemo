// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//
// Pixel shader for billboarding particles
// 
//
// with ATI hardware, uniform variable MUST be used by output          
// variables. That's why factor_size is used by gl_FrontColor    
//
// !!!!!Attribute variable CAN'T be modified (ex: gl_Color)!!!!!!
//
// ============================================================================
// texture
uniform float alpha;           // alpha color factor                   
uniform float factor_size;     // texture size factor                  
uniform int   use_point;       // false=texture, true=point
uniform int   perspective;     // false=orthographic, true=perspective

// colormap
uniform vec3 colormap[100]; // rgb colormap
uniform int ncmap;          // colormap length
uniform float powalpha;     // alpha power value
uniform int reverse_cmap;   // use reversed color map ?

// physical values
uniform int data_phys_valid;// Is data phys valid ? 
uniform float data_phys_min; // minimum physical value
uniform float data_phys_max; // maximum physical value

// attribute values for all the vertex
attribute float a_sprite_size; // a different value for each particles 
attribute float a_phys_data;   // physical data value for each particles

// varying variable
//varying vec4 col;

// functions declaration
vec4 computeColor();

// ============================================================================
void main()                                                            
{           
  vec4 col=vec4(0.0,0.0,0.0,0.0);
  // compute color
  if (data_phys_valid==1) {
    col=computeColor();
  } else {
    col = vec4(gl_Color.r,gl_Color.g,gl_Color.b,gl_Color.a);   
  }
  
  // compute texture size
  vec4 vert = gl_Vertex;
  if (use_point==1) { // use point = same size whatever the distance from observer
    float pointSize =  a_sprite_size/a_sprite_size*factor_size;
    gl_PointSize = max(2., pointSize*1.);
  } 
  else {           // use texture, size change according to the distance
    float pointSize =  a_sprite_size*factor_size;                          
    vec3 pos_eye = vec3 (gl_ModelViewMatrix * vert);                  
    if (perspective==1) {
      gl_PointSize = max(0.00001, pointSize / (1.0 - pos_eye.z));        
    } else {
      gl_PointSize = max(0.00001, pointSize - pos_eye.z + pos_eye.z);        
    }
  }
  gl_TexCoord[0] = gl_MultiTexCoord0;                                
  gl_Position = ftransform();
//  if (1==0) {
//    gl_FrontColor =  vec4( gl_Color.r+col.x +float(factor_size)*0. + float(use_point)*0.,          
//                           gl_Color.g+col.y                                             ,         
//                           gl_Color.b+col.z                                             ,         
//                           (gl_Color.a+col.w) * alpha);
//  } else {
    gl_FrontColor =  vec4( col.x ,          
                           col.y ,         
                           col.z ,         
                           col.w * alpha);   
  //}
}

// ============================================================================
vec4 computeColor() {
  vec4 col;
  
  if (data_phys_valid==1 && a_phys_data>0.0) {
    float logpri=log(a_phys_data);
    float log_rho=0.0;
            
    if ( (logpri >= data_phys_min) &&
         (logpri <= data_phys_max) &&
         (data_phys_max - data_phys_min) != 0.0) {
      log_rho = (logpri  - data_phys_min) / ( data_phys_max - data_phys_min);
    }
    // use float to avoid too many casting
    float fcindex;
    float fncmap=float(ncmap);
    
    if (reverse_cmap==0) { // normal colormap
      fcindex = log_rho*(fncmap-1.);
    } else {             // reverse colormap
      fcindex = fncmap-1.-(log_rho*(fncmap-1.));
    }
    int cindex;
    cindex=int(min(fcindex,fncmap-1.));
    cindex=int(max(0.,fcindex));
    //cindex=int(max(0.,cindex));
    col.x = colormap[cindex].x;    // red
    col.y = colormap[cindex].y;    // green
    col.z = colormap[cindex].z;    // blue
    if (log_rho>0.0)
      col.w = pow(log_rho,powalpha); // alpha
    else
      col.w = 0.;
  } else {
    col = vec4(gl_Color.r,gl_Color.g,gl_Color.b,1.0);
  }  
  return col;
}
// ============================================================================

