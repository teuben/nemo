// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
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
// ============================================================================
// texture
uniform float alpha;           // alpha color factor                   
uniform float factor_size;     // texture size factor                  
uniform int   use_point;       // false=texture, true=point

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

// functions declaration
vec4 computeColor();

// ============================================================================
void main()                                                            
{                                                                      
  // compute color
  if (data_phys_valid) {
    vec4 col=computeColor();
    gl_Color.r = col.x;
    gl_Color.g = col.y;
    gl_Color.b = col.z;
    gl_Color.a = col.w;
  }
  
  // compute texture size
  vec4 vert = gl_Vertex;
  if (use_point) { // use point = same size whatever the distance from observer
    float pointSize =  a_sprite_size/a_sprite_size*factor_size;
    gl_PointSize = max(2., pointSize*1.);
  } 
  else {           // use texture, size change according to the distance
    float pointSize =  a_sprite_size*factor_size;                          
    vec3 pos_eye = vec3 (gl_ModelViewMatrix * vert);                   
    gl_PointSize = max(0.00001, pointSize / (1.0 - pos_eye.z));        
  }
  gl_TexCoord[0] = gl_MultiTexCoord0;                                
  gl_Position = ftransform();                                        
  gl_FrontColor =  vec4(gl_Color.r + float(factor_size)*0. + use_point*0.,          
                        gl_Color.g                                       ,         
                        gl_Color.b                                       ,         
                        gl_Color.a * alpha);
}

// ============================================================================
vec4 computeColor() {
  vec4 col;
  
  if (data_phys_valid && a_phys_data>0.0) {
    float logpri=log(a_phys_data);
    float log_rho=0.0;
    
    int cindex;
    
    if ( (logpri >= data_phys_min) &&
         (logpri <= data_phys_max) &&
         (data_phys_max - data_phys_min) != 0.0) {
      log_rho = (logpri  - data_phys_min) / ( data_phys_max - data_phys_min);
    }
    
    if (!reverse_cmap) { // normal colormap
      cindex = log_rho*(ncmap-1);
    } else {             // reverse colormap
      cindex = ncmap-1-log_rho*(ncmap-1);
    }
    cindex=min(cindex,ncmap-1);
    cindex=max(0,cindex);
    col.x = colormap[cindex].x;    // red
    col.y = colormap[cindex].y;    // green
    col.z = colormap[cindex].z;    // blue
    if (log_rho>0)
      col.w = pow(log_rho,powalpha); // alpha
    else
      col.w = 0.;
  } else {
    col = (gl_Color.r,gl_Color.g,gl_Color.b,1.0);
    //col = (1.0,1.0,1.0,1.0);
  }  
  return col;
}
// ============================================================================

