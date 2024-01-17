/*
 * YAPP: Yet Another Plotting Package.
 *
 *  Using the OpenGL gflw interface.  See https://www.glfw.org/
 *
 *  Ubuntu:     sudo apt-get install libglfw3-dev
 *  microbrew:  brew tap homebrew/versions && brew install glfw3
 */

#include <stdinc.h>
#include <yapp.h>
#include <GLFW/glfw3.h>

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
 dprintf(0,"[YAPP_GLFW (experimental)\n");

 int major, minor, revision;
 glfwGetVersion(&major, &minor, &revision);
 dprintf(0,"Running against GLFW %i.%i.%i\n", major, minor, revision);
 
 GLFWwindow* window;
 if (!glfwInit())
   return -1;
 window = glfwCreateWindow(1000, 1000, "yapp_glfw3", NULL, NULL);
 if (!window) {
   glfwTerminate();
   return -1;
 }
 glfwMakeContextCurrent(window);

#if 1 
 while (!glfwWindowShouldClose(window)) {
   /* Render here */
   //glClear(GL_COLOR_BUFFER_BIT);
   /* Swap front and back buffers */
   glfwSwapBuffers(window);
   /* Poll for and process events */
   glfwPollEvents();
 }
 glfwTerminate();
#endif 
 return 0;
}

int plswap() 
{ return 0;}

real plxscale(real x, real y)
{ return 0;}

real plyscale(real x, real y)
{ return 0;}

int plltype(int lwid, int lpat)
{ return 0;}

int plline(real x, real y)
{ return 0;}

int plmove(real x, real y)
{ return 0;}

int plpoint(real x, real y)
{ return 0;}

int plcircle(real x, real y, real r)
{ return 0;}

int plcross(real x, real y, real s)
{ return 0;}

int plbox(real x, real y, real s)
{ return 0;}

int pljust(int jus)
{ return 0;}

int pltext(string msg, real x, real y, real hgt, real ang)
{ return 0;}

int plflush() 
{ return 0;}

int plframe()
{ return 0;}

int plstop()
{ return 0;}

void plcolor(int color)
{ }

int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin,
	      real cell,real fmin,real fmax,real findex, real blankval)
{ return 0;}

int pl_contour(real *frame,int nx,int ny, int nc, real *c)
{ return 0;}

int pl_screendump(string fname)
{ return 0;}

int pl_getpoly(float *x,float *y,int n)
{ return 0;}

int pl_cursor(real *x,real *y, char *c)
{ return 0;}

int pl_interp(string cmd)
{ return 0;}
