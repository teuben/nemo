/*
N-Body integrator using leapfrog scheme

C99 version, based on C++ version
Jeroen Bédorf
4-Dec-2010

Compile as:

g++ -O4 nbabel.c  -o nbabel
or
gcc -std=c99 nbabel.c -o nbabel -lm


Run as:

more input128 | ./nbabel

*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>


typedef double real;
typedef real real3[3];

//Global variables
int n;          //Number of particles 
real *m;        //Masses array
real3 *r;       //Positions
real3 *v;       //Velocities
real3 *a;       //Accelerations
real3 *a0;      //Prev accelerations


void acceleration()
{
  for(int i=0; i < n; i++)
  {
    a[i][0] = a[i][1] = a[i][2] = 0;  //Reset acceleration
  }
  
  //For each star
  for(int i=0; i < n; i++)
  {
    real rij[3];
    
    //For each remaining star
    //for(int j=i+1; j < n; j++)
    for(int j=0; j < n; j++)
    {
      if(i != j) //prevent self interaction
      {
        //Distance between the two stars
        rij[0] = r[i][0] - r[j][0];
        rij[1] = r[i][1] - r[j][1];
        rij[2] = r[i][2] - r[j][2];
        
        real RdotR = (rij[0]*rij[0]) + (rij[1]*rij[1]) + (rij[2]*rij[2]);
        real apre  = 1.0 / sqrt(RdotR*RdotR*RdotR);
        
        //Update acceleration
        a[i][0] -= m[j]*apre*rij[0];
        a[i][1] -= m[j]*apre*rij[1];
        a[i][2] -= m[j]*apre*rij[2];
      }//end i!=j
    }//end for j      
  }//end for i
}//end acceleration


//Update positions
void updatePositions(real dt)
{
  for(int i=0; i < n; i++)
  {
    //Update the positions, based on the calculated accelerations and velocities
    a0[i][0] = a[i][0]; a0[i][1] = a[i][1]; a0[i][2] = a[i][2];
    
    //for each axis (x/y/z)
    r[i][0] += dt*v[i][0] + 0.5*dt*dt*a0[i][0];
    r[i][1] += dt*v[i][1] + 0.5*dt*dt*a0[i][1];
    r[i][2] += dt*v[i][2] + 0.5*dt*dt*a0[i][2];
  }
}

//Update velocities based on previous and new accelerations
void updateVelocities(real dt)
{
  //Update the velocities based on the previous and old accelerations
  for(int i=0; i < n; i++)
  {     
    v[i][0] += 0.5*dt*(a0[i][0]+a[i][0]);
    v[i][1] += 0.5*dt*(a0[i][1]+a[i][1]);
    v[i][2] += 0.5*dt*(a0[i][2]+a[i][2]);
    
    a0[i][0] = a[i][0]; a0[i][1] = a[i][1]; a0[i][2] = a[i][2];
  }  
}
  
//Compute the energy of the system, 
//contains an expensive O(N^2) part which can be moved to the acceleration part
//where this is already calculated
void energies(real *EKin, real *EPot)
{
  real rij[3];
  real tempKin = 0;
  real tempPot = 0;
  
      
  //Kinetic energy
  for(int i=0; i < n; i++)
  {
    real tempV = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
    tempKin   += 0.5*m[i]*tempV;      
  }
  
  //Potential energy
  for(int i=0; i < n; i++)
  {
    for(int j=i+1; j < n; j++)
    {
      //Distance between the two stars
      rij[0] = r[i][0] - r[j][0];
      rij[1] = r[i][1] - r[j][1];
      rij[2] = r[i][2] - r[j][2];

      real tempRij = (rij[0]*rij[0]) + (rij[1]*rij[1]) + (rij[2]*rij[2]);
      tempPot     -= m[i]*m[j]/sqrt(tempRij);
    }//end for j
  }//end for i
  
  *EKin = tempKin;
  *EPot = tempPot;
}//end energies


int main() {

  int nAlloc = 16384; //Allocated memory size
  
  m  = (real* )malloc(sizeof(real )*nAlloc);
  r  = (real3*)malloc(sizeof(real3)*nAlloc);
  v  = (real3*)malloc(sizeof(real3)*nAlloc);
  a  = (real3*)malloc(sizeof(real3)*nAlloc);
  a0 = (real3*)malloc(sizeof(real3)*nAlloc);
  
  real mass;
  int dummy;
   
  real rx,ry,rz,vx,vy,vz;
  
  int count = 0;
  while(scanf("%d %lf %lf %lf %lf %lf %lf %lf\n", &dummy, &mass, 
        &rx, &ry, &rz, &vx, &vy, &vz) != EOF)
  {
    m[count]    = mass;
    r[count][0] = rx; 
    r[count][1] = ry;
    r[count][2] = rz;
    v[count][0] = vx; 
    v[count][1] = vy;
    v[count][2] = vz;
    //printf("%d %d %f %f %f %f %f %f %f\n", count++, dummy, m, 
    //    rx, ry, rz, vx, vy, vz);
    count++;
    
    //See if we need to allocate more memory
    if(count == nAlloc)
    {
      nAlloc += 1024;
      m       = (real *)realloc(m, sizeof(real )*nAlloc);
      r       = (real3*)realloc(r, sizeof(real3)*nAlloc);
      v       = (real3*)realloc(r, sizeof(real3)*nAlloc);
      a       = (real3*)realloc(r, sizeof(real3)*nAlloc);
      a0      = (real3*)realloc(r, sizeof(real3)*nAlloc);      
    }    
  }
  
  //Set number of items in the system
  n = count;
      
  //Compute initial energu of the system
  real kinEnergy, potEnergy, totEnergy0;
  energies(&kinEnergy, &potEnergy);
  printf("Energies: %lf %lf %lf \n", kinEnergy+potEnergy, kinEnergy, potEnergy);
  totEnergy0 = kinEnergy+potEnergy;
  
  //Start time, end time and simulation step
  real t        = 0.0;
  real tend     = 1.0;
  real dt       = 1e-3;
  int k         = 0;

  //Initialize the accelerations
  acceleration();

  //Start the main loop
  while (t<tend)
  {
    //Update positions based on velocities and accelerations
    updatePositions(dt);
    
    //Get new accelerations
    acceleration();

    //Update velocities
    updateVelocities(dt);
         
    t += dt;
    k += 1;
    if (k%10==0)
    {
      energies(&kinEnergy, &potEnergy);
      real totEnergy = kinEnergy+potEnergy;
      real dE        = (totEnergy-totEnergy0)/totEnergy0;
      printf("t= %lg E=: %lg %lg %lg dE = %lg\n",t, totEnergy, kinEnergy, potEnergy, dE); 
    }    
  } //end while

  return 1;
} //end program

    

