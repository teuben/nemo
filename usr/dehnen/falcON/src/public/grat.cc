//-----------------------------------------------------------------------------+
//                                                                             |
// grat.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// implementing nbdy/grat.h                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/grat.h>
#include <public/Pi.h>
#include <public/nums.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // parameters controlling code                                                
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  const int MAX_TREE_DEPTH = 50;                   // should be plenty          
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary stuff for class grav_mac                                         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //===========================================================================#
  // class nbdy::InvertZ                                                       |
  //                                                                           |
  // methods for inverting                                                     |
  //                                                                           |
  //         theta^(P+2)/(1-theta)^2 * y^a = 1                                 |
  //                                                                           |
  // for  1/theta(y)                                                           |
  //                                                                           |
  //===========================================================================#
  class InvertZ {
  private:
    static const uint N = 1000, N1=N-1;            // size of tables            
    const    unsigned P;                           // expansion order           
    const        real A,hA,sA;                     // parameters                
    real             *Z,*Y;                        // tables                    
    //--------------------------------------------------------------------------
    real z(const real y) const {                   // z(y) = 1/theta - 1        
      if(y < Y[ 0]) return pow(y,hA);
      if(y > Y[N1]) return pow(y,sA);
      return polev(y,Y,Z,N);
    }    
  public:
    //--------------------------------------------------------------------------
    InvertZ(const real a,                          // I: power a                
	    const uint p=3u) :                     //[I: power P]               
      P   ( p ),
      A   ( a ),
      hA  ( half * A ),
      sA  ( A/(P+2.) )
    {
      MemoryCheck(Z = new real[N]); 
      MemoryCheck(Y = new real[N]);
      register double z,iA=1./A,
	zmin = 1.e-4,
	zmax = 1.e4,
	lmin = log(zmin),
	dlz  = (log(zmax)-lmin)/double(N1);
      for(register int i=0; i!=N; ++i) {
	z    = std::exp(lmin+i*dlz);
	Z[i] = z;
	Y[i] = pow(z*z*pow(1+z,P),iA);
      }
    }
    //--------------------------------------------------------------------------
    ~InvertZ() {
      delete[] Z;
      delete[] Y;
    }
    //--------------------------------------------------------------------------
    real invtheta(const real y) const {
      return one + z(y);
    }    
  };
  //============================================================================
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class nbdy::grav_mac                                                         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
grav_mac::grav_mac(const MAC_type mc,
		   const real     t0,
		   const uint     p) :
  MAC  ( mc ),
  P    ( p ),
  TH0  ( min(one,abs(t0)) ),
  iTH0 ( one/TH0)
{
  switch(MAC) {
  case const_theta:
    IZ = 0;
    break;
  case theta_of_M:
    // th^(p+2)    M  (d-2)/d   th0^(p+2)
    // --------  (---)        = ---------
    // (1-th)^2   M0            (1-th0)^2
#if NDIM==2
    IZ = 0;
#else
    MemoryCheck(IZ = new InvertZ(third,P));
#endif
    break;
  case theta_of_M_ov_r:
    // th^(p+2)    Q  (d-2)/(d-1)   th0^(p+2)               M  
    // --------  (---)            = ---------  with  Q := -----
    // (1-th)^2   Q0                (1-th0)^2             r_max
#if NDIM==2
    IZ = 0;
#else
    MemoryCheck(IZ = new InvertZ(half,P));
#endif
    break;
  case theta_of_M_ov_rq:
    // th^(p+2)    S     th0^(p+2)                M   
    // --------  (---) = ---------  with  S := -------
    // (1-th)^2   S0     (1-th0)^2             r_max^2
    MemoryCheck(IZ = new InvertZ(one,P));
    break;
  }
}
//------------------------------------------------------------------------------
void grav_mac::reset(const MAC_type mc,
		     const real     t0, 
		     const uint     p) {
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
  if(MAC != mc || P != p) {
    if(IZ) delete IZ;
    MAC  = mc;
    P    = p;
    switch(MAC) {
    case const_theta:
      IZ = 0;
      break;
    case theta_of_M:
#if NDIM==2
      IZ = 0;
#else
      MemoryCheck(IZ = new InvertZ(third,P));
#endif
      break;
    case theta_of_M_ov_r:
#if NDIM==2
      IZ = 0;
#else
      MemoryCheck(IZ = new InvertZ(half,P));
#endif
      break;
    case theta_of_M_ov_rq:
      MemoryCheck(IZ = new InvertZ(one,P));
      break;
    }
  }
}
//------------------------------------------------------------------------------
void grav_mac::reset_theta(const real t0)
{
  TH0  = min(one,abs(t0));
  iTH0 = one/TH0;
}
//------------------------------------------------------------------------------
void grav_mac::set_rcrit(const grav_tree* T) const {
  switch(MAC) {
  case const_theta:
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(iTH0);
    break;
  case theta_of_M: {
#if NDIM==2
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(iTH0);
#else
    register real 
      M0 = mass(T->root()),
      iF = pow(square(1-TH0)/pow(TH0,P+2u), 3u) / M0;
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(IZ->invtheta(mass(Ci)*iF));
#endif
  } break;
  case theta_of_M_ov_r: {
#if NDIM==2
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(iTH0);
#else
    register int  i  = 0;
    register real Q0 = mass(T->root()) / rmax(T->root());
    register real *Q = new real[T->N_cells()]; MemoryCheck(Q);
    LoopCellsDown(grav_tree,T,Ci) {
      Q[i] = mass(Ci)/rmax(Ci);
      if(Q[i] > Q0) Q0 = Q[i];
      ++i;
    }
    register real iF = square(square(1-TH0)/pow(TH0,P+2u)) / Q0;
    i = 0;
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(IZ->invtheta(iF*Q[i++]));
    delete[] Q;
#endif
  } break;
  case theta_of_M_ov_rq: {
    register int  i  = 0;
    register real S0 = mass(T->root()) / square(rmax(T->root()));
    register real *S = new real[T->N_cells()]; MemoryCheck(S);
    LoopCellsDown(grav_tree,T,Ci) {
      S[i] = mass(Ci)/square(rmax(Ci));
      if(S[i] > S0) S0 = S[i];
      ++i;
    }
    register real iF = square(1-TH0)/pow(TH0,P+2u) / S0;
    i = 0;
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(IZ->invtheta(iF*S[i++]));
    delete[] S;
  } break;
  }
}
//------------------------------------------------------------------------------
grav_mac::~grav_mac()
{
  if(IZ) delete IZ;
}
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary stuff for class grav_tree                                        
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  typedef grav_tree::cell_iterator cell_iter;        // cell iterator           
  typedef grav_tree::soul_iterator soul_iter;        // cell iterator           
  //----------------------------------------------------------------------------
  inline real bmax(const grav_cell*const&C)
    // This routines returns the distance from the cell's cofm to its most      
    // distant corner.                                                          
  {
    register real bq=zero;                         // initialize return variable
    for(register indx i=0; i<NDIM; i++)            // loop dimensions           
      bq += square(radius(C)+abs(cofm(C)[i]-center(C)[i])); // add in square    
    return sqrt(bq);                               // return sqrt of total      
  }
  //----------------------------------------------------------------------------
  template<bool REUSE>
  inline void eval_basic_source(const grav_tree*const&T)
    // This routines computes for all cells of tree T:                          
    //  - mass, 1/mass, cofm, rmax (needed for almost everything)               
    // To be used AFTER souls have been updated                                 
  {
    register real  mon;                            // monopole                  
    register vect  dip;                            // dipole                    
    LoopCellsUp(grav_tree,T,Ci) {                  // loop cells upwards >      
      mon = zero;                                  //   reset monopole          
      dip = zero;                                  //   reset dipole            
      LoopCellKids(cell_iter,Ci,c) {               //   loop sub-cells c >      
	mon += mass(c);                            //       sum up monopole     
	dip += cofm(c);                            //       sum up dipole       
      }                                            //   <                       
      LoopSoulKids(cell_iter,Ci,s) {               //   loop sub-souls s >      
	mon += mass(s);                            //       sum up monopole     
	dip += mass(s)*cofm(s);                    //       sum up dipole       
      }                                            //   <                       
      Ci->mass() = mon;                            //   set mass                
      Ci->imass()= (mon==zero)? zero:one/mon;      //   set 1/mass              
      Ci->cofm() = dip;                            //   set dipole              
    }                                              // <                         
    register real dmax, x, Xq;                     // d_max, auxiliary vars     
    LoopCellsUp(grav_tree,T,Ci) {                  // loop cells upwards >      
      Ci->cofm() *= imass(Ci);                     //   cofm = dipole/mass      
      dmax = zero;                                 //   reset dmax              
      LoopSoulKids(cell_iter,Ci,s)                 //   loop sub-souls s >      
	update_max(dmax,dist_sq(cofm(s),cofm(Ci)));//     update d_max^2 <      
      if(has_soul_kids(Ci)) dmax = sqrt(dmax);     //   d_max due to sub-souls  
      LoopCellKids(cell_iter,Ci,c) {               //   loop sub-cells c >      
	Xq = dist_sq(cofm(c),cofm(Ci));            //     distance^2            
	x  = dmax - rmax(c);                       //     d=distance + r_i      
	if(zero > x || Xq>x*x)                     //     IF(d>d_max)      >    
	  dmax = sqrt(Xq) + rmax(c);               //       set d_max = d  <    
      }                                            //   <                       
      Ci->rmax() = REUSE? dmax:min(dmax,bmax(Ci)); //   assign r_max            
    }                                              // <                         
  }
  //----------------------------------------------------------------------------
  inline void update_and_pass_flags(const grav_tree*const&T,
				    int                  &ns,
				    int                  &nc)
  {
    register int n=0;
    if       (T->use_sbodies()) {                  // 1. loops souls, get flags 
      LoopSouls(grav_tree,T,Si) {
	Si->copy_sink_flag(T->my_sbodies());
	if(is_sink(Si)) ++n;
      }
#ifdef ALLOW_MPI
    } else if(T->use_pbodies()) { 
      LoopSouls(grav_tree,T,Si) {
	Si->copy_sink_flag(T->my_pbodies());
	if(is_sink(Si)) ++n;
      }
#endif
    } else {
      LoopSouls(grav_tree,T,Si) {
	Si->copy_sink_flag(T->my_flags());
	if(is_sink(Si)) ++n;
      }
    }
    ns = n;                                        // # soul sinks              
    n  = 0;
    LoopCellsUp(grav_tree,T,Ci) {                  // 2. loops cells, set flags 
      Ci->reset_sink_flag();
      LoopCellKids(cell_iter,Ci,c) Ci->add_sink_flag(c);
      LoopSoulKids(cell_iter,Ci,s) Ci->add_sink_flag(s);
      if(is_sink(Ci)) n++;
    }
    nc = n;                                        // # cell sinks              
  }
  //----------------------------------------------------------------------------
  template<int ORDER> inline void evaluate_poles (const grav_tree*const&);
  template<int ORDER> inline void normalize_poles(const grav_tree*const&);
  //----------------------------------------------------------------------------
  template<> inline void evaluate_poles<3>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    SYM2(M2);                                      //   macro in tens.h         
    LoopCellsUp(grav_tree,T,Ci) {                  //   loop tree upwards >     
      M2 = zero;                                   //     reset M2 = 0          
      LoopSoulKids(cell_iter,Ci,s) {               //     loop over soul kids > 
	Xi  = cofm(s)-cofm(Ci);                    //       Xi = Yi - Z         
	M2.add_outer_prod(Xi,mass(s));             //       M2 += Xi2 Mi2       
      }                                            //     <                     
      LoopCellKids(cell_iter,Ci,c) {               //     loop over cell kids > 
	M2 += quad(c);                             //       M2 += Mi2           
	Xi  = cofm(c)-cofm(Ci);                    //       Xi = Zi - Z         
	M2.add_outer_prod(Xi,mass(c));             //       M2 += Xi2 Mi2       
      }                                            //     <                     
      Ci->quad() = M2;                             //     assign quadrupole     
    }                                              //   <                       
  }
  template<> inline void normalize_poles<3>(const grav_tree*const&T) {
    LoopCellsDown(grav_tree,T,Ci)                  // loop cells                
      Ci->quad() *= if2*imass(Ci);                 //   normalize quadrupole <  
  }
#if P_ORDER > 3
  template<> inline void evaluate_poles<4>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    SYM2(M2); SYM2(X2);                            //   macro in tens.h         
    SYM3(M3);                                      //   macro in tens.h         
    LoopCellsUp(grav_tree,T,Ci) {                  //   loop tree upwards >     
      M2 = zero;                                   //     reset M2 = 0          
      M3 = zero;                                   //     reset M3 = 0          
      LoopSoulKids(cell_iter,Ci,s) {               //     loop over soul kids > 
	Xi = cofm(s)-cofm(Ci);                     //       Xi  = Yi - Z        
	M2+= X2.outer_prod(Xi,mass(s));            //       M2 += Xi2 Mi        
	M3.add_outer_prod(X2,Xi);                  //       M3 += Xi3 Mi        
      }                                            //     <                     
      LoopCellKids(cell_iter,Ci,c) {               //     loop over cell kids > 
	M2+= quad(c);                              //       M2 += Mi2           
	M3+= octo(c);                              //       M3 += Mi3           
	Xi = cofm(c)-cofm(Ci);                     //       Xi  = Zi - Z        
	M2+= X2.outer_prod(Xi,mass(c));            //       M2 += Xi2 Mi        
	M3.add_outer_prod(X2,Xi);                  //       M3 += Xi3 Mi        
	M3.add_symm_prod(quad(c),Xi);              //       M3 += Xi % Mi2      
      }                                            //     <                     
      Ci->quad() = M2;                             //     assign quadrupole     
      Ci->octo() = M3;                             //     assign octopole       
    }                                              //   <                       
  }
  template<> inline void normalize_poles<4>(const grav_tree*const&T) {
    LoopCellsDown(grav_tree,T,Ci) {                // loop cells                
      Ci->quad() *= if2*imass(Ci);                 //   normalize quadrupole    
      Ci->octo() *= if3*imass(Ci);                 //   normalize octopole      
    }
  }
#if P_ORDER > 4
  template<> inline void evaluate_poles<5>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    SYM2(M2); SYM2(X2);                            //   macro in tens.h         
    SYM3(M3); SYM3(X3);                            //   macro in tens.h         
    SYM4(M4);                                      //   macro in tens.h         
    LoopCellsUp(grav_tree,T,Ci) {                  //   loop tree upwards >     
      M2 = zero;                                   //     reset quadrupole      
      M3 = zero;                                   //     reset octopole        
      M4 = zero;                                   //     reset hexadecupole    
      LoopSoulKids(cell_iter,Ci,s) {               //     loop over soul kids > 
	Xi = cofm(s)-cofm(Ci);                     //       Xi  = Yi - Z        
	M2+= X2.outer_prod(Xi,mass(s));            //       M2 += Xi2 Mi        
	M3+= X3.outer_prod(X2,Xi);                 //       M3 += Xi3 Mi        
	M4.add_outer_prod(X3,Xi);                  //       M4 += Xi4 Mi        
      }                                            //     <                     
      LoopCellKids(cell_iter,Ci,c) {               //     loop over cell kids > 
	M2+= quad(c);                              //       M2 += Mi2           
	M3+= octo(c);                              //       M3 += Mi3           
	M4+= hexa(c);                              //       M4 += Mi4           
	Xi = cofm(c)-cofm(Ci);                     //       Xi  = Zi - Z        
	M2+= X2.outer_prod(Xi,mass(c));            //       M2 += Xi2 Mi        
	M3+= X3.outer_prod(X2,Xi);                 //       M3 += Xi3 Mi        
	M4.add_outer_prod(X3,Xi);                  //       M4 += Xi4 Mi        
	M3+= X3.outer_prod(quad(c),three*Xi);      //       M3 += 3 Xi Mi2      
	M4.add_outer_prod(X3,two*Xi);              //       M4 += 6 Xi2 Mi2     
	M4.add_outer_prod(octo(c),four*Xi);        //       M4 += 4 Xi Mi3      
      }                                            //     <                     
      Ci->quad() = M2;                             //     assign quadrupole     
      Ci->octo() = M3;                             //     assign octopole       
      Ci->hexa() = M4;                             //     assign hexadecupole   
    }                                              //   <                       
  }
  template<> inline void normalize_poles<5>(const grav_tree*const&T) {
    LoopCellsDown(grav_tree,T,Ci) {                // loop cells                
      Ci->quad() *= if2*imass(Ci);                 //   normalize quadrupole    
      Ci->octo() *= if3*imass(Ci);                 //   normalize octopole      
      Ci->hexa() *= if4*imass(Ci);                 //   normalize hexadecupole  
    }
  }
#endif
#endif
}                                                  // END: namespace nbdy       
#ifdef ALLOW_INDI
# include <proper/grat.ind>
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_tree                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_tree::set_cell_source() {
  if(CELL_SOURCE==0) NbdyError("no memory");
  register real*ci=CELL_SOURCE;                    // pter to cell's source     
  LoopMyCellsDown(Ci) {                            // loops cells         >     
    Ci->SOURCE = ci;                               //   give memory to cell     
    ci += nCS;                                     //   pter to next memory     
  }                                                // <                         
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_cell_source() {
  if(CELL_SOURCE) delete[] CELL_SOURCE;            // delete old allocation     
  MemoryCheck(CELL_SOURCE=new real[N_cells()*nCS]);// allocate memory           
  set_cell_source();                               // set cells: pter to memory 
}
//------------------------------------------------------------------------------
inline void grav_tree::set_cell_coeffs() {
  if(CELL_COEFFS==0) NbdyError("no memory");
  register real*ci=CELL_COEFFS;                    // pter to cell's coeffs     
  LoopMyCellsDown(Ci)                              // loops cells         >     
    if(is_sink(Ci)) {                              //   IF(soul is sink)       >
      Ci->COEFFS = ci;                             //     set cell: pter to sink
      ci += nCC;                                   //     pter to next memory   
    } else                                         //   < ELSE                 >
      Ci->COEFFS = 0;
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_cell_coeffs() {
  if(CELL_COEFFS) delete[] CELL_COEFFS;            // delete old allocation     
  MemoryCheck(CELL_COEFFS = new real[Ncs*nCC]);    // allocate memory           
  set_cell_coeffs();                               // set cells: pter to memory 
}
//------------------------------------------------------------------------------
inline void grav_tree::set_soul_sinkpt() {
  if(SOUL_SINKPT==0) NbdyError("no memory");
  register real*si=SOUL_SINKPT;                    // pter to soul's source     
  LoopMySouls(Si)                                  // loops souls              >
    if(is_sink(Si)) {                              //   IF(soul is sink)       >
      Si->SINKPT = si;                             //     set soul: pter to sink
      si += nSS;                                   //     pter to next memory   
    } else                                         //   < ELSE                 >
      Si->SINKPT = 0;
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_soul_sinkpt() {
  if(SOUL_SINKPT) delete[] SOUL_SINKPT;            // delete old allocation     
  MemoryCheck(SOUL_SINKPT = new real[Nss*nSS]);    // allocate memory           
  set_soul_sinkpt();                               // set souls: pter to memory 
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const sbodies*const&b,        // I: sbodies                
#ifdef ALLOW_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,n,MAX_TREE_DEPTH ),
#ifdef ALLOW_INDI
  SOFT        ( s ),
#endif
  M           ( 0 ),
#ifdef ALLOW_INDI
  EP          ( 0 ),
#endif
  nCS         ( 
#ifdef ALLOW_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  SOUL_SINKPT ( 0 ),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 )
{
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const sbodies*const&b,        // I: sbodies                
		     vect          const&xmin,     // I: x_min                  
		     vect          const&xmax,     // I: x_max                  
#ifdef ALLOW_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,xmin,xmax,n,MAX_TREE_DEPTH ),
#ifdef ALLOW_INDI
  SOFT        ( s ),
#endif
  M           ( 0 ),
#ifdef ALLOW_INDI
  EP          ( 0 ),
#endif
  nCS         ( 
#ifdef ALLOW_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  SOUL_SINKPT ( 0 ),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 )
{
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
#ifdef ALLOW_MPI
//------------------------------------------------------------------------------
grav_tree::grav_tree(const pbodies*const&b,        // I: pbodies                
#ifdef ALLOW_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,n,MAX_TREE_DEPTH ),
#ifdef ALLOW_INDI
  SOFT        ( s ),
#endif
  M           ( 0 ),
#ifdef ALLOW_INDI
  EP          ( 0 ),
#endif
  nCS         ( 
#ifdef ALLOW_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  SOUL_SINKPT ( 0 ),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 )
{
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const pbodies*const&b,        // I: pbodies                
		     vect          const&xmin,     // I: x_min                  
		     vect          const&xmax,     // I: x_max                  
#ifdef ALLOW_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,xmin,xmax,n,MAX_TREE_DEPTH ),
#ifdef ALLOW_INDI
  SOFT        ( s ),
#endif
  M           ( 0 ),
#ifdef ALLOW_INDI
  EP          ( 0 ),
#endif
  nCS         ( 
#ifdef ALLOW_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  SOUL_SINKPT ( 0 ),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 )
{
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
#endif
//------------------------------------------------------------------------------
grav_tree::grav_tree(const int      *f,            // I: array with flags       
		     const areal    *x[NDIM],      // I: arrays with x,y,z      
		     const areal    *m,            // I: array with masses      
#ifdef ALLOW_INDI
		     const areal    *ep,           // I: array with eps_i       
#endif
		     const uint      nb,           // I: size of arrays         
#ifdef ALLOW_INDI
		     const soft_type s,            //[I: global/individual]     
#endif
		     const int       n) :          //[I: N_crit]                
  base_tree   ( f,x,nb,n,MAX_TREE_DEPTH ),
#ifdef ALLOW_INDI
  SOFT        ( s ),
#endif
  M           ( m ),
#ifdef ALLOW_INDI
  EP          ( ep ),
#endif
  nCS         ( 
#ifdef ALLOW_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  SOUL_SINKPT ( 0 ),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 )
{
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(M);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
void grav_tree::rebuild(const int       Nc,        //[I: N_crit]                
			const int       Nu) {      //[I: N_cut for re_grow]     
  if(SOUL_SINKPT) { delete[] SOUL_SINKPT; SOUL_SINKPT=0; } // free memory       
  if(CELL_SOURCE) { delete[] CELL_SOURCE; CELL_SOURCE=0; } // free memory       
  if(Nu) base_tree::rebuild(Nu,Nc,MAX_TREE_DEPTH); // rebuild base tree   OR    
  else   base_tree::build  (   Nc,MAX_TREE_DEPTH); // build   base tree         
  if     (use_sbodies()) LoopMySouls(Si) Si->set_mass(my_sbodies());
#ifdef ALLOW_MPI
  else if(use_pbodies()) LoopMySouls(Si) Si->set_mass(my_pbodies());
#endif
  else                   LoopMySouls(Si) Si->set_mass(M);
  reset_cell_source();                             // set memory: cell's source 
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
void grav_tree::reuse() {
  if(SOUL_SINKPT) { delete[] SOUL_SINKPT; SOUL_SINKPT=0; } // free memory       
  if     (use_sbodies()) LoopMySouls(Si) Si->set_mass_and_pos(my_sbodies());
#ifdef ALLOW_MPI
  else if(use_pbodies()) LoopMySouls(Si) Si->set_mass_and_pos(my_pbodies());
#endif
  else                   LoopMySouls(Si) Si->set_mass_and_pos(M,my_pos());
  eval_basic_source<true>(this);
}
//------------------------------------------------------------------------------
void grav_tree::prepare_density(bool re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update & pass; count sinks
  if(re_use_mem) set_soul_sinkpt();                // give sinkpt to sink souls 
  else         reset_soul_sinkpt();                // give sinkpt to sink souls 
}
//------------------------------------------------------------------------------
void grav_tree::prepare_grav_exact(
#ifdef ALLOW_INDI
				   real Ns, real ex, uint Nr, real fe,
#endif
				   bool re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update & pass; count sinks
#ifdef ALLOW_INDI
  if(SOFT==individual)                             // IF individual softening   
    update_adjust_and_pass_eph(this,Ns,ex,Nr,fe);  //   get eph_i = eps_i/2     
#endif
  if(re_use_mem) set_soul_sinkpt();                // give sinkpt to sink souls 
  else         reset_soul_sinkpt();                // give sinkpt to sink souls 
  LoopMySouls(Si) if(is_sink(Si)) Si->reset_srce();// reset souls' grav source  
}
//------------------------------------------------------------------------------
void grav_tree::prepare_grav_approx(const grav_mac*MAC,
				    bool           give_coeffs,
#ifdef ALLOW_INDI
				    real           Ns, 
				    real           ex,
				    uint           Nr,
				    real           fe,
#endif
				    bool           re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update & pass; count sinks
#ifdef ALLOW_INDI
  if(SOFT==individual)                             // IF individual softening   
    update_adjust_and_pass_eph(this,Ns,ex,Nr,fe);  //   get eph_i = eps_i/2     
#endif
  if(re_use_mem) set_soul_sinkpt();                // give sinkpt to sink souls 
  else         reset_soul_sinkpt();                // give sinkpt to sink souls 
  LoopMySouls(Si) if(is_sink(Si)) Si->reset_srce();// reset souls' grav source  
  if(give_coeffs)                                  // IF(cell coeffs to give)  >
    if(re_use_mem) set_cell_coeffs();              //  give coeffs to sink cells
    else         reset_cell_coeffs();              //  give coeffs to sink cells
  evaluate_poles <P_ORDER>(this);                  // compute the multipoles    
  normalize_poles<P_ORDER>(this);                  // normalize poles           
  MAC->set_rcrit(this);                            // set r_crit for all cells  
}
//------------------------------------------------------------------------------
void grav_tree::prepare_neighbour_counting(
					   const real* SIZE,
					   bool re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update & pass; count sinks
  if(re_use_mem) set_soul_sinkpt();                // give sinkpt to sink souls 
  else         reset_soul_sinkpt();                // give sinkpt to sink souls 
#ifdef ALLOW_INDI
  if(SOFT == individual) {                         // IF(individual sizes)     >
    if       (use_sbodies()) {
      LoopMySouls(Si) if(is_sink(Si)) {
	Si->copy_eph(my_sbodies());                //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
#ifdef ALLOW_MPI
    } else if(use_pbodies()) { 
      LoopMySouls(Si) if(is_sink(Si)) {
	Si->copy_eph(my_pbodies());                //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
#endif
    } else {
      LoopMySouls(Si) if(is_sink(Si)) {
	Si->copy_eph(my_eps());                    //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
    }
  } else                                           // < ELSE(global size)      >
#endif
    LoopMySouls(Si) if(is_sink(Si)) {              //   loop sink souls        >
      Si->num()   = 0u;                            //     reset counter         
    }                                              //   <                       
  register real smax;                              // cell size                 
  LoopMyCellsUp(Ci) {                              // loop cells upwards       >
    smax = zero;                                   //   reset cell size         
#ifdef ALLOW_INDI
    if(SIZE == 0) {
      LoopSoulKids(cell_iter,Ci,s) if(is_sink(s)){ //   loop soul kids         >
	update_max(smax, sqrt(dist_sq(pos(Ci),pos(s))) + size(s));//update size 
      }                                            //   <                       
    } else {
#endif
      LoopSoulKids(cell_iter,Ci,s) if(is_sink(s)){ //   loop soul kids         >
	update_max(smax, sqrt(dist_sq(pos(Ci),pos(s))) + *SIZE); // update size 
      }                                            //   <                       
#ifdef ALLOW_INDI
    }
#endif
    LoopCellKids(cell_iter,Ci,c) if(is_sink(c)){   //   loop cell kids         >
      update_max(smax, sqrt(dist_sq(pos(Ci),pos(c))) + size(c)); // update size 
    }                                              //   <                       
    Ci->size() = smax;                             //   set cell size           
  }                                                // <                         
}
////////////////////////////////////////////////////////////////////////////////
