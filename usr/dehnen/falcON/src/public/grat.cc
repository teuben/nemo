//-----------------------------------------------------------------------------+
//                                                                             |
// grat.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
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
#include <public/iact.h>
#include <public/grav.h>
#include <public/Pi.h>
#include <public/nums.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
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
      sA  ( A/(P+2.) ),
      Z   ( falcON_New(real,N) ),
      Y   ( falcON_New(real,N) )
    {
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
#if falcON_NDIM==2
    IZ = 0;
#else
    IZ = falcON_Memory(new InvertZ(third,P));
#endif
    break;
  case theta_of_M_ov_r:
    // th^(p+2)    Q  (d-2)/(d-1)   th0^(p+2)               M  
    // --------  (---)            = ---------  with  Q := -----
    // (1-th)^2   Q0                (1-th0)^2             r_max
#if falcON_NDIM==2
    IZ = 0;
#else
    IZ = falcON_Memory(new InvertZ(half,P));
#endif
    break;
  case theta_of_M_ov_rq:
    // th^(p+2)    S     th0^(p+2)                M   
    // --------  (---) = ---------  with  S := -------
    // (1-th)^2   S0     (1-th0)^2             r_max^2
    IZ = falcON_Memory(new InvertZ(one,P));
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
#if falcON_NDIM==2
      IZ = 0;
#else
      IZ = falcON_Memory(new InvertZ(third,P));
#endif
      break;
    case theta_of_M_ov_r:
#if falcON_NDIM==2
      IZ = 0;
#else
      IZ = falcON_Memory(new InvertZ(half,P));
#endif
      break;
    case theta_of_M_ov_rq:
      IZ = falcON_Memory(new InvertZ(one,P));
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
#if falcON_NDIM==2
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(iTH0);
#else
    register real 
      M0 = mass(T->root()),
      iF = pow(square(1-TH0)/pow(TH0,P+2u), 3u) / M0;
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(IZ->invtheta(mass(Ci)*iF));
#endif
  } break;
  case theta_of_M_ov_r: {
#if falcON_NDIM==2
    LoopCellsDown(grav_tree,T,Ci) Ci->set_rcrit(iTH0);
#else
    register int  i  = 0;
    register real Q0 = mass(T->root()) / rmax(T->root());
    register real *Q = falcON_New(real,T->N_cells());
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
    register real *S = falcON_New(real,T->N_cells());
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
    for(register int d=0; d!=Ndim; ++d)            // loop dimensions           
      bq += square(radius(C)+abs(cofm(C)[d]-center(C)[d])); // add in square    
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
	Si->copy_active_flag(T->my_sbodies());
	if(is_active(Si)) ++n;
      }
#ifdef falcON_MPI
    } else if(T->use_pbodies()) { 
      LoopSouls(grav_tree,T,Si) {
	Si->copy_active_flag(T->my_pbodies());
	if(is_active(Si)) ++n;
      }
#endif
    } else if(T->use_barrays()) {
      LoopSouls(grav_tree,T,Si) {
	Si->copy_active_flag(T->my_barrays());
	if(is_active(Si)) ++n;
      }
    } else
      falcON_Error("tree has neither bodies nor array data");
    ns = n;                                        // # active souls            
    n  = 0;
    LoopCellsUp(grav_tree,T,Ci) {                  // 2. loops cells, set flags 
      Ci->reset_active_flag();
      LoopCellKids(cell_iter,Ci,c) Ci->add_active_flag(c);
      LoopSoulKids(cell_iter,Ci,s) Ci->add_active_flag(s);
      if(is_active(Ci)) n++;
    }
    nc = n;                                        // # active cells            
  }
  //----------------------------------------------------------------------------
  template<int ORDER> inline void evaluate_poles (const grav_tree*const&);
  template<int ORDER> inline void normalize_poles(const grav_tree*const&);
  //----------------------------------------------------------------------------
  template<> inline void evaluate_poles<3>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    falcON_SYM2(M2);                               //   macro in tens.h         
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
#if falcON_ORDER > 3
  template<> inline void evaluate_poles<4>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    falcON_SYM2(M2); falcON_SYM2(X2);              //   macro in tens.h         
    falcON_SYM3(M3);                               //   macro in tens.h         
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
#if falcON_ORDER > 4
  template<> inline void evaluate_poles<5>(const grav_tree*const&T) {
    register vect Xi;                              // distance vector           
    falcON_SYM2(M2); falcON_SYM2(X2);              //   macro in tens.h         
    falcON_SYM3(M3); falcON_SYM3(X3);              //   macro in tens.h         
    falcON_SYM4(M4);                               //   macro in tens.h         
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
#ifdef falcON_INDI
# include <proper/grat_ind.cc>
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::grav_tree                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void grav_tree::set_cell_source() {
  if(CELL_SOURCE==0) falcON_Error("no memory");
  register real*ci=CELL_SOURCE;                    // pter to cell's source     
  LoopMyCellsDown(Ci) {                            // loops cells         >     
    Ci->SOURCE = ci;                               //   give memory to cell     
    ci += nCS;                                     //   pter to next memory     
  }                                                // <                         
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_cell_source() {
  if(CELL_SOURCE) delete[] CELL_SOURCE;            // delete old allocation     
  CELL_SOURCE = falcON_New(real,N_cells()*nCS);    // allocate memory           
  set_cell_source();                               // set cells: pter to memory 
}
//------------------------------------------------------------------------------
inline void grav_tree::set_cell_coeffs(bool const&all) {
  if(CELL_COEFFS==0) falcON_Error("no memory");
  register real*ci=CELL_COEFFS;                    // pter to cell's coeffs     
  if(all)                                          // IF all souls              
    LoopMyCellsDown(Ci) {                          //   LOOP cells              
      Ci->COEFFS = ci;                             //     set cell: pter->sink  
      ci += nCC;                                   //     pter to next memory   
    }                                              //   END LOOP                
  else                                             // ELSE (only active)        
    LoopMyCellsDown(Ci)                            //   LOOP cells              
      if(is_active(Ci)) {                          //     IF(soul is active)    
	Ci->COEFFS = ci;                           //       set cell: pter->sink
	ci += nCC;                                 //       pter to next memory 
      } else                                       //     ELSE                  
	Ci->COEFFS = 0;                            //       cell: pter->sink = 0
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_cell_coeffs(bool const&all) {
  if(CELL_COEFFS) delete[] CELL_COEFFS;            // delete old allocation     
  CELL_COEFFS = falcON_New(real,Ncs*nCC);          // allocate memory           
  set_cell_coeffs(all);                            // set cells: pter to memory 
}
//------------------------------------------------------------------------------
inline void grav_tree::set_soul_sinkpt(bool const&all) {
  if(SOUL_SINKPT==0) falcON_Error("no memory");
  register real*si=SOUL_SINKPT;                    // pter to soul's source     
  if(all)                                          // IF all souls              
    LoopMySouls(Si) {                              //   LOOP souls              
      Si->SINKPT = si;                             //     set soul: pter->sink  
      si += nSS;                                   //     pter to next memory   
    }                                              //   END LOOP                
  else                                             // ELSE (only active)        
    LoopMySouls(Si)                                //   LOOP souls              
      if(is_active(Si)) {                          //     IF(soul is active)    
	Si->SINKPT = si;                           //       set soul: pter->sink
	si += nSS;                                 //       pter to next memory 
      } else                                       //     ELSE                  
	Si->SINKPT = 0;                            //       soul: pter->sink = 0
}
//------------------------------------------------------------------------------
inline void grav_tree::reset_soul_sinkpt(bool const&all) {
  if(SOUL_SINKPT) delete[] SOUL_SINKPT;            // delete old allocation     
  SOUL_SINKPT = falcON_New(real,Nss*nSS);          // allocate memory           
  set_soul_sinkpt(all);                            // set souls: pter to memory 
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const sbodies*const&b,        // I: sbodies                
#ifdef falcON_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,n,Default::MaxDepth ),
#ifdef falcON_INDI
  SOFT        ( s ),
#endif
  nCS         ( 
#ifdef falcON_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  Ncoeffs     ( 0u),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 ),
  SOUL_SINKPT ( 0 )
{
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const sbodies*const&b,        // I: sbodies                
		     vect          const&xmin,     // I: x_min                  
		     vect          const&xmax,     // I: x_max                  
#ifdef falcON_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,xmin,xmax,n,Default::MaxDepth ),
#ifdef falcON_INDI
  SOFT        ( s ),
#endif
  nCS         ( 
#ifdef falcON_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  Ncoeffs     ( 0u),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 ),
  SOUL_SINKPT ( 0 )
{
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
#ifdef falcON_MPI
//------------------------------------------------------------------------------
grav_tree::grav_tree(const pbodies*const&b,        // I: pbodies                
#ifdef falcON_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,n,Default::MaxDepth ),
#ifdef falcON_INDI
  SOFT        ( s ),
#endif
  nCS         ( 
#ifdef falcON_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  Ncoeffs     ( 0u),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 ),
  SOUL_SINKPT ( 0 )
{
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
grav_tree::grav_tree(const pbodies*const&b,        // I: pbodies                
		     vect          const&xmin,     // I: x_min                  
		     vect          const&xmax,     // I: x_max                  
#ifdef falcON_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,xmin,xmax,n,Default::MaxDepth ),
#ifdef falcON_INDI
  SOFT        ( s ),
#endif
  nCS         ( 
#ifdef falcON_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  Ncoeffs     ( 0u),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 ),
  SOUL_SINKPT ( 0 )
{
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
#endif
//------------------------------------------------------------------------------
grav_tree::grav_tree(const barrays*const&b,        // I: body arrays            
#ifdef falcON_INDI
		     const soft_type     s,        //[I: global/individual]     
#endif
		     const int           n) :      //[I: N_crit]                
  base_tree   ( b,n,Default::MaxDepth ),
#ifdef falcON_INDI
  SOFT        ( s ),
#endif
  nCS         ( 
#ifdef falcON_INDI
	        SOFT==global? grav_cell::N_eph() : 
#endif
	        grav_cell::N_tot() ),
  Ncoeffs     ( 0u),
  CELL_SOURCE ( 0 ),
  CELL_COEFFS ( 0 ),
  SOUL_SINKPT ( 0 )
{
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  reset_cell_source();                             // set cell's memory         
  LoopMySouls(Si) Si->set_mass(b);                 // set masses                
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
void grav_tree::rebuild(const int       Nc,        //[I: N_crit]                
			const int       Nu) {      //[I: N_cut for re_grow]     
  report REPORT("grav_tree::rebuild(%d,%d)",Nc,Nu);
  if(SOUL_SINKPT) { delete[] SOUL_SINKPT; SOUL_SINKPT=0; } // free memory       
  if(CELL_SOURCE) { delete[] CELL_SOURCE; CELL_SOURCE=0; } // free memory       
  if(Nu) base_tree::rebuild(Nu,Nc,Default::MaxDepth); // rebuild base tree  OR  
  else   base_tree::build  (   Nc,Default::MaxDepth); // build   base tree      
  LoopMyCellsDown(Ci) Ci->COEFFS = 0;              // reset cell coeffs         
  if     (use_sbodies()) LoopMySouls(Si) Si->set_mass(my_sbodies());
#ifdef falcON_MPI
  else if(use_pbodies()) LoopMySouls(Si) Si->set_mass(my_pbodies());
#endif
  else if(use_barrays()) LoopMySouls(Si) Si->set_mass(my_barrays());
  else falcON_Error("tree has neither bodies nor arrays for data");
  reset_cell_source();                             // set memory: cell's source 
  eval_basic_source<false>(this);                  // set basic source props    
}
//------------------------------------------------------------------------------
void grav_tree::reuse() {
  report REPORT("grav_tree::reuse()");
  if(SOUL_SINKPT) { delete[] SOUL_SINKPT; SOUL_SINKPT=0; } // free memory       
  if     (use_sbodies()) LoopMySouls(Si) Si->set_mass_and_pos(my_sbodies());
#ifdef falcON_MPI
  else if(use_pbodies()) LoopMySouls(Si) Si->set_mass_and_pos(my_pbodies());
#endif
  else if(use_barrays()) LoopMySouls(Si) Si->set_mass_and_pos(my_barrays());
  else falcON_Error("tree has neigher bodies nor arrays for data");
  eval_basic_source<true>(this);
}
//------------------------------------------------------------------------------
void grav_tree::prepare_density(bool re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update &pass; count active
  if(re_use_mem) set_soul_sinkpt();                // active souls: give sinkpt 
  else         reset_soul_sinkpt();                // active souls: give sinkpt 
}
//--------------------------------------------------------------------------
inline void grav_tree::update_grav_eps(bool const&U,
				       bool const&all)
{
  if(use_sbodies())
#ifdef falcON_INDI
    if(U)
      if(all)
	LoopMySouls(Si) {
	  Si->update_grav(my_sbodies());
	  Si->update_eps (my_sbodies());
        }
      else
	LoopMySouls(Si) { if(is_active(Si)) {
	  Si->update_grav(my_sbodies());
	  Si->update_eps (my_sbodies());
	} }
    else
#endif
      if(all)
	LoopMySouls(Si) Si->update_grav(my_sbodies());
      else
        LoopMySouls(Si) { if(is_active(Si)) 
	  Si->update_grav(my_sbodies());
        }
#ifdef falcON_MPI
  else if(use_pbodies())
#ifdef falcON_INDI
    if(U)
      if(all)
	LoopMySouls(Si) {
	  Si->update_grav(my_pbodies());
	  Si->update_eps (my_pbodies());
        }
      else
	LoopMySouls(Si) { if(is_active(Si)) {
	  Si->update_grav(my_pbodies());
	  Si->update_eps (my_pbodies());
	} }
    else
#endif
      if(all)
	LoopMySouls(Si) Si->update_grav(my_pbodies());
      else
        LoopMySouls(Si) { if(is_active(Si)) 
	  Si->update_grav(my_pbodies());
        }
#endif
  else if(use_barrays())
#ifdef falcON_INDI
    if(U)
      if(all)
	LoopMySouls(Si) {
	  Si->update_grav(my_barrays());
	  Si->update_eps (my_barrays());
        }
      else
	LoopMySouls(Si) { if(is_active(Si)) {
	  Si->update_grav(my_barrays());
	  Si->update_eps (my_barrays());
	} }
    else
#endif
      if(all)
	LoopMySouls(Si) Si->update_grav(my_barrays());
      else
        LoopMySouls(Si) { if(is_active(Si)) 
	  Si->update_grav(my_barrays());
        }
  else
    falcON_Error("tree has neither bodies nor arrays for data");
}
//------------------------------------------------------------------------------
inline void grav_tree::prepare_grav_exact(
					  bool const&al,
#ifdef falcON_INDI
					  real const&Ns,
					  real const&em,
					  real const&ex,
					  uint const&Nr,
					  real const&fe,
#endif
					  bool const&re_use_mem)
{
  if(al) {
    Nss = N_souls();
    Ncs = N_cells();
  } else
    update_and_pass_flags(this,Nss,Ncs);           // update &pass; count active
  const bool all = al || Nss==N_souls();
#ifdef falcON_INDI
  if(SOFT==individual)                             // IF individual softening   
    update_adjust_and_pass_eph(this,all,Ns,em,ex,Nr,fe);// get eph_i = eps_i/2  
#endif
  if(re_use_mem) set_soul_sinkpt(all);             // souls: give sinkpt        
  else         reset_soul_sinkpt(all);             // souls: give sinkpt        
  if(all) LoopMySouls(Si) Si->reset_srce();
  else    LoopMySouls(Si) if(is_active(Si)) Si->reset_srce();
}
//------------------------------------------------------------------------------
void grav_tree::exact_gravity(kern_type  const&KERNEL,
			      grav_stat *const&STATS,
			      real       const&EPS,
			      bool       const&al
#ifdef falcON_INDI
			     ,real       const&Nsoft,
			      uint       const&Nref,
			      real       const&emin,
			      real       const&efac
#endif
			      )
{
  prepare_grav_exact(al,
#ifdef falcON_INDI
		     Nsoft,emin,EPS,Nref,efac,
#endif
		     0);
  const bool all = al || Nss==N_souls();
  if(N_active_cells()==0)
    return warning("[grav_tree::exact_gravity()]: nobody active");
  STATS->reset();
  if(all) {
    grav_iact_all K(STATS,EPS,KERNEL
#ifdef falcON_INDI
		    ,SOFT
#endif
		    );
    K.direct_summation(root());
    LoopMySouls(Si) Si->normalize_grav();
  } else {
    grav_iact K(STATS,EPS,KERNEL
#ifdef falcON_INDI
		,SOFT
#endif
		);
    K.direct_summation(root());
    LoopMySouls(Si) if(is_active(Si)) Si->normalize_grav();
  }
  update_grav_eps(
#ifdef falcON_INDI
		  Nsoft,
#else
		  false,
#endif
		  all);
}
//------------------------------------------------------------------------------
inline
void grav_tree::prepare_grav_approx(const grav_mac* const&MAC,
				    bool            const&al,
				    bool            const&give_coeffs,
#ifdef falcON_INDI
				    real            const&Ns, 
				    real            const&em,
				    real            const&ex,
				    uint            const&Nr,
				    real            const&fe,
#endif
				    bool            const&re_use_mem)
{
  report REPORT("grav_tree::prepare_grav_approx()");
  if(al) {
    Nss = N_souls();
    Ncs = N_cells();
  } else
    update_and_pass_flags(this,Nss,Ncs);           // update &pass; count active
  const bool all = al || Nss==N_souls();
#ifdef falcON_INDI
  if(SOFT==individual)                             // IF individual softening   
    update_adjust_and_pass_eph(this,all,Ns,em,ex,Nr,fe);// get eph_i = eps_i/2  
#endif
  if(re_use_mem) set_soul_sinkpt(all);             // souls: give sinkpt        
  else         reset_soul_sinkpt(all);             // souls: give sinkpt        
  if(all) LoopMySouls(Si) Si->reset_srce();
  else    LoopMySouls(Si) if(is_active(Si)) Si->reset_srce();
  if(give_coeffs)                                  // IF(cell coeffs to give)   
    if(re_use_mem) set_cell_coeffs();              //   active cells: give coeff
    else         reset_cell_coeffs();              //   active cells: give coeff
  evaluate_poles <falcON_ORDER>(this);             // compute the multipoles    
  normalize_poles<falcON_ORDER>(this);             // normalize poles           
  MAC->set_rcrit(this);                            // set r_crit for all cells  
}
//------------------------------------------------------------------------------
void grav_tree::approx_gravity(const grav_mac*const&GMAC,
			       kern_type  const&KERNEL,
			       grav_stat *const&STATS,
			       real       const&EPS,
			       bool       const&al,
			       bool       const&split,
#ifdef falcON_INDI
			       real       const&Nsoft,
			       uint       const&Nref,
			       real       const&emin,
			       real       const&efac,
#endif
			       const int        direct[4])
{
  report REPORT("grav_tree::approx_gravity()");
  prepare_grav_approx(GMAC,al,                     // prepare for grav          
#ifdef _OPENMP
		      (!split)? 1 :                // only if openmp & !split   
#endif
		      0,
#ifdef falcON_INDI
		      Nsoft,emin,EPS,Nref,efac,
#endif
		      false);
  const bool all = al || Nss==N_souls();
  if(N_active_cells()==0)                          //   DONE if nobody active   
    return warning("[grav_tree::approx_gravity()]: nobody active");
  STATS->reset();                                  //   reset iaction statistics

#ifdef _OPENMP
  if(!split)                                       // IF openmp && !splitting   
    if(all) {                                      //   IF all assumed active   
      grav_iact_all GK(STATS,EPS,KERNEL,           //     init gravity kernel   
#  ifdef falcON_INDI
		       SOFT,
#  endif
		       direct,false);
      SelfInteractorP<grav_iact_all> MI(&GK,root(),depth());
      MI.interact_parallel();                      //     interaction phase     
      GK.evaluate(root());                         //     evaluation phase      
      Ncoeffs = N_active_cells();                  //     # coeffs used         
    } else {                                       //   ELSE                    
      grav_iact GK(STATS,EPS,KERNEL,               //     init gravity kernel   
#  ifdef falcON_INDI
		   SOFT,
#  endif
		   direct,false);
      SelfInteractorP<grav_iact> MI(&GK,root(),depth());
      MI.interact_parallel();                      //     interaction phase     
      GK.evaluate(root());                         //     evaluation phase      
      Ncoeffs = N_active_cells();                  //     # coeffs used         
    }
  else {                                           // ELSE (no openmp or split) 
#endif
    report REPORT2("interaction & evaluation");
    if(all) {                                      //   IF all are active       
      register uint NP = split? 4+N_cells()/8 : N_cells();
                                                   //     initial size: C_i pool
      grav_iact_all_s GK(STATS,EPS,NP,KERNEL,      //     init gravity kernel   
#ifdef falcON_INDI
			 SOFT,
#endif
			 direct,false);
      MutualInteractor<grav_iact_all_s> MI(&GK,split? depth()-1 : depth());
                                                   //   init mutual interactor  
      if(split) {                                  //     IF splitting          
	LoopCellKids(cell_iter,root(),c1) {        //       LOOP cell kids c1   
	  report::info("interaction");
	  MI.cell_self(c1);                        //         self-iaction c1   
	  LoopCellSecd(cell_iter,root(),c1+1,c2)   //         LOOP kids c2>c1   
	    MI.cell_cell(c1,c2);                   //           interaction c1,2
	  LoopSoulKids(cell_iter,root(),s2)        //         LOOP soul kids s  
	    MI.cell_soul(c1,s2);                   //           interaction c1,s
	  report::info("evaluation");
	  GK.evaluate(c1);                         //         evaluation phase  
	}                                          //       END LOOP            
	LoopSoulKids(cell_iter,root(),s1) {        //       LOOP soul kids s1   
	  report::info("interaction");
	  LoopSoulSecd(cell_iter,root(),s1+1,s2)   //         LOOP kids s2>s1   
	    GK.interact(s1,s2);                    //           interaction s1,2
	  report::info("evaluation");
	  s1->normalize_grav();                    //         evaluation phase  
	}                                          //       END LOOP            
      } else {                                     //     ELSE                  
	report::info("interaction");
	MI.cell_self(root());                      //       interaction phase   
	report::info("evaluation");
	GK.evaluate(root());                       //       evaluation phase    
      }                                            //     ENDIF                 
      Ncoeffs = GK.coeffs_used();                  //     remember # coeffs used
    } else {                                       //   ELSE: not all are active
      register uint NP = split? 4+N_active_cells()/8 : N_active_cells();
                                                   //     initial size: C_i pool
      grav_iact_s GK(STATS,EPS,NP,KERNEL,          //     init gravity kernel   
#ifdef falcON_INDI
		     SOFT,
#endif
		     direct,false);
      MutualInteractor<grav_iact_s> MI(&GK,split? depth()-1 : depth());
                                                   //     init mutual interactor
      if(split) {                                  //     IF splitting          
	LoopCellKids(cell_iter,root(),c1) {        //       LOOP cell kids c1   
	  if(is_active(c1)) {                      //        IF active s1:      
	    report::info("interaction");
	    MI.cell_self(c1);                      //         self-iaction c1   
	    LoopCellSecd(cell_iter,root(),c1+1,c2) //         LOOP kids c2>c1   
	      MI.cell_cell(c1,c2);                 //           interaction c1,2
	    LoopSoulKids(cell_iter,root(),s2)      //         LOOP soul kids s  
	      MI.cell_soul(c1,s2);                 //           interaction c1,s
	    report::info("evaluation");
	    GK.evaluate(c1);                       //         evaluation phase  
	  } else {                                 //        ELSE: inactive c1  
	    report::info("interaction");
	    LoopCellSecd(cell_iter,root(),c1+1,c2) //         LOOP kids c2>c1   
	      if(is_active(c2))MI.cell_cell(c1,c2);//           interaction c1,2
	    LoopSoulKids(cell_iter,root(),s2)      //         LOOP soul kids s  
	      if(is_active(s2))MI.cell_soul(c1,s2);//           interaction c1,s
	    report::info("no evaluation");
	  }                                        //        ENDIF              
	}                                          //       END LOOP            
	LoopSoulKids(cell_iter,root(),s1) {        //       LOOP soul kids s1   
	  if(is_active(s1)) {                      //        IF active s1:      
	    report::info("interaction");
	    LoopSoulSecd(cell_iter,root(),s1+1,s2) //         LOOP kids s2>s1   
	      GK.interact(s1,s2);                  //           interaction s1,2
	    report::info("evaluation");
	    s1->normalize_grav();                  //         evaluation phase  
	  } else {                                 //        ELSE: inactive s1  
	    report::info("interaction");
	    LoopSoulSecd(cell_iter,root(),s1+1,s2) //         LOOP kids s2>s1   
	      if(is_active(s2)) GK.interact(s1,s2);//           interaction s1,2
	    report::info("no evaluation");
	  }                                        //        ENDIF              
	}                                          //       END LOOP            
      } else {                                     //     ELSE                  
	report::info("interaction");
	MI.cell_self(root());                      //       interaction phase   
	report::info("evaluation");
	GK.evaluate(root());                       //       evaluation phase    
      }                                            //     ENDIF                 
      Ncoeffs = GK.coeffs_used();                  //     remember # coeffs used
    }                                              //   ENDIF                   
#ifdef _OPENMP
  }                                                // ENDIF                     
#endif
  update_grav_eps(                                 // update bodies             
#ifdef falcON_INDI
		  Nsoft,
#else
		  false,
#endif
		  all);
}
//------------------------------------------------------------------------------
inline void grav_tree::prepare_count_neighbours(
						real const &SIZE,
						bool const &re_use_mem)
{
  update_and_pass_flags(this,Nss,Ncs);             // update &pass; count active
  if(re_use_mem) set_soul_sinkpt();                // active souls: give sinkpt 
  else         reset_soul_sinkpt();                // active souls: give sinkpt 
#ifdef falcON_INDI
  if(SOFT == individual) {                         // IF(individual sizes) THEN 
    if       (use_sbodies()) {
      LoopMySouls(Si) if(is_active(Si)) {
	Si->copy_eph(my_sbodies());                //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
#ifdef falcON_MPI
    } else if(use_pbodies()) { 
      LoopMySouls(Si) if(is_active(Si)) {
	Si->copy_eph(my_pbodies());                //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
#endif
    } else if(use_barrays()) {
      LoopMySouls(Si) if(is_active(Si)) {
	Si->copy_eph(my_barrays());                //     get size == eps       
	Si->sizeq() = square(size(Si));            //     set size^2            
	Si->num()   = 0u;                          //     reset counter         
      }
    } else
      falcON_Error("tree has neither bodies nor arrays for data");
  } else                                           // ELSE(global size)         
#endif
    LoopMySouls(Si) if(is_active(Si)) {            //   LOOP active souls       
      Si->num()   = 0u;                            //     reset counter         
    }                                              // ENDIF                     
  register real smax;                              // cell size                 
  LoopMyCellsUp(Ci) {                              // LOOP cells upwards        
    smax = zero;                                   //   reset cell size         
#ifdef falcON_INDI
    if(SOFT == individual) {
      LoopSoulKids(cell_iter,Ci,s) if(is_active(s)){ // LOOP active soul kids   
	update_max(smax, sqrt(dist_sq(pos(Ci),pos(s))) + size(s));//update size 
      }                                            //   END LOOP                
    } else {
#endif
      LoopSoulKids(cell_iter,Ci,s) if(is_active(s)){ // LOOP active soul kids   
	update_max(smax, sqrt(dist_sq(pos(Ci),pos(s))) + SIZE); // update size  
      }                                            //   END LOOP                
#ifdef falcON_INDI
    }
#endif
    LoopCellKids(cell_iter,Ci,c) if(is_active(c)){ //   LOOP active cell kids   
      update_max(smax, sqrt(dist_sq(pos(Ci),pos(c))) + size(c)); // update size 
    }                                              //   END LOOP                
    Ci->size() = smax;                             //   set cell size           
  }                                                // END LOOP                  
}
//------------------------------------------------------------------------------
void grav_tree::count_neighbours(real const&EPS)
{
  if(use_sbodies() && !my_sbodies()->has(io::n) ||
#ifdef falcON_MPI
     use_pbodies() && !my_pbodies()->has(io::n) ||
#endif
     use_barrays() && !my_barrays()->has(io::n))
    falcON_ErrorF("nobody has memory for num","grat_tree::count_neighbours()");
#ifdef falcON_INDI
  switch(SOFT) {
  case individual: {
    prepare_count_neighbours();
    if(N_active_souls()==0)
      return warning("[grav_tree::count_neighbours()]: nobody active");
    neighbour_counter<grav_tree,individual> count;
    MutualInteractor<neighbour_counter<grav_tree,individual> > 
    MI(&count,depth());
    MI.cell_self(root());
  } break;
  case global: {
#endif
    if(EPS == zero)
      return warning("[grav_tree::count_neighbours()]: eps=0 -> no neighbours");
    prepare_count_neighbours(EPS);
    if(N_active_souls()==0)
      return warning("[grav_tree::count_neighbours()]: nobody active");
    neighbour_counter<grav_tree,global> count(EPS);
    MutualInteractor<neighbour_counter<grav_tree,global> > 
    MI(&count,depth());
    MI.cell_self(root());
#ifdef falcON_INDI
  } break;
  }
#endif
  if     (use_sbodies())
    LoopMySouls(Si) { if(is_active(Si)) Si->update_num(my_sbodies()); }
#ifdef falcON_MPI
  else if(use_pbodies())
    LoopMySouls(Si) { if(is_active(Si)) Si->update_num(my_pbodies()); }
#endif
  else if(use_barrays())
    LoopMySouls(Si) { if(is_active(Si)) Si->update_num(my_barrays()); }
  else
    falcON_Error("tree has neither bodies nor arrays for data");
}
////////////////////////////////////////////////////////////////////////////////
