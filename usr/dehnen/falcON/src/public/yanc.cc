//-----------------------------------------------------------------------------+
//                                                                             |
// yanc.cc                                                                     |
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
//                        YET ANOTHER N-BODY CODE                              |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <yanc.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>

#ifdef unix
#  define GIVE_HOST_INFO
extern "C" {
#  include <unistd.h>
#  include <pwd.h>
}
#endif
#include <public/nbdy.h>
#include <public/nmio.h>
#include <public/ionl.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
namespace nbdy{
  //===========================================================================+
  //                                                                           |
  // auxiliary constants and functions for the installation of class ndata     |
  //                                                                           |
  //============================================================================
#ifdef ALLOW_INDI
  inline basic_nbody::soft_type soften(const indx SOFTEN) {
    switch(SOFTEN) {
    case 1:  return basic_nbody::individual_fixed;
    case 2:  return basic_nbody::individual_adaptive;
    default: return basic_nbody::global_fixed; 
    }
  }
#endif
  //----------------------------------------------------------------------------
  inline kern_type kern(const indx KERN) {
    switch(KERN % 10) { 
    case  0: return p0;
    case  1: return p1;
    case  2: return p2;
    case  3: return p3;
    case  9: return newton;
    default: warning("kernel unknown, defaulting to P1");
      return Default::kernel;
    }
  }
  //----------------------------------------------------------------------------
  const char* opts[] = {"start_data",      //  0
			"theta",           //  1
			"hgrow",           //  2  replaces reuse                
			"Ncrit",           //  3
			"root_center",     //  4  not used, abolished           
			"reuse",           //  5  not used, replaced by hgrow   
			"softening",       //  6 used is ALLOW_INDI is defined  
			"Nsoft",           //  7 used is ALLOW_INDI is defined  
			"eps",             //  8
			"eps_given",       //  9
			"kernel",          // 10
			"Grav",            // 11
			"Nref",            // 12 used is ALLOW_INDI is defined  
			"facc",            // 13
			"fpot",            // 14
			"fcom",            // 15
			"time",            // 16
			"hmin",            // 17
			"Nsteps",          // 18
			"give_pot",        // 19
			"give_rho",        // 20
			"data_format",     // 21
			"N",               // 22
			"##","#",          // 23, 24
			"pot_file",        // 25 not used, abolished            
			"r_max",           // 26 not used, abolished            
			"data_file"};      // 27
  const unsigned nopt=28;
  //----------------------------------------------------------------------------
#ifdef ALLOW_INDI
  inline nbdy::uint check_eps_range(const real EPS,
				    nbdy::bodies *BODIES)
  {
    register nbdy::uint alright=0u;
    LoopBodies(bodies,BODIES,Bi)
      if(nbdy::eps(Bi) > EPS) {                       // eps out of range ?     
	Bi.eps() = EPS;                               // put it in range        
	alright++;
      }
    return alright;
  }
#endif
  //===========================================================================+
  //                                                                           |
  // class nbdy::ndata                                                         |
  //                                                                           |
  // auxiliary for the installation of class yanc                              |
  //                                                                           |
  //============================================================================
  class ndata {
#ifdef ALLOW_NEMO
    int                    ND;               // # NEMO devices                  
    nemo_out              *NEMO;             // NEMO output devices             
  public:
    real                   CPU_I;            // initial CPU time                
#endif
  public:
    std::string            FILE;             // input file                      
    pot_provider          *PP;               // pot provider                    
    const extpot          *PEX;              // external potential              
    real                   THETA;            // opening angle                   
    int                    HGROW;            // grow fresh tree every 2^h esteps
    int                    NCRIT;            // max # bodies in unsplit cell    
#ifdef ALLOW_INDI
    real                   NSOFT;            // # bodies in softening spheres   
    uint                   NREF;             // # bodies for density estimate   
    indx                   SOFTEN;           // softening type                  
    bool                   EPS_GIVEN;        // eps_i are given on input        
#endif
    real                   EPS;              // eps / eps_max                   
    indx                   KERN;             // (0/1/2/3/9) soft kernel         
    double                 G,iG;             // constant of gravity, reciprocal 
    real                   FACC;             // factor f_acc for time-stepping  
    real                   FPOT;             // factor f_pot for time-stepping  
    real                   FCOM;             // factor f_com for time-stepping  
    real                   TINI;             // initial simulation time         
    real                   TIME;             // simulation time                 
    io                     WRITE;            // what to write to output         
    int                    HMIN;             // tau_min = 1/2^h0                
    indx                   NSTEPS;           // number of tau                   
    indx                   FORMAT;           // data format: 0/1 : ascii/binary 
    bodies                *BODIES;           // bodies themselves               
  private:
    ndata();                                 // not implemented                 
    ndata(const ndata&);                     // not implemented                 
    //--------------------------------------------------------------------------
    io bodydatabits() const {
      register io need = io::mxvpaf;         // basic data bits                 
#ifdef ALLOW_INDI
      if(SOFTEN) need |= io::e;              // need eps_i                      
#endif
      if(NSTEPS>1) need |= io::l;            // need levels                     
      return need;
    }
    //--------------------------------------------------------------------------
  private:
    void read_yanc        ()                     // read input using yanc format
    {
      std::ifstream in;
      open_error(in,FILE.c_str());
      char option[100];
      register uint iopt = 0;
      register uint Ns,Nbody;
      do {
	if(in)
	  skip_char(in,'#') >> option;
	else {
	  std::cerr<<"### YANC: end of input seen before \"start_data\"\n";
	  return;
	}
	for(iopt=0; iopt<nopt && std::strcmp(option,opts[iopt]); iopt++);
	switch(iopt) {
	case  0:                   SwallowRestofLine(in); break;
	case  1: in >> THETA;      SwallowRestofLine(in); break;
	case  2: in >> HGROW;      SwallowRestofLine(in); break;
	case  3: in >> NCRIT;      SwallowRestofLine(in); break;
	case  4: std::cerr<<"### YANC parameter \"root_center\" deprecated\n";
	                           SwallowRestofLine(in); break;
	case  5: std::cerr<<"### YANC parameter \"reuse\" deprecated;"
			  <<" use \"hgrow\" instead\n";
	                           SwallowRestofLine(in); break;
#ifdef ALLOW_INDI
	case  6: in >> SOFTEN;     SwallowRestofLine(in); break;
	case  7: in >> NSOFT;      SwallowRestofLine(in); break;
#else
	case  6: std::cerr<<"### YANC parameter \"softening\" disabled\n";
	                           SwallowRestofLine(in); break;
	case  7: std::cerr<<"### YANC parameter \"Nsoft\" disabled\n";
	                           SwallowRestofLine(in); break;
#endif
	case  8: in >> EPS;        SwallowRestofLine(in); break;
	case  9: 
#ifdef ALLOW_INDI
	         EPS_GIVEN = true;
#else
		 std::cerr<<"### YANC parameter \"eps_given\" disabled\n";
#endif
	                           SwallowRestofLine(in); break;
	case 10: in >> KERN;       SwallowRestofLine(in); break;
	case 11: in >> G; iG=1./G; SwallowRestofLine(in); break;
#ifdef ALLOW_INDI
	case 12: in >> NREF;       SwallowRestofLine(in); break;
#else
	case 12: std::cerr<<"### YANC parameter \"Nref\" disabled\n";
	                           SwallowRestofLine(in); break;
#endif
	case 13: in >> FACC;       SwallowRestofLine(in); break;
	case 14: in >> FPOT;       SwallowRestofLine(in); break;
	case 15: in >> FCOM;       SwallowRestofLine(in); break;
	case 16: in >> TIME;       SwallowRestofLine(in); break;
	case 17: in >> HMIN;       SwallowRestofLine(in); break;
	case 18: in >> Ns;         SwallowRestofLine(in);
	         if(Ns>=1) NSTEPS = Ns;
	         else std::cerr<<"### YANC: Nsteps>0 required\n"; break;
	case 19: WRITE |= io::p;   SwallowRestofLine(in); break;
	case 20: WRITE |= io::r;   SwallowRestofLine(in); break;
	case 21: in >> FORMAT;     SwallowRestofLine(in); break;
	case 22: in >> Nbody;      SwallowRestofLine(in); break;
	case 23: case 24:          SwallowRestofLine(in); break;
	case 25: std::cerr<<"### YANC parameter \"pot_file\" deprecated\n";
			           SwallowRestofLine(in); break;
	case 26:                   SwallowRestofLine(in); break;
	case 27: in >> FILE;       read_yanc();                 return;
	default: std::cerr<<"### YANC: unknown parameter \""<<option
			  <<"\" ignored\n";
	                           SwallowRestofLine(in); break;
	}
      } while(iopt);
      if(BODIES) BODIES->reset(Nbody,bodydatabits());
      else       MemoryCheck(BODIES = new bodies(Nbody,bodydatabits()));
#ifdef ALLOW_INDI
#define WRITEOUT(WR) EPS_GIVEN? WR | io::e : WR
#else
#define WRITEOUT(WR) WR
#endif
      if(FORMAT == 0)
	BODIES->read_yanc_ascii (in,WRITEOUT(io::mxv));
      else
	BODIES->read_yanc_binary(in,WRITEOUT(io::mxv));
      if(! in) {
	std::cerr<<"### YANC: too few bodies in input\n";
	delete BODIES;
	BODIES = 0;
	return;
      }
#ifdef ALLOW_INDI
      if(SOFTEN && EPS_GIVEN) {
	uint a = check_eps_range(EPS,BODIES);
	if(NSOFT>zero && a>0) std::cerr<<"### YANC: "<<a<<" eps_i > eps_max\n";
      }
#endif
      if(G != 1.0 && G != 0.)                      // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= G;                          //     set M->G*M           <
    }
    //--------------------------------------------------------------------------
  private:
    void write_yanc_header(                       // write header in yanc format
			   std::ostream&out,      // I/O:  stream to write to   
			   const char *file,      // I:    name of file written 
#ifdef ALLOW_INDI
			   const bool  eps_i,     // I:    add param eps_given  
#endif
			   const char *prog,      // I:    name of calling main 
			   const indx  format)    // I:    format type          
      const 
    {
#ifdef GIVE_HOST_INFO
      char host[100]; gethostname(host,100);
      time_t now = ::time(0);
#endif
      out.setf(std::ios::left);
      out<<std::setfill(' ');
      out<<"## -------------------------------"
	"-----------------------------------\n"
	 <<"##  file          \""<<file<<"\"\n"
#ifdef GIVE_HOST_INFO
	 <<"##  generated on  "  <<ctime(&now);
      if(prog) out<<"##  with program  \""<<prog<<"\"\n";
      out<<"##  by user       \""<<(getpwuid(geteuid())->pw_name)<<"\"\n"
	 <<"##  on machine    \""<<host<<"\"\n";
#else
      ;
      if(prog) out<<"##  generated with program  \""<<prog<<"\"\n";
#endif
      if(!FILE.empty())
	out<<"##\n"
	   <<"##  original data read from file \""<<FILE<<"\"\n";
      out<<"##  initial time = "<<TINI<<"\n"
	 <<"## -------------------------------"
	"-----------------------------------\n";
      // tree issues
      out<<"# theta          "<<std::setw(15)<<THETA;
      if(THETA<zero) out<<" constant opening angle\n";
      else           out<<" mass-dependent opening angle; theta_min\n";
      out<<"# hgrow          "<<std::setw(15)<<HGROW
	 <<" grow fresh tree every 2^"<<HGROW<<"smallest steps\n"
	 <<"# Ncrit          "<<std::setw(15)<<NCRIT
	 <<" max # of bodies in un-splitted cell\n";
      // softening issues
#ifdef ALLOW_INDI
      out<<"# softening      "<<std::setw(15)<<SOFTEN;
      switch(SOFTEN) {
      case 1:
	out<<" individually fixed softening lengths\n";
	break;
      case 2:
	out<<" individually adaptive softening lengths with\n"
	   <<"# Nsoft          "<<std::setw(15)<<NSOFT
	   <<" # bodies in softening sphere\n"
	   <<"# eps            "<<std::setw(15)<<EPS
	   <<" maximum eps_i\n";
	break;
      default:
	out<<" globally fixed softening length:\n"
	   <<"# eps            "<<std::setw(16)<<EPS<<"\n";
      }
#else
      out<<"# eps            "<<std::setw(16)<<EPS<<"\n";
#endif
      // time-stepping issues
      if(NSTEPS==1) 
	out<<"# Nsteps         1               (ordinary leap-frog)\n"
	   <<"# hmin           "<<std::setw(15)<<HMIN<<" with step tau "
	   <<tau_min()<<"\n";
      else {
	out<<"# facc           "<<std::setw(15)<<FACC;
	if(FACC) out<<"tau <= "<<FACC<<" / |acc|\n"; else out<<"\n";
	out<<"# fpot           "<<std::setw(15)<<FPOT;
	if(FPOT) out<<"tau <= "<<FPOT<<" / |pot|\n"; else out<<"\n";
	out<<"# fcom           "<<std::setw(15)<<FCOM;
	if(FCOM) out<<"tau <= "<<FCOM<<" * sqrt|pot| / |acc|\n"; else out<<"\n";
	out<<"# hmin           "<<std::setw(15)<<HMIN<<" tau_min = 2^"<<(-HMIN)
	   <<" = "<<tau_min()<<"\n";
	out<<"# Nsteps         "<<std::setw(15)<<NSTEPS<<" tau_max = 2^"
	   <<(NSTEPS-HMIN-1)
	   <<" = "<<tau_max()<<"\n";
      }
      // data issues
      out<<"# time           "<<std::setw(15)<<TIME<<" simulation time\n"
	 <<"# N              "<<std::setw(15)<<BODIES->N_bodies()
	 <<" number of bodies\n"
	 <<"# Grav           "<<std::setw(15)<<G<<"\n"
	 <<"# data_format    "<<std::setw(15)<<format;
      if(format==0) out<<" bodies in yanc ascii format\n";
      else          out<<" bodies in yanc binary format\n";
      if(format==0) {
	if(WRITE & io::p)
	  out<<"# give_pot                         potential given on output\n";
	if(WRITE & io::r)
	  out<<"# give_rho                         density given on output\n";
      }
#ifdef ALLOW_INDI
      if(eps_i)
	out<<"# eps_given                        eps_i given with data below\n";
#endif
      out<<"## -----------------------------------"
	"--------------------------------\n"
	 <<"# start_data"<<std::endl;
    }
    //--------------------------------------------------------------------------
#ifdef ALLOW_NEMO
    //--------------------------------------------------------------------------
  public: void open_nemo(const int  d=0,
			 const char*file=0, 
			 const bool resume=false)
    {
      if(d >= ND) error("[open_nemo()]: nemo device does not exist");
      if(resume) {
	if(file==0 || FILE == std::string(file))
	  (NEMO+d)->open_to_append(FILE.c_str());
	else {
	  (NEMO+d)->open(file);
	  (NEMO+d)->write_history();
	}
      } else if(file) {
	if(FILE == std::string(file) && FILE != "-")
	  error("out==in; use option resume instead");
	(NEMO+d)->open(file);
	(NEMO+d)->write_history();
      } else
	warning("cannot open nemo output");
    }
    //--------------------------------------------------------------------------
  public: void close_nemo  (const int d=0)
    {
      if(d < ND) (NEMO+d)->close();
    }
    //--------------------------------------------------------------------------
  public: int nemo_devices() const { return ND; }
    //--------------------------------------------------------------------------
  public: bool nemo_is_open(const int d=0) const
    {
      if(d >= ND) return false;
      return (NEMO+d)->is_open();
    }
    //--------------------------------------------------------------------------
  public:
    void read_nemo(                                // read from NEMO snapshot   
		   const io   want  =io::mxv,      //[I: what to read]          
		   const bool resume=false)        //[I: resume old simul?]     
    {
      const io support = want | bodydatabits();
      if(!BODIES) 
	{ MemoryCheck(BODIES=new bodies(0,support)); }// initialize bodies      
      else BODIES->reset(BODIES->N_bodies(),support); // or reset bodies        
      nemo_in NEMO(FILE.c_str());                  // open nemo input           
      NEMO.read_history();                         // read history              
      register io got;                             // what did we get?          
      if(resume) {                                 // IF (resuming old sim)    >
	do {                                       //   do: read particles     >
	  NEMO.open_set(nemo_io::snap);            //     open nemo snapshot    
	  BODIES->read_nemo_particles(NEMO,got,&TIME,want);// read particles    
	  if(NEMO.is_present(nemo_io::diags)) {    //     IF(diags set)      >  
	    NEMO.open_set(nemo_io::diags);         //       open diags set      
	    if(NEMO.is_present(nemo_io::cputime))  //       IF(cpu time there)  
	      CPU_I=NEMO.read(nemo_io::cputime)*60;//         read CPU [min]    
	    NEMO.close_set(nemo_io::diags);        //       close diags set     
	  }                                        //     <                     
	  NEMO.close_set(nemo_io::snap);           //     close nemo snapshot   
	} while(NEMO.is_present(nemo_io::snap));   //   < while more to be read 
      } else {                                     // < ELSE (not resuming) >   
	  NEMO.open_set(nemo_io::snap);            //   open nemo snapshot      
	  BODIES->read_nemo_particles(NEMO,got,&TIME,want);// read particles    
	  NEMO.close_set(nemo_io::snap);           //   close nemo snapshot     
      }
#ifdef ALLOW_INDI
      if(got & io::e) EPS_GIVEN = true;            // set EPS_GIVEN             
#endif
      if(G != 1.0 && G != 0.)                      // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= G;                          //     set M->G*M           <
    }
    //--------------------------------------------------------------------------
  public:
    void write_nemo (const io          w=io::mxv,  //[I: what to write out]     
		     const basic_nbody*NBODY=0,    //[I: nbody for diags O]     
		     const int         d    =0)    //[I: index of nemo stream]  
    {
      if(d >= ND)
	error("[write_nemo()]: nemo device does not exist");
      nemo_out *OUT=NEMO+d;
      if(!OUT->is_open())
	error("nemo output stream not opened; cannot write snapshot");
      if(!BODIES) {
	warning("no body exist; cannot write snapshot");
	return;
      }
      if(w & io::m && G != 1.0 && G != 0.)         // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= iG;                         //     set G*M->M           <
      register uint        i,j;
      OUT->open_set(nemo_io::snap);                // open a new nemo snapshot  
      BODIES->write_nemo_particles(*OUT, &TIME,    //   write out particles     
#ifdef ALLOW_INDI
				   SOFTEN? w | io::e : w
#else
				   w
#endif
				   );
      if(w & io::m && G != 1.0 && G != 0.)         // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= G;                          //     set G*M->M           <
      if(NBODY) {                                  //   IF diags output wantd  >
	OUT->open_set(nemo_io::diags);             //     open diagnostics set  
	OUT->single_vec(0)=NBODY->total_energy();  //       copy total energy   
	OUT->single_vec(1)=NBODY->kin_energy();    //       copy kin energy     
	OUT->single_vec(2)=NBODY->pot_energy();    //       copy pot energy     
	OUT->write(nemo_io::energy);               //       write energies      
	for(i=0;i!=NDIM;++i) for(j=0;j!=NDIM;++j)  //     loop dims twice    > >
	  OUT->single_mat(i,j) = NBODY->kin_energy(i,j);//   copy K_ij       < <
	OUT->write(nemo_io::KinT);                 //     write K_ij            
	for(i=0;i!=NDIM;++i) for(j=0;j!=NDIM;++j)  //     loop dims twice    > >
	  OUT->single_mat(i,j) = NBODY->pot_energy(i,j);//   copy W_ij       < <
	OUT->write(nemo_io::PotT);                 //     write W_ij            
	for(i=0;i!=NDIM;++i) for(j=0;j!=NDIM;++j)  //     loop dimensions    > >
	  OUT->single_mat(i,j) = NBODY->as_angmom(i,j); //   copy A_ij       < <
	OUT->write(nemo_io::AmT);                  //     write A_ij            
	for(i=0;i!=NDIM;++i) {                     //     loop dims            <
	  OUT->single_phs(0,i)=NBODY->center_of_mass()[i];// copy c-of-m pos    
	  OUT->single_phs(1,i)=NBODY->total_momentum()[i];// copy c-of-m vel    
	}                                          //     <                     
	OUT->write(nemo_io::cofm);                 //     write c-of-m (x,v)    
	OUT->write(nemo_io::cputime,               //     write accum CPU [min] 
		      NBODY->cpu_total()/60.);     //                           
	OUT->close_set(nemo_io::diags);            //     close diagnostics set 
      }                                            //   <                       
      OUT->close_set(nemo_io::snap);               // close nemo snapshot       
      OUT->reset();                                // reset nemo output         
    }
    //--------------------------------------------------------------------------
  public:
    ndata(const char*  const&file,
	  bool         const&nemo,
	  double       const&th,
	  int          const&hg,
	  int          const&nc,
	  double       const&ep,
	  int          const&kn,
#ifdef ALLOW_INDI
	  int          const&sf,
	  double       const&ns,
	  int          const&nr,
#endif
	  int          const&hm,
	  int          const&nl,
	  double       const&fa,
	  double       const&fp,
	  double       const&fc,
	  const extpot*const&pt,
	  const bool        resume   =false,        // only applies if nemo=true
	  const io          read_more=io::o,        // only applies if nemo=true
	  const int         nd       =2,            // # nemo devices           
	  double const&     gr       =1.0)
      :
      ND        (nemo? nd : 0),
      NEMO      (new nemo_out[ND]),
      CPU_I     (zero),
      FILE      (file),
      PEX       (pt),
      THETA     (real(th)),
      HGROW     (hg<0? 0 : hg),
      NCRIT     (uint(max(0,nc))),
#ifdef ALLOW_INDI
      NSOFT     (real(abs(ns))),
      NREF      (uint(max(0,nr))),
      SOFTEN    (sf),
      EPS_GIVEN (false),
#endif
      EPS       (real(ep)),
      KERN      (indx(max(0,kn))),
      G         (gr),
      iG        (1.0/G),
      FACC      (real(fa)),
      FPOT      (real(fp)),
      FCOM      (real(fc)),
      TIME      (zero),
      WRITE     (io::mxv),
      HMIN      (hm),
      NSTEPS    (indx(max(0,nl))),
      FORMAT    (0),
      BODIES    (0)
    {
      if(ND) { MemoryCheck(NEMO); }
      if(nemo) {
	register io  need  = io::mxv | read_more;
	read_nemo(need, resume);
      } else
	read_yanc();
      TINI = TIME;
#ifdef ALLOW_INDI
      if(SOFTEN == basic_nbody::individual_adaptive && HGROW)
	error("inidividual adaptive softening lengths must not be used"
	      "with hgrow != 0");
#endif
    }
#endif
    //--------------------------------------------------------------------------
  public:
    ndata(const char  *file,
	  pot_provider*pp=0) :
      FILE      (file),
      PP        (pp),
      THETA     (Default::theta),
      HGROW     (0),
      NCRIT     (Default::Ncrit),
#ifdef ALLOW_INDI
      NSOFT     (zero),
      NREF      (32u),
      SOFTEN    (0),
      EPS_GIVEN (false),
#endif
      EPS       (zero),
      KERN      (1),
      G         (1.0),
      iG        (1.0),
      FACC      (0.01),
      FPOT      (zero),
      FCOM      (zero),
      TIME      (zero),
      WRITE     (io::mxv),
      HMIN      (6),
      NSTEPS    (indx(1)),
      FORMAT    (0),
      BODIES    (0)
    {
      read_yanc();
      TINI = TIME;
#ifdef ALLOW_INDI
      if(SOFTEN == basic_nbody::individual_adaptive && HGROW)
	error("inidividual adaptive softening lengths must not be used"
	      "with hgrow != 0");
#endif
    }
    //--------------------------------------------------------------------------
  public:  ~ndata()
    {
#ifdef ALLOW_NEMO
      for(register int d=0; d!=ND; ++d) close_nemo(d);
      delete[] NEMO;
#endif
      if(BODIES) delete BODIES;
    }
    //--------------------------------------------------------------------------
  public:
    void describe(                                 // describe parameter        
		  std::ostream&out,                // I: output stream          
		  const char  *name,               // I: name of main           
		  const char  *output)             // I: base name: output files
    {
      if(!okay()) {
	out<<"#  ERROR occured during input of parameters and bodies\n";
	return;
      }
      out<<
	"# =======================================================================";
      if(name) out<<"\n# \""<<name<<"\"";
      out<<"\n#  parameters:"
	// tree issues
	 <<"\n#   theta                "<<THETA
	 <<"\n#   hgrow                "<<HGROW
	 <<"\n#   Ncrit                "<<NCRIT;
      // softening issues
#ifdef ALLOW_INDI
      switch(SOFTEN) {
      case 1:
	out<<"\n#   softening            individually fixed";
	if(EPS_GIVEN) out<<"\n#   eps_i                given with data";
	else          out<<"\n#   eps_i                to be initialized using"
			 <<"\n#   Nsoft                "<<NSOFT
			 <<"\n#   eps_max              "<<EPS;
	break;
      case 2:
	out<<"\n#   softening            individually adaptive"
	   <<"\n#   Nsoft                "<<NSOFT
	   <<"\n#   eps_max              "<<EPS;
	break;
      default:
	out<<"\n#   global eps           "<<EPS;
	break;
      }
#else
      out<<"\n#   global eps           "<<EPS;
#endif
      out<<"\n#   kernel               "<<nbdy::describe(kern_type(KERN));
      // time-stepping issues
      if(NSTEPS==1)
	out<<"\n#   time step            "<<tau_min();
      else {
	out<<"\n# facc           "<<std::setw(15)<<FACC;
	if(FACC) out<<"tau <= "<<FACC<<" / |acc|";
	out<<"\n# fpot           "<<std::setw(15)<<FPOT;
	if(FPOT) out<<"tau <= "<<FPOT<<" / |pot|";
	out<<"\n# fcom           "<<std::setw(15)<<FCOM;
	if(FCOM) out<<"tau <= "<<FCOM<<" * sqrt(-pot) / |acc|";
	out<<"\n#   tau_min              "<<tau_min()
	   <<"\n#   tau_max              "<<tau_max();
      }
      // data issues
      out<<"\n#   time                 "<<TIME
	 <<"\n#   N                    "<<BODIES->N_bodies()
	 <<"\n#   Grav                 "<<G
	 <<"\n#   input file           "<<FILE;
      if(output)
	out<<"\n#   output files         "<<output<<".XXXX";
      if(FORMAT)  out<<"\n#   data format          yanc binary";
      else        out<<"\n#   data format          yanc ascii";
      if(WRITE & io::p) out<<"\n#   potential            written on output";
      if(WRITE & io::r) out<<"\n#   density              written on output";
      out<<std::endl;
    }
    //--------------------------------------------------------------------------
  public:
    static void describe_yanc_format(std::ostream& o) {
      o<<
	" Format of YANC I/O files\n"
	" ========================\n"
	" The file starts with a header with one line per parameter.  The syntax\n"
	" of these lines is \"[#] PARAMETER VALUE(S) COMMENT\", where PARAMETER is\n"
	" a character string, see below, while  VALUE is the  numerical value of\n"
	" of the associated parameter. COMMENT is anything between VALUE and EOL\n"
	" and is ignored, as is an optional \"#\" at  the beginning of the  header\n"
	" lines. The meaning of the various parameters is as follows.\n\n"

	" PARAMETER   VALUE default meaning\n"
	" ---------------------------------------------------------------------\n"
	" ##                      this line is (quietly) ignored\n"
	" data_file    FILE       ignore rest of file and use FILE instead\n"
	" theta          th  0.5  opening angle  < 0 -> th=|VALUE| = const\n"
	"                                        > 0 -> th=th(M), th_min=TH\n"
	" hgrow          hg    0  grow fresh tree every 2^hg smallest steps\n"
	" Ncrit          Nc    6  max # bodies in un-splitted cell\n"
#ifdef ALLOW_INDI
	" softening       S    0  determines usage of softening lengths\n"
	"                           0 -> globally fixed\n"
	"                           1 -> individually fixed (NOT TESTED)\n"
	"                           2 -> individually adapted (NOT TESTED)\n"
	" Nsoft          Ns    0  # bodies in softening spheres (for S>0)\n"
#endif
	" eps            ep    0  fixed / maximum softening length\n"
#ifdef ALLOW_INDI
	" eps_given               individual eps_i are given on input\n"
#endif
	" kernel         kr    1  P_kr (P_0=Plummer) softening kernel\n"
	" Grav           G     1  Newton's constant of gravity\n"
#ifdef ALLOW_INDI
	" Nref           Nf   16  # bodies in density estimate with scheme=1\n"
#endif
	" facc           fa 0.01  factor in time-stepping  tau < fa/acc\n"
	" fpot           fp    0  factor in time-stepping  tau < fp/pot\n"
	" fcom           fc    0  factor in time-stepping  tau < fc*sqrt(pot)/acc\n"
	" hmin           hm    6  tau_min = (1/2)^hm\n"
	" Nsteps         Nt    1  # time steps (1 -> leap-frog)\n"
	" give_pot        g    0  0/1/2/3 give no/N-body pot/external pot/both\n"
	" give_rho                give mass density on output\n"
	" time            t    0  initial simulation time\n"
	" N               n    0  number of bodies to be read\n"
	" data_format   for    0  0/1 -> body data are in ascii/binary\n"
	" start_data              this terminates the header\n\n"
	" A parameter already  set may be superseeded  by another line.  A line\n"
	" with unknown PARAMETER is notified but otherwise ignored.\n\n"
	" The body data consist of an integer, whose bits indicate which data\n"
	" are actually given (see below) plus N lines (ASCII/binary) containing\n"
	" these data in this order (if given):\n"
	"   mass         (1    float) :  bit 0,  value   1\n"
	"   position     (NDIM floats):  bit 1,  value   2\n"
	"   velocity     (NDIM floats):  bit 2,  value   4\n"
#ifdef ALLOW_INDI
	"   eps_i        (1    float) :  bit 3,  value   8\n"
#endif
	"   potential    (1    float) :  bit 4,  value  16\n"
	"   acceleration (NDIM floats):  bit 5,  value  32\n"
	"   density      (1    float) :  bit 6,  value  64\n"
	"   auxiliary    (1    float) :  bit 7,  value 128\n"
	"   external pot (1    float) :  bit 8,  value 256\n";
    }
    //-------------------------------------------------------------------------+
  public:
    bool write_yanc_ascii (                        // write yanc format, ascii  
			   const char* file,       // I:  name: output file     
			   const char* prog) const // I:  name: calling program 
    {
      std::ofstream out;
      if(!open(out,file)) {
	warning("[write_yanc_ascii()]: cannot output");
	return false;
      }
      if(!okay()) {
	warning("[write_yanc_ascii()]: not okay(); cannot output");
	return false;
      }
      write_yanc_header(out,file,
#ifdef ALLOW_INDI
			SOFTEN,
#endif
			prog,0);
      if(WRITE && io::m && G != 1.0 && G != 0.)    // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= G;                          //     set M->G*M           <
      BODIES->write_yanc_ascii(out,WRITEOUT(WRITE));
      if(WRITE && io::m && G != 1.0 && G != 0.)    // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= iG;                         //     set G*M->M           <
      if(out) return true;
      else    return false;
    }
    //-------------------------------------------------------------------------+
  public:
    bool write_yanc_binary(                        // write yanc format, binary 
			   const char* file,       // I:  name: output file     
			   const char* prog) const // I:  name: calling program 
    {
      std::ofstream out;
      if(!open(out,file)) {
	std::cerr<<"YANC: write_yanc_binary(): cannot output\n";
	return false;
      }
      if(!okay()) {
	out<<" YANC: not okay(), cannot write_yanc_binary()\n";
	return false;
      }
      write_yanc_header(out,file,
#ifdef ALLOW_INDI
			SOFTEN,
#endif
			prog,1);
      if(WRITE && io::m && G != 1.0 && G != 0.)    // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= G;                          //     set M->G*M           <
      BODIES->write_yanc_binary(out,WRITEOUT(WRITE));
      if(WRITE && io::m && G != 1.0 && G != 0.)    // IF(G != 1)                
	LoopBodies(sbodies,BODIES,Bi)              //   loop bodies            >
	  Bi.mass() *= iG;                         //     set G*M->M           <
      if(out) return true;
      else    return false;
    }
    //-------------------------------------------------------------------------+
    // informators                                                              
    bool              okay ()       const { return BODIES!=0; }
    indx        const&Nsteps()      const { return NSTEPS; }
    uint        const&Nbodies()     const { return BODIES->N_bodies(); }
    int         const&hgrow()       const { return HGROW; }
    int         const&Ncrit()       const { return NCRIT; }
    std::string const&data_file()   const { return FILE; }
    real        const&eps()         const { return EPS; }
    real              tau_min()     const { return pow(half,HMIN); }
    real              tau_max()     const { return pow(half,HMIN+1-NSTEPS); }
    io          const&writing()     const { return WRITE; }
#ifdef ALLOW_INDI
    indx              softening()   const { return SOFTEN; }
    real        const&Nsoft()       const { return NSOFT; }
    uint        const&Nref()        const { return NREF; }
#endif
    indx              kernel()      const { return KERN; }
    double      const&Grav  ()      const { return G; }
    real        const&time  ()      const { return TIME; }
    real        const&initial_time()const { return TINI; }
    real        const&theta()       const { return THETA; }
    real        const&facc()        const { return FACC; }
    real        const&fpot()        const { return FPOT; }
    real        const&fcom()        const { return FCOM; }
  };
}
//=============================================================================#
//                                                                             |
// class nbdy::yanc                                                            |
//                                                                             |
//=============================================================================#
void yanc::describe(std::ostream&out,
		    const char  *name,
		    const char  *outf)
{
  MYNDATA->describe(out,name,outf);
  out<<
    "#   total energy:        "<<MYNBODY->total_energy()<<
    "\n# ==========================================================="
    "============\n#\n";
}
//------------------------------------------------------------------------------
bool yanc::write_yanc_ascii(const char* file, const char* prog)
{
#ifdef ALLOW_INDI
  if(MYNDATA->WRITE & io::r) MYNBODY->estimate_mass_densities();
#endif
  return MYNDATA->write_yanc_ascii(file,prog);
}
//------------------------------------------------------------------------------
bool yanc::write_yanc_binary(const char* file, const char* prog) const
{
  return MYNDATA->write_yanc_binary(file,prog);
}
//------------------------------------------------------------------------------
void yanc::describe_yanc_format(std::ostream& o)
{
  return ndata::describe_yanc_format(o);
}
//------------------------------------------------------------------------------
void yanc::full_step() 
{
  if(! MYNDATA->okay() )
    error("[yanc::full_step()]: something wrong, perhaps invalid data?");
  MYNBODY->full_step();
  MYNDATA->TIME = MYNBODY->time();
}
//------------------------------------------------------------------------------
void yanc::stats (std::ostream&o) const
{
  o<<" "; MYNBODY->stats(o);
}
//------------------------------------------------------------------------------
void yanc::stats_head  (std::ostream&o) const
{
  o<<"#"; MYNBODY->stats_head(o);
  o<<"#"; MYNBODY->stats_line(o);
}
//------------------------------------------------------------------------------
#define									       \
FUNC(TYPE,NAME) TYPE yanc::NAME() const { return TYPE(MYNDATA->NAME()); }
FUNC(bool,   okay)
FUNC(int,    Nsteps)
FUNC(int,    Nbodies)
FUNC(int,    hgrow)
FUNC(int,    Ncrit)
FUNC(io,     writing)
#ifdef ALLOW_INDI
FUNC(int,    softening)
FUNC(int,    Nref)
FUNC(float,  Nsoft)
#endif
FUNC(int,    kernel)
FUNC(double, Grav)
FUNC(float,  time)
FUNC(float,  initial_time)
FUNC(float,  theta)
FUNC(float,  facc)
FUNC(float,  fpot)
FUNC(float,  fcom)
FUNC(float,  eps)
FUNC(float,  tau_min)
FUNC(float,  tau_max)
#undef FUNC
float       yanc::cpu_total() const { return MYNBODY->cpu_total(); }
const char* yanc::data_file() const { return MYNDATA->data_file().c_str(); }
//------------------------------------------------------------------------------
#ifdef ALLOW_NEMO
//------------------------------------------------------------------------------
int yanc::nemo_devices() const {
  return MYNDATA->nemo_devices();
}
//------------------------------------------------------------------------------
bool yanc::nemo_is_open(const int d) const {
  return MYNDATA->nemo_is_open(d);
}
//------------------------------------------------------------------------------
void yanc::open_nemo (const int d, const char *file, const bool resume) {
  MYNDATA->open_nemo(d,file,resume);
}
//------------------------------------------------------------------------------
void yanc::close_nemo(const int d) { MYNDATA->close_nemo(d); }
//------------------------------------------------------------------------------
void yanc::write_nemo(const io w, const int d, const bool di)
{
#ifdef ALLOW_INDI
  if(w & io::r) MYNBODY->estimate_mass_densities();
#endif
  MYNDATA->write_nemo(w,di? MYNBODY : 0,d);
}
//------------------------------------------------------------------------------
void yanc::describe_nemo(std::ostream&out,
			 const char  *com)
{
  if(!okay())
    error("something wrong with N-body data, perhaps not read yet?");
  out<<"#"; MYNBODY->stats_line(out);
  out<<"# \""<<com<<"\"\n#\n";
#ifdef GIVE_HOST_INFO
  char host[100]; gethostname(host,100);
  time_t now = ::time(0);
  out<<"# run at  "<<ctime(&now)
     <<"#     by  \""<<(getpwuid(geteuid())->pw_name)<<"\"\n"
     <<"#     on  \""<<host<<"\"\n#\n";
#endif
  out.flush();
}
//------------------------------------------------------------------------------
yanc::yanc(const char     *in,                      // I: data input file       
	   const bool      nemo,                    // I: do nemo I/O?          
	   const double    th,                      // I: theta                 
	   const int       hg,                      // I: hgrow                 
	   const int       nc,                      // I: Ncrit                 
	   const double    ep,                      // I: eps                   
	   const int       kn,                      // I: kern                  
	   const int       hm,                      // I: hmin                  
	   const int       nl,                      // I: Nsteps                
	   const double    fa,                      // I: f_a                   
	   const double    fp,                      //[I: f_p]                  
	   const double    fc,                      //[I: f_c]                  
#ifdef ALLOW_INDI
	   const double    ns,                      //[I: Nsoft]                
	   const int       nr,                      //[I: Nref]                 
	   const int       sf,                      //[I: softening: global]    
#endif
	   const double    gr,                      //[I: grav]                 
	   const bool      resume,                  //[I: resume old sim (nemo)]
	   const extpot   *pt,                      //[I: external potential]   
	   const io        read_more,               //[I: what else to read?]   
	   const int       nd,                      //[I: # nemo output streams]
	   const int       dir[4]   ) :             //[I: direct sum control]   
  MYNDATA( new ndata(in,nemo,th,hg,nc,ep,kn,
#ifdef ALLOW_INDI
		     sf,ns,nr,
#endif
		     hm,nl,fa,fp,fc,pt,resume,read_more,nd,gr)),
  MYNBODY( MYNDATA->NSTEPS > 1 ?
	   static_cast<basic_nbody*>(
	     new BlockStepCode(MYNDATA->BODIES,
			       MYNDATA->EPS,
			       MYNDATA->TIME,
			       MYNDATA->HMIN+1-MYNDATA->NSTEPS,
			       MYNDATA->NSTEPS,
			       MYNDATA->FACC,
			       MYNDATA->FPOT,
			       MYNDATA->FCOM,
			       MYNDATA->HGROW,
			       MYNDATA->THETA,
			       MYNDATA->NCRIT,
			       kern(MYNDATA->KERN),
#ifdef ALLOW_INDI
			       soften(MYNDATA->SOFTEN),
			       MYNDATA->NSOFT,
			       MYNDATA->NREF,
			       two,
#endif
			       MYNDATA->PEX,
			       dir) )
	   :
	   static_cast<basic_nbody*>(
	     new LeapFrogCode( MYNDATA->BODIES,
			       MYNDATA->EPS,
			       MYNDATA->TIME,
			       MYNDATA->HMIN,
			       MYNDATA->HGROW,
			       MYNDATA->THETA,
			       MYNDATA->NCRIT,
			       kern(MYNDATA->KERN),
#ifdef ALLOW_INDI
			       soften(MYNDATA->SOFTEN),
			       MYNDATA->NSOFT,
			       MYNDATA->NREF,
			       two,
#endif
			       MYNDATA->PEX,
			       dir) )
	   )
{
  MemoryCheck(MYNDATA);
  MemoryCheck(MYNBODY);
  if(resume) MYNBODY->reset_cpu_total(MYNDATA->CPU_I);
}
#endif
//------------------------------------------------------------------------------
yanc::yanc(const char  *in,                         // I: input file            
	   pot_provider*pv) :                       //[I: provider for ext pot] 
  MYNDATA( new ndata(in,pv) ),
  MYNBODY( MYNDATA->NSTEPS > 1 ?
	   static_cast<basic_nbody*>(
	     new BlockStepCode(MYNDATA->BODIES,
			       MYNDATA->EPS,
			       MYNDATA->TIME,
			       MYNDATA->HMIN+1-MYNDATA->NSTEPS,
			       MYNDATA->NSTEPS,
			       MYNDATA->FACC,
			       MYNDATA->FPOT,
			       MYNDATA->FCOM,
			       MYNDATA->HGROW,
			       MYNDATA->THETA,
			       MYNDATA->NCRIT,
			       kern(MYNDATA->KERN),
#ifdef ALLOW_INDI
			       soften(MYNDATA->SOFTEN),
			       MYNDATA->NSOFT,
			       MYNDATA->NREF,
			       two,
#endif
			       MYNDATA->PEX,
			       Default::direct) )
	   :
	   static_cast<basic_nbody*>(
	     new LeapFrogCode( MYNDATA->BODIES,
			       MYNDATA->EPS,
			       MYNDATA->TIME,
			       MYNDATA->HMIN,
			       MYNDATA->HGROW,
			       MYNDATA->THETA,
			       MYNDATA->NCRIT,
			       kern(MYNDATA->KERN),
#ifdef ALLOW_INDI
			       soften(MYNDATA->SOFTEN),
			       MYNDATA->NSOFT,
			       MYNDATA->NREF,
			       two,
#endif
			       MYNDATA->PEX,
			       Default::direct) )
	   )
{
  MemoryCheck(MYNDATA);
  MemoryCheck(MYNBODY);
}
//------------------------------------------------------------------------------
yanc::~yanc()
{
  delete MYNBODY;
  delete MYNDATA;
}
//------------------------------------------------------------------------------
const bodies* yanc::mybodies() const
{
  return MYNDATA->BODIES;
}
//------------------------------------------------------------------------------
