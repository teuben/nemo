//
// nbody.C
//
// Version 1999/2/9 Jun Makino
//
// Change in the calling format for apply_vf, from function name to
// add explicit address operator. This was necessary to pass g++ 2.8.1
// 
// Version 1.1 Jun Makino : take parameters from command line argumensts
//                          Coding for 1.1 started on 1998/12/31
// Version 1.0 Jun Makino : First publisized version, parameters
//                          read in from the standard input

#ifndef NOGRAPHICS
#ifndef FINAL_PLOT
#define GRAPHICS
#endif
#endif


#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#include  <cstdlib>
#include  <math.h>
#include  <iostream>
#include  <fstream>

using namespace std;

#ifdef GRAPE5
#include <gp5util.h>
#endif /* GRAPE5 */

#ifndef NOGRAPHICS
#include <cpgplot.h>
#endif
#define real double
#include "vector.h"
#include "nbody_particle.h"
#define NBODY


typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef nbody_VF_ptr real_VF_ptr;
typedef nbody_RF_ptr real_RF_ptr;
typedef nbody_RRF_ptr real_RRF_ptr;

extern "C" double cpusec();
int  pgetopt(int argc, char ** argv,  char * optstr);
void pskipopt();
void pnextarg(char **argv);



static real_system pb;
static int snapin_flag = 0;
static int snapout = 0;
static real  dt = 0.015625;
static real  dtsnapout = 1;
static int outlogstep = 1;
static real tend = 10;
static real eps = 0.025;
static real theta = 0.75;
static int ncrit = 1024;
static int collision_flag = 0;
static int relx_flag = 0;
static int relv_flag = 0;
static vector relpos, relv;
static real pos_scale = 1;
static real vel_scale = 1;
static char finame[255];
static char foname[255];
static int snapout_acc_flag = 0;
static int snapio_binary_flag = 0;
static int command_exec_flag = 0;
static char command_name[255];
static real mmin = 0.0;
static real xmin = 0.0, xmax = 0.0;
static int autoscale = 1;
static int node_div_crit = 8;

#ifdef GB99
static real dt0 = 0.002;
static real eps0 = 0.025;
#endif // GB99

#ifdef P2M2
static int full_dof = 1;
static int me_order = 1;
static real ppscale = 0.2;
#endif // P2M2

#ifdef TESTCODE
static    vector evalpos = vector(1.0, 0.0, 0.0);
#endif // TESTCODE


void real_particle::read(istream& ss)
{
    ss 	>> index ;
    real dum;
    ss >>pos >>vel >> mass>> acc_gravity >> phi_gravity >> dum;
}

void real_system::read(istream& s)
{
    cerr << "enter read: ";PRL(nsize);
    s>>n >> eps2_for_gravity >> use_self_gravity >> theta_for_tree >> ncrit_for_tree;
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }

    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }
    for(int i = 0; i<n; i++){
	(pb+i)->read(s);
    }
    cerr << "exit  read: ";PRL(nsize);
}

#ifdef IONEMO

extern "C" void io_nemo(char *, char *, ...);

static int *ionemo_np = NULL;
static real *ionemo_tp = NULL;
static real *ionemo_mp = NULL;
static real (*ionemo_xp)[3] = NULL;
static real (*ionemo_vp)[3] = NULL;
static real *ionemo_phip = NULL;
static real (*ionemo_accp)[3] = NULL;

void real_system::btos(char *fname)
{
    io_nemo(fname, "double, n, t, x, v, m, read, info",
	    &ionemo_np, &ionemo_tp, &ionemo_xp, &ionemo_vp, &ionemo_mp);
    io_nemo(fname,"close");

    n = *ionemo_np;
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }
    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }

    if (ionemo_tp == NULL) {
	cerr << "btos: no time tag in input file " << fname << endl;
	cerr << "set time to 0.0" << endl;
	time = time0 = 0.0;
	ionemo_tp = new real[1];
    }
    else {
	time0 = *ionemo_tp;
    }
    for(int i = 0; i < n; i++){
	vector vtmp;
	(pb+i)->set_mass(ionemo_mp[i]);
	(pb+i)->set_index(i);
	vtmp = vector(ionemo_xp[i][0], ionemo_xp[i][1], ionemo_xp[i][2]);
	(pb+i)->set_pos(vtmp);
	vtmp = vector(ionemo_vp[i][0], ionemo_vp[i][1], ionemo_vp[i][2]);
	(pb+i)->set_vel(vtmp);
    }
    if (!snapout)
    {
	delete ionemo_mp;
	delete ionemo_xp;
	delete ionemo_vp;
    }
}

void real_system::stob(char *fo, char *fi)
{
    *ionemo_np = n;
    *ionemo_tp = time;

    for (int i = 0; i < n; i++) {
	ionemo_mp[i] = (pb+i)->get_mass();
    }
    for (int i = 0; i < n; i++) {
	vector vtmp = (pb+i)->get_pos();
	for (int k = 0; k < 3; k++) {
	    ionemo_xp[i][k] = vtmp[k];
	}
    }
    for (int i = 0; i < n; i++) {
	vector vtmp = (pb+i)->get_vel();
	for (int k = 0; k < 3; k++) {
	    ionemo_vp[i][k] = vtmp[k];
	}
    }
    for (int i = 0; i < n; i++) {
	if (ionemo_phip == NULL) {
	    ionemo_phip = new real[n];
	}
	ionemo_phip[i] = (pb+i)->get_phi_gravity();
    }

    if (snapout_acc_flag == 1) {
	if (ionemo_accp == NULL) {
	    ionemo_accp = new real[n][3];
	}
	for (int i = 0; i < n; i++) {
	    vector vtmp = (pb+i)->get_acc_gravity();
	    for (int k = 0; k < 3; k++) {
		ionemo_accp[i][k] = vtmp[k];
	    }
	}
    }

    if (snapout_acc_flag == 1) {
	io_nemo(fo, "double, n, t, x, v, m, p, a, h, save, info",
		&ionemo_np, &ionemo_tp, &ionemo_xp, &ionemo_vp,
		&ionemo_mp, &ionemo_phip, &ionemo_accp, fi);
    }
    else {
	io_nemo(fo, "double, n, t, x, v, m, p, h, save, info",
		&ionemo_np, &ionemo_tp, &ionemo_xp, &ionemo_vp,
		&ionemo_mp, &ionemo_phip, fi);
    }
    io_nemo(fo,"close");
}


#endif // IONEMO

void real_system::atos(istream& s)
{
    cerr << "enter atos: ";PRL(nsize);
    int nddum;
    s>>n >> nddum >> time0;
	
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }

    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }
    for(int i = 0; i<n; i++){
	real mtmp;
	s>>mtmp;
	(pb+i)->set_mass(mtmp);
	(pb+i)->set_index(i);
	
    }
    for(int i = 0; i<n; i++){
	vector ptmp;
	s>>ptmp;
	(pb+i)->set_pos(ptmp);
    }
    for(int i = 0; i<n; i++){
	vector vtmp;
	s>>vtmp;
	(pb+i)->set_vel(vtmp);
    }
    cerr << "exit  atos: ";PRL(nsize);
}


void real_system::stoa(ostream& ss)
{
    cerr << "enter  stoa" << endl;

    ss.precision(12);
    ss <<n <<endl;
    ss << "3\n";
    ss <<  time << endl;
	
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_mass() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_pos() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_vel() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_phi_gravity() <<endl;}
    if (snapout_acc_flag == 1) {
	for(int i = 0; i<n; i++){ss<<(pb+i)->get_acc_gravity() <<endl;}
    }

    cerr << "exit  stoa" << endl;
}

void real_particle::write(ostream& ss)
{
    ss << index << " " ;
    ss << pos << " "<< vel << " "<< mass << " "  <<acc_gravity << " " <<  phi_gravity <<endl;
}

void real_system::write(ostream& s)
{
    s<<n << " " << eps2_for_gravity <<" " << use_self_gravity
     << " " << theta_for_tree << " " <<ncrit_for_tree <<endl;;
    for(int i = 0; i<n; i++){
	(pb+i)->write(s);
    }
}



void real_particle::dump()
{
    PRL(index);
    PRC(pos);  PRC(vel);    PRL(acc_gravity);
    PRC(mass);  PRL(phi_gravity);
}

void real_system::dump()
{
    PRL(n);
    for(int i = 0; i<n; i++){
	(pb+i)->dump();
    }
}

const real piinv = 1.0/M_PI;
real real_particle::kinetic_energy()
{
    return (0.5*vel*vel)*mass;
}

real real_particle::energy()
{
    return kinetic_energy();
}


real real_system::kinetic_energy()
{
    real sum = 0;
    for(int i = 0; i<n; i++){
	sum+= (pb+i)->kinetic_energy();
    }
    return sum;
}




real real_system::energy()
{
    real sum = kinetic_energy();
    if(use_self_gravity){
	for(int i = 0; i<n; i++){
	    sum+= (pb+i)->get_phi_gravity()*(pb+i)->get_mass()
#ifdef REAL_GRAVITY		
		*0.5;
#else
	    ;
#endif	    
	}
    }
    return sum;
}


void real_system::apply_vf(real_VF_ptr f)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)();
    }
}
	
void real_system::apply_vf(real_RF_ptr f, real r)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)(r);
    }
}
	
void real_system::apply_vf(real_RRF_ptr f, real r, real r2)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)(r, r2);
    }
}


void real_system::integrate(real dt)
{
    apply_vf(&real_particle::predict,dt);
    calculate_gravity();
    apply_vf(&real_particle::correct,dt);
}    
    
const vector uniform_random_position_in_sphere(double power_index)
{
    vector x;
    do{
	for(int i = 0; i<3;i++) x[i] = drand48()*2-1;
    }while (x*x >= 1);
    x *=  pow(x*x, 3.0/(power_index+3)-1);
    return x;
}

void real_system::calculate_cmterms()
{
    pos = 0.0;
    vel = 0.0;
    mass = 0.0;
    real_particle * p = pb;
    for(int i = 0; i<n;i++){
	pos += p->get_pos()*p->get_mass();
	vel += p->get_vel()*p->get_mass();
	mass += p->get_mass();
	p++;
    }
//    cerr << "CM ";PRC(pos/mass); PRL(vel/mass);
    cout << "CM "<< (pos/mass) << "\nCMV " <<(vel/mass)<<endl;
}
    

void real_system::create_uniform_sphere(int nbody, real power_index, real r0)
{
    PRC(nbody); PRL(power_index);
    n = nbody;
    pb = new real_particle[n];
    nsize = n;
    real_particle * p = pb;
    for(int i = 0; i<n;i++){
	p->set_pos(uniform_random_position_in_sphere(power_index)*r0);
	p->set_vel(0);
	p->set_mass(1.0/n);
	p->set_index(i);
	p++;
    }
}

void copy_nbody_particles(real_system * source,
			  real_system * destination)
{

    destination->pb = new real_particle[source->nsize];
    real_particle * pbs = source->pb;
    real_particle * pbd = destination->pb;
    for(int i = 0; i<source->n;i++){
	*pbd = *pbs; 
	pbd++;
	pbs++;
    }
}

int compare_particle_pointer(real_particle * * p1, real_particle * *p2)
{
    int i1 = (*p1)->get_index();
    int i2 = (*p2)->get_index();
    if (i1 > i2){
	return 1;
    }else if (i1 == i2){
	return 0;
    }else{
	return -1;
    }
}

void initgraph()
{
#ifdef GRAPHICS    
    if(cpgopen("/XWINDOW") != 1) exit(EXIT_FAILURE);
    cpgask(0);
    cpgvsiz(1.0,6.0,1.0,6.0);
#endif    
}

void real_particle::plot(real parm)
{
#ifdef GRAPHICS    
    int iparm = (int)parm;
    cpgpt1(pos[0],pos[1],-1);
#endif
}
void real_system::plot(real t)
{
#ifdef GRAPHICS    
    real xmax = plot_xmax;
    cpgbbuf();
    cpgsci(1);
    cpgslw(2);
    cpgenv(-xmax, xmax, -xmax, xmax,  1, 0);
    char label[255];
    sprintf(label, "T = %g N=%d", t,n);
    cpglab("x","y", label);
    cpgsci(1);
    int psize = 4;
    if (n < 1000){
	psize *= (int)pow(1000.0/n,0.25);
    }
    cpgslw(psize);
    apply_vf(&real_particle::plot, 0.0);
    cpgebuf();
#endif
}

void real_system::make_collision(vector relpos, vector relv)
{
    real_particle *new_pb;
    new_pb = new real_particle[n*2];
    int i;
    for( i=0; i<n;i++){
	*(new_pb+i) = *(pb+i);
	*(new_pb+i+n) = *(pb+i);
    }
    delete [] pb;
    pb = new_pb;
    for( i=0; i<n;i++){
	(pb+i)->inc_pos(0.5*relpos);
	(pb+i)->inc_vel(0.5*relv);
	(pb+i)->set_index(i);
    }
    for( i=n; i<n*2;i++){
	(pb+i)->inc_pos(-0.5*relpos);
	(pb+i)->inc_vel(-0.5*relv);
	(pb+i)->set_index(i);
    }
    n *= 2;
    nsize *= 2;
    PR(n); PRL(nsize);
}

void
parse_flags(int argc, char ** argv)
{
    extern char *poptarg;
    int c;

#ifdef GRAPE5
#define SCALEOPTS "M:R:"
#else // no GRAPE5
#define SCALEOPTS
#endif // GRAPE5 ""

#if defined(GB99)
    char* param_string = "abi:o:d:D:l:T:e:E:f:t:n:N:w:cx:v:s:S:h"SCALEOPTS;
#elif defined(P2M2)
    char* param_string = "abi:o:d:D:l:T:e:t:n:N:w:cx:v:s:S:p:r:q:hFX:"SCALEOPTS;
#else
    char* param_string = "abi:o:d:D:l:T:e:t:n:N:w:cx:v:s:S:hX:"SCALEOPTS;
#endif

#undef SCALEOPTS


    while ((c = pgetopt(argc, argv, param_string)) != -1){
        switch(c) {
#ifdef GRAPE5
	case 'M': mmin = atof(poptarg);
	    break;
	case 'R': autoscale = 0;
	    xmax = atof(poptarg);
	    xmin = -xmax;
	    break;
#endif // GRAPE5
	case 'a': snapout_acc_flag = 1;
	    break;
#ifdef IONEMO
	case 'b': snapio_binary_flag = 1;
	    break;
#endif // IONEMO
	case 'i': strcpy(finame,poptarg);
	    snapin_flag = 1;
	    break;
	case 'o': strcpy(foname,poptarg);
	    snapout = 1;
	    break;
	case 'd': dt = atof(poptarg);
	    break;
	case 'D': dtsnapout = atof(poptarg);
	    break;
	case 'l': outlogstep = atoi(poptarg);
	    break;
	case 'T': tend = atof(poptarg);
	    break;
	case 'e': eps = atof(poptarg);
	    break;
#ifdef GB99
	case 'E': eps0 = atof(poptarg);
	    break;
	case 'f': dt0 = atof(poptarg);
	    break;
#endif
	case 't': theta = atof(poptarg);
	    break;
	case 'n': ncrit = atoi(poptarg);
	    break;
	case 'N': node_div_crit = atoi(poptarg);
	    break;
	case 'w': pb.plot_xmax = atof(poptarg);
	    break;
	case 'c': collision_flag = 1;
	    break;
	case 'x': relx_flag = 1;
	    relpos = vector(atof(poptarg),
			    atof(poptarg+1),
			    atof(poptarg+2));
	    pskipopt();pskipopt();
	    break;
	case 'v': relv_flag = 1;
	    relv = vector(atof(poptarg),
			  atof(poptarg+1),
			  atof(poptarg+2));
	    pskipopt();pskipopt();
	    break;
	case 's': pos_scale = atof(poptarg);
	    break;
	case 'S': vel_scale = atof(poptarg);
	    break;
	case 'X': strcpy(command_name,poptarg);
	    command_exec_flag = 1;
	    break;
#ifdef P2M2
	case 'F': full_dof = 0;
	    break;
	case 'p': me_order = atoi(poptarg);
	    break;
	case 'r': ppscale = atof(poptarg);
	    break;
#endif // P2M2
#ifdef TESTCODE
	case 'q': 
	    real rr, tt, pp;

	    rr = atof(poptarg);
	    pnextarg(argv);
	    tt = atof(poptarg);
	    pnextarg(argv);
	    pp = atof(poptarg);
	    tt = tt/180.0*M_PI;
	    pp = pp/180.0*M_PI;
//	    PR(rr);PR(tt);PRL(pp);
	    evalpos = vector(rr*sin(tt)*cos(pp),
			     rr*sin(tt)*sin(pp),
			     rr*cos(tt));

	    break;
#endif
	case 'h':		      
	    cerr << "list of options\n";
	    cerr << "-a        output accelaration to snapshot output file \n";
#ifdef IONEMO
	    cerr << "-b        assume NEMO structured binary format for snapshot I/O\n";
#endif // IONEMO
	    cerr << "-i        name of snapshot input file       (no default)\n";
	    cerr << "-o        name of snapshot output file      (default: no output)\n";
	    cerr << "-D        time interval for snapshot output (default: 1)\n";
	    cerr << "-l        interval for log output (default: 1: all step)\n";
	    cerr << "-T        time to stop integration          (default: 10)\n";
#ifdef GB99
	    cerr << "-f        initial timestep\n";
	    cerr << "-d        final timestep\n";
	    cerr << "-E        initial softening parameter\n";
	    cerr << "-e        final softening parameter\n";
#else
	    cerr << "-d        timestep (default: 0.015625)\n";
	    cerr << "-e        softening parameter (default: 0.025)\n";
#endif
	    cerr << "-t        opening angle theta               (default: 0.75)\n";
	    cerr << "-n        ncrit for Barnes' vectrization    (default: 1024)\n";
	    cerr << "          (ONLY used with GRAPE/HARP implementation)\n";
	    cerr << "-N        cell division criterion    (default: 8 or 8 * # of pseudo particle)\n";
	    cerr << "          divide cell until it contains >N particles\n";
	    cerr << "          smaller memory requirement for larger N\n";
	    cerr << "-w        window size for PGPLOT snapshot plot (default: 10)\n";
	    cerr << "-c        flag for collision run\n";
	    cerr << "-x        relative position vector for collision run (no default)\n";
	    cerr << "-v        relative velocity vector for collision run (no default)\n";
	    cerr << "-s        scale factor for position scaling (default: 1)\n";
	    cerr << "-S        scale factor for velocity scaling (default: 1)\n";
#ifdef P2M2
	    cerr << "-p        order of multipole expansion (default: 1)\n";
	    cerr << "          (ONLY used with P2M2 implementation)\n";
	    cerr << "-r        sphere radius on which pseudo particles are distributed (default: 0.2)\n";
	    cerr << "          (ONLY used with P2M2 implementation)\n";
	    cerr << "-F        never use full degree of freedom of pseudo particle even when p < 3\n";
	    cerr << "          (ONLY used with P2M2 implementation)\n";
#endif // P2M2
#ifdef GRAPE5
	    cerr << "-M        set min mass to be resolved (default: mass of min mass particle)\n";
	    cerr << "-R        set the min and max value of coordinate to -R and R, respectively\n";
	    cerr << "          (default: the corrdinate is automatically rescaled every timestep\n";
	    cerr << "           so that the root cell is included in the range of [xmin,xmax])\n";
	    cerr << "          (ONLY used with GRAPE-5 implementation)\n";
#endif // GRAPE5
#ifdef TESTCODE
	    cerr << "-q        position(r, theta[degree], phi[degree]) at which force is evaluated \n";
#endif // TESTCODE
	    cerr << "-X        command to be executed each time snapshot are output (default: no output) \n";
	    cerr << "          snapshot file name and its ID (counting from 0) is passed to the command.\n";
	    cerr << "-h        print this help\n";
	    exit(1);
	    break;
	default:
	    cerr << "unknown option: -" << *(char *)&c << " abort." << endl;
	    exit(1);
	}
    }
}


#ifndef TESTCODE

int
main(int argc, char ** argv)
{
    ifstream fsnapin;
    ofstream fsnapout;
    real  tsnapout = 0;
    foname[0] = '?';
    foname[1] = '\0';
    pb.plot_xmax = 10;

    parse_flags(argc, argv);
    if (snapin_flag == 0){
	cerr << "Snapshot input file required (-i)\n";
	exit(1);
    }

#ifdef PAR_TRAV
    extern void create_force_calculator(void);
    create_force_calculator();
#endif


    PR(dt); PR(dtsnapout); PRL(tend);
    cerr << "snapin = " << finame << " snapout = " << foname <<endl;

#ifdef P2M2
    pb.load_design(me_order, ppscale, full_dof);
    cerr << "expansion_order = " << pb.p2m2.order
	 << " spherical_design_order = " << pb.p2m2.spherical_design
	 << " num_of_pseudo_particle = " << pb.p2m2.npp
	 << " sphere_scale = " << pb.p2m2.ppscale
	 <<endl;
    cerr.precision(17);
//    cerr.setf(ios::scientific, ios::floatfield);
    for (int i = 0; i < pb.p2m2.npp; i++) {
	cerr << pb.p2m2.pppos0[i][0] << " "
	     << pb.p2m2.pppos0[i][1] << " "
	     << pb.p2m2.pppos0[i][2] << endl;
    }
    cerr.precision(6);
#endif // P2M2

#ifdef IONEMO
    if (snapio_binary_flag) {
	pb.btos(finame);
    }
    else
#endif // IONEMO
    {
	fsnapin.open(finame,ios::in);
	pb.atos(fsnapin);
    }
    fsnapin.close();
    if (pb.time0 > tend) {
	cerr << "start time " << pb.time0 << " larger than end time"
	     << tend << ". abort." << endl << endl;
	exit(1);
    }
    cerr << "n= " << pb.n << endl;
    pb.use_self_gravity = 1;

    PRL(node_div_crit);
    cerr << "eps= " << eps << " theta=" << theta << " ncrit=" <<ncrit <<endl;
#ifdef GB99
    cerr << "eps0= " << eps0 << endl;
    cerr << "dt0= " << dt0 << endl;
#endif
    pb.node_div_crit = node_div_crit;
    pb.eps2_for_gravity = eps*eps;
    pb.theta_for_tree = theta;
    pb.ncrit_for_tree = ncrit;
    if(collision_flag == 1){
	if (relx_flag == 0){
	    cerr << "relative position required (x option)\n";
	    exit(1);
	}
	if (relv_flag == 0){
	    cerr << "relative velocity required (v option)\n";
	    exit(1);
	}
	pb.make_collision(relpos, relv);
	cerr << "Collision relp= " << relpos << " relv=" << relv <<endl;
    }else{
	pb.apply_vf(&real_particle::scale_pos,pos_scale);
	pb.apply_vf(&real_particle::scale_vel,vel_scale);
	cerr << "posscale= " << pos_scale << " velscale=" << vel_scale <<endl;
    }

    cerr << endl;
#ifdef GRAPHICS    
    initgraph();
#endif

#ifdef GRAPE5
    g5_open();
    pb.set_autoscale(autoscale);

    real m0, x0, x1;
    if (mmin != 0.0) {
	pb.set_massmin(mmin);
	m0 = mmin;
    }
    else {
	pb.set_massmin();
	m0 = pb.get_massmin();
    }
    if (xmin != 0.0) {
	x0 = xmin;
	x1 = xmax;
    }
    else {
	x0 = -32.0;
	x1 = 32.0;
    }

    g5_set_range(x0, x1, m0);
    PRL(pb.get_autoscale());
    cerr << "g5_set_range("
	 << x0 << ", "
	 << x1 << ", "
	 << m0 << ")" << endl;
    g5_close();
#endif /* GRAPE5 */

#ifdef GB99 // variable eps & dt

    real zval = 24.0; // red shift at T = 0
    real tmpzz = pow(1.0/(zval+1), 3.0/2.0);
    real t0 = tmpzz/(1.0-tmpzz)*tend;
    real eps1 = (zval+1)*pow((pb.time+t0)/(tend+t0), 2.0/3.0)*eps0;
    pb.eps2_for_gravity = eps1 * eps1;
    cerr << "current eps= " << sqrt(pb.eps2_for_gravity) << endl;

    real dt1 = dt0;
    tsnapout = dtsnapout;
    cerr << "initial calculate_garvity() start at " << cpusec() << endl;
    pb.calculate_gravity();
    cerr << "initial calculate_garvity() end at " << cpusec() << endl;
    real E0 = pb.energy();
    PRL(E0);
    pb.plot(0.0);

    int step = 0;
    while (pb.time < tend) {

	if (pb.eps2_for_gravity < eps * eps)
	{
	    real eps1 = (zval+1)*pow((pb.time+t0)/(tend+t0), 2.0/3.0)*eps0;
	    pb.eps2_for_gravity = eps1 < eps ? eps1 * eps1: eps * eps;
	}
	else
	{
	    pb.eps2_for_gravity = eps * eps;
	}
	cerr << "current eps: " << sqrt(pb.eps2_for_gravity) << endl;

	if (dt1 < dt)
	{
	    dt1 = dt0*(pb.time+t0)/t0;
	}
	else
	{
	    dt1 = dt;
	}
	pb.time += dt1;
	pb.integrate(dt1);
	cerr << "pb.time: " << pb.time << " dt1: " << dt1 << " step: " << step << endl;

 	if (step%outlogstep == outlogstep - 1){
 	    real KE = pb.kinetic_energy();
 	    real E = pb.energy();
 	    real Eerr = fabs((E0-E)/E0);
 	    real Q = KE/(KE-E);
 	    real T = pb.time;
// 	    PR(T); PR(KE); PR(Q); PR(E); PRL(Eerr);
 	    cerr << "T= " << T << " KE= " << KE << " Q= " << Q
		 << " E= " << E << " Eerr= " << Eerr
		 << " Eabserr= " << fabs(E0-E) << endl;
 	    pb.plot(pb.time);
 	    pb.calculate_cmterms();
 	}
	if (tsnapout < pb.time){
	    if (snapout)
	    {
		static int nout = 0;
		ofstream fs;
 		char fn[255];
		char cmd[255];

 		sprintf(fn, foname, nout);

#ifdef IONEMO
		if (snapio_binary_flag) {
		    pb.stob(fn, finame);
		    cerr << "no binary I/O supported. abort." << endl;
		    exit(1);
		}
		else
#endif // IONEMO
		{
		    fs.open(fn, ios::out); 
		    pb.stoa(fs);
		    fs.close();
		}
		nout++;

		cerr << "time: " << pb.time << " out: " << fn << endl;
	    }
	    tsnapout += dtsnapout;
	}
	cerr << "CPU sec = " <<cpusec() << endl <<endl;
	step++;
    }

#else // constant eps & dt

    int nstep = (int)((tend-pb.time0)/dt+0.1);
    dt = (tend-pb.time0)/nstep;
    int outsnapstep = (int) (dtsnapout/dt+0.1);
    cerr << "initial calculate_garvity() start at " << cpusec() << endl;
    pb.calculate_gravity();
    cerr << "initial calculate_garvity() end at " << cpusec() << endl;
    real E0 = pb.energy();
    PRL(E0);
    pb.plot(0.0);

    // calculate force at T=0 and quit
    if (snapout && dtsnapout == 0.0)
    {
	snapout_acc_flag = 1;
	static int nout = 0;
	ofstream fs;
	char fn[255];

	sprintf(fn, foname, nout);

#ifdef IONEMO
	if (snapio_binary_flag) {
	    pb.stob(fn, finame);
	}
	else
#endif // IONEMO
	{
	    fs.open(fn, ios::out); 
	    pb.stoa(fs);
	    fs.close();
	}

	cerr << "out: " << fn << endl;
	cerr << "CPU sec = " <<cpusec() << endl << endl;
	exit(1);
    }

    for(int i=0;i<nstep;i++){

	pb.time = pb.time0+(i+1)*dt;
	pb.integrate(dt);

 	cerr << "Exit integrate, cpu = " <<cpusec() << endl;
 	if ( i%outlogstep == outlogstep - 1){
 	    real KE = pb.kinetic_energy();
 	    real E = pb.energy();
 	    real Eerr = ((E0-E)/E0);
 	    real Q = KE/(KE-E);
 	    real T = pb.time;
// 	    PR(T); PR(KE); PR(Q); PR(E); PRL(Eerr);
 	    cout << "T= " << T << " KE= " << KE << " Q= " << Q
		 << " E= " << E << " Eerr= " << Eerr
		 << " Eabserr= " << (E0-E) << endl;
 	    pb.plot(pb.time);
 	    pb.calculate_cmterms();
 	}
	if (i % outsnapstep == outsnapstep - 1){
	    if (snapout)
	    {
		static int nout = 0;
		ofstream fs;
 		char fn[255];
		char cmd[255];

 		sprintf(fn, foname, nout);

#ifdef IONEMO
		if (snapio_binary_flag) {
		    pb.stob(fn, finame);
		}
		else
#endif // IONEMO
		{
		    fs.open(fn, ios::out); 
		    pb.stoa(fs);
		    fs.close();
		}

		if (command_exec_flag) {
		    sprintf(cmd, "%s %s %03d", command_name, fn, nout);
		    cerr << "execute " << cmd << endl << endl;
		    system(cmd);
		}
		nout++;
		cerr << "time: " << pb.time << " out: " << fn << endl;
	    }
	}
	cerr << "CPU sec = " <<cpusec() << endl <<endl;
    } // i loop

#endif // GB99
    return 0;
}

#else // TESTCODE

#include "BHtree.h"
bhnode * get_bn(void);
#ifdef P2M2
void set_cm_quantities_for_default_tree(P2m2 *p2m2);
#else // no P2M2
void set_cm_quantities_for_default_tree();
#endif // P2M2
void clear_tree_counters();
void print_tree_counters();

void real_system::load_homo_particle(int n)
{
    cerr << "enter load_home_particle: "; PR(n); PRL(nsize);

    real mtmp = 1.0/n;

    this->n = n;
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }
    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }
    for(int i = 0; i<n; i++){
	vector ptmp;
	ptmp = vector(drand48()-0.5,
		      drand48()-0.5,
//		      drand48()*0.5-0.25);
		      drand48()-0.5);
	(pb+i)->set_mass(mtmp);
	(pb+i)->set_pos(ptmp);
    }
}

void real_system::calculate_gravity_direct_at(vector evalpos, vector &acc, real &phi)
{
    acc = 0.0;
    phi = 0.0;
    for (int i = 0; i < n; i++) {
	vector pos = (pb+i)->get_pos();
	real mass = (pb+i)->get_mass();
	vector dr = evalpos-pos;
	real r = sqrt(square(dr)+eps2_for_gravity);

	acc += -mass/r/r/r*dr;
	phi += -mass/r;
    }
}


#include "p2m2.h"

void real_system::calculate_gravity_classicalME_at(vector evalpos, real &phi)
{
    phi = 0.0;
    for (int i = 0; i < n; i++) {
	vector pos = (pb+i)->get_pos();
	real mi = (pb+i)->get_mass();
	real ri = abs(pos);
	real rj = abs(evalpos);
	real cosg = pos*evalpos/ri/rj;
	real pln[16];
	plgndr0(p2m2.order+1, cosg, pln);
	for (int l = 0; l <= p2m2.order; l++) {
	    phi += -mi/rj*pow(ri/rj, l)*pln[l];
	}
    }
}

void real_system::calculate_gravity_at(vector evalpos, vector &acc, real &phi)
{
    cerr << "Enter setup_tree, cpu = " << cpusec() << endl;    
    setup_tree(&p2m2);
    cerr << "Enter set_cm_quantities_for_default_tree, cpu = " << cpusec() << endl;    
    set_cm_quantities_for_default_tree(&p2m2);
    clear_tree_counters();

    apply_vf(&(real_particle::clear_acc_phi_gravity));

    acc = 0.0;
    phi = 0.0;

    bhnode *bnlocal;
    real largetheta2 = 1e8;

    cerr << "Enter accumulate_force_from_tree, cpu = " << cpusec() << endl;    
    bnlocal = get_bn();
    bnlocal->accumulate_force_from_tree(evalpos, eps2_for_gravity, largetheta2,
					acc, phi ,&p2m2);
    
    cerr << "Exit evaluate_gravity, cpu = " << cpusec() << endl;    
}


main(int argc, char ** argv)
{
    ifstream fsnapin;
    ofstream fsnapout;
    real  tsnapout = 0;
    foname[0] = '?';
    foname[1] = '\0';
    pb.plot_xmax = 10;
    parse_flags(argc, argv);

    cerr << "---- Test code ----" << endl << endl;
    cerr << "snapin = " << finame << " snapout = " << foname <<endl;

#ifdef P2M2
    pb.load_design(me_order, ppscale, full_dof);
    cerr << "expansion_order = " << pb.p2m2.order
	 << " spherical_design_order = " << pb.p2m2.spherical_design
	 << " num_of_pseudo_particle = " << pb.p2m2.npp
	 << " sphere_scale = " << pb.p2m2.ppscale
	 <<endl;
#if 0
    cerr.precision(17);
    for (int i = 0; i < pb.p2m2.npp; i++) {
	cerr << pb.p2m2.pppos0[i][0] << " "
	     << pb.p2m2.pppos0[i][1] << " "
	     << pb.p2m2.pppos0[i][2] << endl;
    }
    cerr.precision(6);
#endif
#endif // P2M2

#ifdef IONEMO
    if (snapio_binary_flag) {
	pb.btos(finame);
    }
    else
#endif // IONEMO
    {
	fsnapin.open(finame,ios::in);
	pb.atos(fsnapin);
	fsnapin.close();
    }
    cerr << "n= " << pb.n << endl;
    pb.use_self_gravity = 1;

    cerr << "eps= " << eps << " theta=" << theta << " ncrit=" <<ncrit <<endl;
    pb.eps2_for_gravity = eps*eps;
    pb.theta_for_tree = theta;
    pb.ncrit_for_tree = ncrit;

    vector acc, acc0;
    real phi, phi0, phi1;
    real l;
    cerr.precision(22);

    ofstream fs;
    if (snapout)
    {
	fs.open(foname, ios::out); 
    }

#define ERR_DIST (1)
#define ERR_RAD (0)
#define ERR_N (0)

#if (ERR_DIST+ERR_RAD+ERR_N) != 1
	incorrect directive
#endif


#if ERR_DIST // err vs distance of evaluation point
    for (real q = 0.001; q < 25.0; q *= 1.2) {
	cerr << "distance: " << q << " p: " << me_order << endl;
#elif ERR_RAD // err vs pp-sphere radius
    for (real s = 0.01; s < 2.0; s *= 1.2) {
	real q = 1.0;
	pb.load_design(me_order, s, full_dof);
	cerr << "rad: " << s << " p: " << me_order << endl;
#elif ERR_N
//    for (int n = 2; n <= 32768; n *= 2) {
    for (int n = 2; n <= 8192; n *= 2) {
	real q = 1.0;
	int ntry = 10;
	real errsum = 0, perrsum = 0;
	real rerrsum = 0, rperrsum = 0;
	for (int i = 0; i < ntry; i++) {
	    cerr << "try: " << i << " p: " << me_order << endl;
	pb.load_homo_particle(n);
#else
	incorrect directive
#endif
	pb.calculate_gravity_direct_at(evalpos*q, acc0, phi0);
	pb.calculate_gravity_at(evalpos*q, acc, phi);

#if ERR_DIST
	pb.calculate_gravity_classicalME_at(evalpos*q, phi1);
#endif

	l = get_bn()->get_length();
	cerr << "l: " << l << endl;
	if (l != 1.0) {
	    cerr << "cell size is not unity" << endl;
//	    exit(1);
	}
	PRL(evalpos);
	PR(acc); PRL(phi);
	PR(acc0); PRL(phi0);
	real err, rerr, perr, rperr, perr1, rperr1;
	err = abs(acc-acc0);
	rerr = abs(acc-acc0)/abs(acc0);
	perr = fabs(phi-phi0);
	rperr = fabs((phi-phi0)/phi0);
	PR(err); PRL(rerr);
#if ERR_DIST
	perr1 = fabs(phi1-phi0);
	rperr1 = fabs((phi1-phi0)/phi0);
//	fs << " " << l/abs(evalpos*q) << " " << rerr << " " << err << endl;	
	fs << " " << abs(evalpos*q)
	   << " " << rperr << " " << perr
	   << " " << rperr1 << " " << perr1 << endl;	
#elif ERR_RAD
	fs << " " << fabs(s) << " " << rerr << " " << err << endl;	
//	fs << " " << fabs(s) << " " << rperr << " " << perr << endl;	
#elif ERR_N
	rerrsum += rerr;
	errsum += err;
	rperrsum += rerr;
	perrsum += err;
	}
	rerrsum /= ntry;
	errsum /= ntry;
	rperrsum /= ntry;
	perrsum /= ntry;
	fs << " " << n << " " << rerrsum << " " << errsum << endl;	
//	fs << " " << n << " " << rperrsum << " " << perrsum << endl;	
#else
	incorrect directive
#endif
    }

    if (snapout) {
	fs.close();
    }

    cerr.precision(6);
    exit(0);
}

#endif // TESTCODE
