/*
 * prototypes for newton0 , derived using
 *     cproto -e -I$NEMOINC -I$NEMOLIB  -D.....<flag>.....
 * where <flags> are one of:
 * 	-DTREE                 
 *    	-DTREE -DQUADPOLE      
 *    	-DREGULARIZATION
 *      -DEXTRAPOLATION 
 *
 * 5-feb-2001    Peter Teuben   -- in an attempt to resurrect "Newton0"
 */

/* binaryin.c */
extern systptr introduce_system(string infile, string *headlineptr, string *announceptr);

/* binaryout.c */
extern void open_binary_output(string outfile, string headline);
extern void close_binary_output(void);
extern void put_full_snapshot(stateptr the_state);
extern void put_diagnostics_snapshot(stateptr the_state);
extern void put_diagnostics(diagptr diags);

/* bodyalgebra.c */
extern bodyptr mk_bodies(int npart);
extern void clear_bodies(bodyptr bod, int npart);
extern bodyptr cp_bodies(bodyptr old_bod, int npart);
extern bodyptr mul_bodies(bodyptr old_bod, int npart, real scale_factor);
extern bodyptr add_bodies(bodyptr bod1, bodyptr bod2, int npart);
extern bodyptr sub_bodies(bodyptr bod1, bodyptr bod2, int npart);
extern void scale_bodies(bodyptr bod, int npart, real scale_factor);
extern void inc_bodies(bodyptr bod, bodyptr d_bod, int npart);
extern void dec_bodies(bodyptr bod, bodyptr d_bod, int npart);
extern ebodyptr mk_ebodies(int npart);
extern void clear_ebodies(ebodyptr ebod, int npart);
extern ebodyptr cp_ebodies(ebodyptr old_ebod, int npart);
extern ebodyptr mul_ebodies(ebodyptr old_ebod, int npart, real scale_factor);
extern ebodyptr add_ebodies(ebodyptr ebod1, ebodyptr ebod2, int npart);
extern ebodyptr sub_ebodies(ebodyptr ebod1, ebodyptr ebod2, int npart);
extern void scale_ebodies(ebodyptr ebod, int npart, real scale_factor);
extern void inc_ebodies(ebodyptr ebod, ebodyptr d_ebod, int npart);
extern void dec_ebodies(ebodyptr ebod, ebodyptr d_ebod, int npart);
extern mbodyptr mk_mbodies(int npart);
extern void clear_mbodies(mbodyptr mbod, int npart);
extern mbodyptr cp_mbodies(mbodyptr old_mbod, int npart);
extern mbodyptr mul_mbodies(mbodyptr old_mbod, int npart, real scale_factor);
extern mbodyptr add_mbodies(mbodyptr mbod1, mbodyptr mbod2, int npart);
extern mbodyptr sub_mbodies(mbodyptr mbod1, mbodyptr mbod2, int npart);
extern void scale_mbodies(mbodyptr mbod, int npart, real scale_factor);
extern void inc_mbodies(mbodyptr mbod, mbodyptr d_mbod, int npart);
extern void dec_mbodies(mbodyptr mbod, mbodyptr d_mbod, int npart);
extern pbodyptr mk_pbodies(int npart);
extern void clear_pbodies(pbodyptr pbod, int npart);
extern pbodyptr cp_pbodies(pbodyptr old_pbod, int npart);
extern pbodyptr mul_pbodies(pbodyptr old_pbod, int npart, real scale_factor);
extern pbodyptr add_pbodies(pbodyptr pbod1, pbodyptr pbod2, int npart);
extern pbodyptr sub_pbodies(pbodyptr pbod1, pbodyptr pbod2, int npart);
extern void scale_pbodies(pbodyptr pbod, int npart, real scale_factor);
extern void inc_pbodies(pbodyptr pbod, pbodyptr d_pbod, int npart);
extern void dec_pbodies(pbodyptr pbod, pbodyptr d_pbod, int npart);
extern cbodyptr mk_cbodies(int npart);
extern void clear_cbodies(cbodyptr cbod, int npart);
extern cbodyptr cp_cbodies(cbodyptr old_cbod, int npart);
extern cbodyptr mul_cbodies(cbodyptr old_cbod, int npart, real scale_factor);
extern cbodyptr add_cbodies(cbodyptr cbod1, cbodyptr cbod2, int npart);
extern cbodyptr sub_cbodies(cbodyptr cbod1, cbodyptr cbod2, int npart);
extern void scale_cbodies(cbodyptr cbod, int npart, real scale_factor);
extern void inc_cbodies(cbodyptr cbod, cbodyptr d_cbod, int npart);
extern void dec_cbodies(cbodyptr cbod, cbodyptr d_cbod, int npart);

/* bodyconversion.c */
extern cbodyptr sel_b_pos(bodyptr sys, int npart);
extern cbodyptr sel_b_vel(bodyptr sys, int npart);
extern cbodyptr sel_b_acc(bodyptr sys, int npart);
extern pbodyptr sel_b_phase(bodyptr sys, int npart);
extern pbodyptr sel_b_phasedot(bodyptr sys, int npart);
extern mbodyptr sel_b_msys(bodyptr sys, int npart);
extern ebodyptr sel_b_esys(bodyptr sys, int npart);
extern cbodyptr sel_b_epos(ebodyptr esys, int npart);
extern cbodyptr sel_b_evel(ebodyptr esys, int npart);
extern pbodyptr sel_b_ephase(ebodyptr esys, int npart);
extern mbodyptr sel_b_emsys(ebodyptr esys, int npart);
extern cbodyptr sel_b_mpos(mbodyptr msys, int npart);
extern cbodyptr sel_b_mvel(mbodyptr msys, int npart);
extern pbodyptr sel_b_mphase(mbodyptr msys, int npart);
extern cbodyptr sel_b_ppos(pbodyptr psys, int npart);
extern cbodyptr sel_b_pvel(pbodyptr psys, int npart);
extern pbodyptr cons_b_pos_vel(cbodyptr csys_pos, cbodyptr csys_vel, int npart);
extern bodyptr cons_b_e_acc(ebodyptr esys, cbodyptr csys, int npart);
extern void renov_b_pos(bodyptr sys, cbodyptr csys, int npart);
extern void renov_b_vel(bodyptr sys, cbodyptr csys, int npart);
extern void renov_b_acc(bodyptr sys, cbodyptr csys, int npart);
extern void renov_b_phase(bodyptr sys, pbodyptr psys, int npart);
extern void renov_b_msys(bodyptr sys, mbodyptr msys, int npart);
extern void renov_b_esys(bodyptr sys, ebodyptr esys, int npart);
extern void renov_b_epos(ebodyptr esys, cbodyptr csys, int npart);
extern void renov_b_evel(ebodyptr esys, cbodyptr csys, int npart);
extern void renov_b_ephase(ebodyptr esys, pbodyptr psys, int npart);
extern void renov_b_emsys(ebodyptr esys, mbodyptr msys, int npart);
extern void renov_b_mpos(mbodyptr msys, cbodyptr csys, int npart);
extern void renov_b_mvel(mbodyptr msys, cbodyptr csys, int npart);
extern void renov_b_mphase(mbodyptr msys, pbodyptr psys, int npart);
extern void renov_b_ppos(pbodyptr psys, cbodyptr csys, int npart);
extern void renov_b_pvel(pbodyptr psys, cbodyptr csys, int npart);
extern void annex_b_pos(bodyptr sys, cbodyptr csys, int npart);
extern void annex_b_vel(bodyptr sys, cbodyptr csys, int npart);
extern void annex_b_acc(bodyptr sys, cbodyptr csys, int npart);
extern void annex_b_phase(bodyptr sys, pbodyptr psys, int npart);
extern void annex_b_msys(bodyptr sys, mbodyptr msys, int npart);
extern void annex_b_esys(bodyptr sys, ebodyptr esys, int npart);
extern void annex_b_epos(ebodyptr esys, cbodyptr csys, int npart);
extern void annex_b_evel(ebodyptr esys, cbodyptr csys, int npart);
extern void annex_b_ephase(ebodyptr esys, pbodyptr psys, int npart);
extern void annex_b_emsys(ebodyptr esys, mbodyptr msys, int npart);
extern void annex_b_mpos(mbodyptr msys, cbodyptr csys, int npart);
extern void annex_b_mvel(mbodyptr msys, cbodyptr csys, int npart);
extern void annex_b_mphase(mbodyptr msys, pbodyptr psys, int npart);
extern void annex_b_ppos(pbodyptr psys, cbodyptr csys, int npart);
extern void annex_b_pvel(pbodyptr psys, cbodyptr csys, int npart);

/* create.c */
extern bodyptr create_plum(int npart);

/* diagnose.c */
extern void diagnostics(stateptr the_state, int init_flag);
extern void final_diagnostics(stateptr the_state);
extern void compute_energy(systptr sys, diagptr diags);
extern real compute_total_energy_only(systptr sys, diagptr diags);
extern void grinding_halt(stateptr the_state);
extern bodyptr cp_NMAX_bodies(bodyptr old_bod, int npart, int n_alloc);

/* differentiate.c */
extern void init_deriv_system(systptr sys, specptr specs);
extern void deriv_system(systptr sys, specptr specs);
extern psystptr d_system(systptr sys, specptr specs, real ds);
extern void push_system(systptr sys, psystptr d_sys);

/* differentiatetree.c -DTREE */
#ifdef TREE
 extern void tree_tree_interaction(systptr sys, specptr specs);
 extern void body_tree_interaction(bodyptr bp, nodeptr root, specptr specs);
#endif

/* differentiatereg.c -DREGULARIZATION */
#ifdef REGULARIZATION 
extern void heggie_mikkola_equations_of_motion(systptr reg_sys, specptr specs);
#endif

/* differentiateext.c -DEXTRAPOLATION */
#ifdef EXTRAPOLATION 
extern void integrate_ext_invidual_timestep(stateptr the_state, real final_time);
#endif

/* integrate.c */
extern void system_step(systptr sys, specptr specs, real ds);
extern void init_system_rev(systptr sys, specptr specs);
extern void system_rev(systptr sys, specptr specs);

/* newton0.c */

/* orbit.c */
extern void rev_engine(stateptr new_state);
extern void evolve(stateptr some_state);

/* out.c */
extern void write_initial_output(stateptr new_state);
extern void write_final_output(stateptr old_state);
extern void maj_out(stateptr the_state, int init_flag);
extern void min_out(stateptr the_state, int init_flag);
extern void init_ascii_out(stateptr new_state);
extern void ascii_out(stateptr the_state);
extern void talkin_talkin(void);
extern void announce(string what_happened, real when_did_it_happen);


/* save.c */
extern void init_save_state(stateptr the_state);
extern void save_state(stateptr the_state);
extern stateptr restore_state(string old_file);

/* soften.c */
extern proc softener(string type_of_softening);

/* statealgebra.c */
extern stateptr mk_state(void);
extern stateptr rm_state(stateptr old_state);
extern specptr mk_specs(void);
extern void rm_specs(specptr old_specs);
extern ctrlptr mk_ctrls(void);
extern void rm_ctrls(ctrlptr old_ctrls);
extern diagptr mk_diags(void);
extern void rm_diags(diagptr old_diags);

/* systemalgebra.c */
extern systptr mk_empty_system(void);
extern systptr mk_system(int npart);
extern void rm_system(systptr old_sys);
extern void clear_system(systptr sys);
extern systptr cp_system(systptr old_sys);
extern systptr mul_system(systptr old_sys, real scale_factor);
extern systptr add_system(systptr sys1, systptr sys2);
extern systptr sub_system(systptr sys1, systptr sys2);
extern void scale_system(systptr sys, real scale_factor);
extern void inc_system(systptr sys, systptr d_sys);
extern void dec_system(systptr sys, systptr d_sys);
extern esystptr mk_empty_esystem(void);
extern esystptr mk_esystem(int npart);
extern void rm_esystem(esystptr old_esys);
extern void clear_esystem(esystptr esys);
extern esystptr cp_esystem(esystptr old_esys);
extern esystptr mul_esystem(esystptr old_esys, real scale_factor);
extern esystptr add_esystem(esystptr esys1, esystptr esys2);
extern esystptr sub_esystem(esystptr esys1, esystptr esys2);
extern void scale_esystem(esystptr esys, real scale_factor);
extern void inc_esystem(esystptr esys, esystptr d_esys);
extern void dec_esystem(esystptr esys, esystptr d_esys);
extern msystptr mk_empty_msystem(void);
extern msystptr mk_msystem(int npart);
extern void rm_msystem(msystptr old_msys);
extern void clear_msystem(msystptr msys);
extern msystptr cp_msystem(msystptr old_msys);
extern msystptr mul_msystem(msystptr old_msys, real scale_factor);
extern msystptr add_msystem(msystptr msys1, msystptr msys2);
extern msystptr sub_msystem(msystptr msys1, msystptr msys2);
extern void scale_msystem(msystptr msys, real scale_factor);
extern void inc_msystem(msystptr msys, msystptr d_msys);
extern void dec_msystem(msystptr msys, msystptr d_msys);
extern psystptr mk_empty_psystem(void);
extern psystptr mk_psystem(int npart);
extern void rm_psystem(psystptr old_psys);
extern void clear_psystem(psystptr psys);
extern psystptr cp_psystem(psystptr old_psys);
extern psystptr mul_psystem(psystptr old_psys, real scale_factor);
extern psystptr add_psystem(psystptr psys1, psystptr psys2);
extern psystptr sub_psystem(psystptr psys1, psystptr psys2);
extern void scale_psystem(psystptr psys, real scale_factor);
extern void inc_psystem(psystptr psys, psystptr d_psys);
extern void dec_psystem(psystptr psys, psystptr d_psys);
extern csystptr mk_empty_csystem(void);
extern csystptr mk_csystem(int npart);
extern void rm_csystem(csystptr old_csys);
extern void clear_csystem(csystptr csys);
extern csystptr cp_csystem(csystptr old_csys);
extern csystptr mul_csystem(csystptr old_csys, real scale_factor);
extern csystptr add_csystem(csystptr csys1, csystptr csys2);
extern csystptr sub_csystem(csystptr csys1, csystptr csys2);
extern void scale_csystem(csystptr csys, real scale_factor);
extern void inc_csystem(csystptr csys, csystptr d_csys);
extern void dec_csystem(csystptr csys, csystptr d_csys);

/* systemconversion.c */
extern csystptr sel_pos(systptr sys);
extern csystptr sel_vel(systptr sys);
extern csystptr sel_acc(systptr sys);
extern psystptr sel_phase(systptr sys);
extern psystptr sel_phasedot(systptr sys);
extern msystptr sel_msys(systptr sys);
extern esystptr sel_esys(systptr sys);
extern csystptr sel_epos(esystptr esys);
extern csystptr sel_evel(esystptr esys);
extern psystptr sel_ephase(esystptr esys);
extern msystptr sel_emsys(esystptr esys);
extern csystptr sel_mpos(msystptr msys);
extern csystptr sel_mvel(msystptr msys);
extern psystptr sel_mphase(msystptr msys);
extern csystptr sel_ppos(psystptr psys);
extern csystptr sel_pvel(psystptr psys);
extern psystptr cons_pos_vel(csystptr csys_pos, csystptr csys_vel);
extern systptr cons_e_acc(esystptr esys, csystptr csys);
extern void renov_pos(systptr sys, csystptr csys);
extern void renov_vel(systptr sys, csystptr csys);
extern void renov_acc(systptr sys, csystptr csys);
extern void renov_phase(systptr sys, psystptr psys);
extern void renov_msys(systptr sys, msystptr msys);
extern void renov_esys(systptr sys, esystptr esys);
extern void renov_epos(esystptr esys, csystptr csys);
extern void renov_evel(esystptr esys, csystptr csys);
extern void renov_ephase(esystptr esys, psystptr psys);
extern void renov_emsys(esystptr esys, msystptr msys);
extern void renov_mpos(msystptr msys, csystptr csys);
extern void renov_mvel(msystptr msys, csystptr csys);
extern void renov_mphase(msystptr msys, psystptr psys);
extern void renov_ppos(psystptr psys, csystptr csys);
extern void renov_pvel(psystptr psys, csystptr csys);
extern void annex_pos(systptr sys, csystptr csys);
extern void annex_vel(systptr sys, csystptr csys);
extern void annex_acc(systptr sys, csystptr csys);
extern void annex_phase(systptr sys, psystptr psys);
extern void annex_msys(systptr sys, msystptr msys);
extern void annex_esys(systptr sys, esystptr esys);
extern void annex_epos(esystptr esys, csystptr csys);
extern void annex_evel(esystptr esys, csystptr csys);
extern void annex_ephase(esystptr esys, psystptr psys);
extern void annex_emsys(esystptr esys, msystptr msys);
extern void annex_mpos(msystptr msys, csystptr csys);
extern void annex_mvel(msystptr msys, csystptr csys);
extern void annex_mphase(msystptr msys, psystptr psys);
extern void annex_ppos(psystptr psys, csystptr csys);
extern void annex_pvel(psystptr psys, csystptr csys);

/* timestep.c */
extern rproc find_timestep_method(string type_of_timestep);
extern real constant_timestep(stateptr some_state);
extern real collision_time_timestep(stateptr some_state);

/* transformext.c */
#ifdef EXTRAPOLATION
extern void setup_ext(systptr initial_sys, specptr initial_specs);
extern void read_ext(systptr sys, specptr specs, real readout_time);
#endif

/* transformreg.c */
#ifdef REGULARIZATION
extern systptr setup_reg(systptr nonreg_sys);
extern systptr mk_reg_version(systptr nonreg_sys);
extern void nonreg_to_reg(systptr nonreg_sys, systptr reg_sys);
extern void reg_to_nonreg(systptr reg_sys, systptr nonreg_sys);
extern void init_ks_matrix(real ks_matrix[4 ][4 ], real vec[4 ]);
extern realptr mk_massmatrix(int npart);
#endif

/* transformtree.c -DTREE */
#ifdef TREE
extern void mk_tree(systptr sys);
#endif

