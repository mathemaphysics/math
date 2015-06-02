/**
 * Variables needed for configuration
 */

/* Basic string constants */
char buf[TOKEN_BUFFER_LENGTH];
char bug[TOKEN_BUFFER_LENGTH];
char *tok[TOKEN_BUFFER_LENGTH];
char *tol[TOKEN_BUFFER_LENGTH];

/* Configuration switches */
int b_points_loaded = 0;
int b_config_loaded = 0;
int b_quad_loaded = 0;
int b_rdn_set = 0;
int b_domain_set = 0;
int b_neig_set = 0;
int b_load_stiff_mat = 0;
int b_load_overl_mat = 0;
int b_save_stiff_mat = 0;
int b_save_overl_mat = 0;
int b_save_proj_mat = 0;
int b_have_stiff_mat = 0;
int b_have_overl_mat = 0;
int b_boundary_proj = 0;
int b_treig_routine_set = 0;
int b_grid_set_dx = 0;
int b_grid_set_x0 = 0;
int b_grid_set_wx = 0;
int b_out_eig_range_set = 0;
int b_use_singular = 1;
int b_nuc_loaded = 0;
int b_use_external_potential = 1;

/* String variables for configuration */
char s_overl_fname[TOKEN_BUFFER_LENGTH] = "amat.out";
char s_stiff_fname[TOKEN_BUFFER_LENGTH] = "bmat.out";
char s_proj_fname[TOKEN_BUFFER_LENGTH] = "wmat.out";
char s_treig_routine[TOKEN_BUFFER_LENGTH] = "treiglr";

/* Integer variables which need to be set */
int i_sing_order = 1;
int i_max_lanczos_steps = 100;

/* Some floating point variables needed */
double f_bndry_tol = 0.05;

