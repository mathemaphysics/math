
/**
 * Kernel function to generate the overlap and Hamiltonian matrices
 */
__global__ void puksham_mbuild_cuda( punity_t *pu, int mdim, int nbp, shape_t *dm, int *nlbase, int **lbase, int *ltg, nuclei_t *nuc, int quadn, double *qpts, double *qwts,
					 long **jap, double **Ap, long **jbp, double **Bp, long cc, long c_spmat_inc,
						int b_use_external_potential, int b_load_overl_mat, int b_load_stiff_mat,
							int b_use_singular, int i_sing_order, int *b_have_stiff_mat, int *b_have_overl_mat )
{
	
}

int main()
{

}

