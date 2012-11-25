extern "C"
{
	int UPCGC7_INIT(long long *handle_ptr, long long *status_ptr, long long *descr_ptr,
		const int *nnzc, const int *niac, const int *niapc, 
		const int *npc, const int *ndegree,
		long long *cuda_iac, const int *iac, long long *cuda_jac, const int *jac,  
		long long *cuda_Ac, long long *cuda_Apc, long long *cuda_xc,
		long long *cuda_dc, long long *cuda_zc, long long *cuda_pc, long long *cuda_qc,
		long long *cuda_scl, long long *cuda_scli, long long *cuda_v, long long *cuda_v0, long long *cuda_v1,
		long long *pl_dc, long long *pl_zc, int *ierr);

	int UPCGC7(int *iter, long long *handle_ptr, long long *status_ptr,	long long *descr_ptr, 
		const int *nnzc, const int *niac, const int *niapc, 
		const int *npc, const int *ndegree,
		long long *cuda_iac, const int *iac, long long *cuda_jac, const int *jac, const int *iuc,
		long long *cuda_Ac, const double *Ac, long long *cuda_Apc, const double *Apc, 
		long long *cuda_xc, double *xc,
		long long *cuda_dc, double *dc, long long *cuda_zc, double *zc,
		long long *cuda_pc, long long *cuda_qc,
		double *palpha, double *pbeta, double *pgamma, 
		long long *cuda_scl, double *scl, long long *cuda_scli, double *scli, long long *cuda_v, double *v, long long *cuda_v0, double *v0, long long *cuda_v1, double *v1, 
		long long *pl_dc, long long *pl_zc,
		double *rho0, double *rho, 
		const int *mxiter,int *icnvg, int *niter, const double *hclose, const double *rclose,
		double *deltax, double *rmax, 
		double *pcat, double *dpt, double *mvt, double *axpyt, double *vvpt, double *misct, double *gputt,
		int *ierr);

	int UPCGC7_FINAL(const int *npc, long long *handle_ptr, 
    long long *cuda_jac, long long *cuda_iac, 
		long long *cuda_Ac, long long *cuda_Apc, long long *cuda_xc,
		long long *cuda_dc, long long *cuda_zc, long long *cuda_pc,	long long *cuda_qc,
		long long *cuda_scl, long long *cuda_scli, long long *cuda_v, long long *cuda_v0, long long *cuda_v1, 
		long long *pl_dc, long long *pl_zc, 
		int *ierr);

}