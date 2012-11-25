#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <limits>
#include <time.h>
/* Using interfaces to cublas and cusparse */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas.h>
/* upcgc interface */
#include "upcgc.h"


extern "C" {

	void cuda_Dvxv(const int *n, double *v1, double *v2, double *v3);

	//int SUPCGJACA(int *ione1, int *ione2, const int *niac, const double *Apc, double *dc, double *zc);

	int SUPCGILU0A(const int *nnzc, const int *niac,
		const int *niapc, const double *Apc, const int *iac,
		const int *jac, const int *iuc, double *dc, double *zc);

	int UPCGC7_INIT(long long *handle_ptr, long long *status_ptr, long long *descr_ptr,
		const int *nnzc, const int *niac, const int *niapc, 
		const int *npc, const int *ndegree,  
		long long *cuda_iac, const int *iac, long long *cuda_jac, const int *jac,  
		long long *cuda_Ac, long long *cuda_Apc, long long *cuda_xc,
		long long *cuda_dc, long long *cuda_zc, long long *cuda_pc, long long *cuda_qc,
		long long *cuda_scl, long long *cuda_scli, long long *cuda_v, long long *cuda_v0, long long *cuda_v1,
		long long *pl_dc, long long *pl_zc, int *ierr)
	{
		void *cuda_jac_gpu, *cuda_iac_gpu;
		void *cuda_Ac_gpu, *cuda_Apc_gpu, *cuda_xc_gpu;
		void *cuda_zc_gpu, *cuda_dc_gpu, *cuda_pc_gpu, *cuda_qc_gpu;
		void *cuda_scl_gpu, *cuda_scli_gpu, *cuda_v_gpu, *cuda_v0_gpu, *cuda_v1_gpu;
		void *pl_dc_gpu, *pl_zc_gpu;

		int cu_status = 0, nnzc1 = *nnzc, niac1 = *niac, niapc1 = *niapc;
		int ndegree1 = *ndegree;

		//add
		int devID;
		int deviceCount;

		//determine the number of GPU devices available
		cudaGetDeviceCount(&deviceCount);
		//make sure that at least one GPU device is available
		if (deviceCount == 0) 
		{
			fprintf(stderr, "UPCGC7_INIT() CUDA error: no devices supporting CUDA.\n");
			exit(-1);
		}
		else 
		{
			devID = 0;
			cudaSetDevice(devID);
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, devID);
			printf("  CUDA device [%d]: %s\n\n", devID, deviceProp.name);
		}

		//create a cusparse object
		cusparseHandle_t cu_handle = 0;
		cu_status = cusparseCreate(&cu_handle);

		if (cu_status != 0)  
		{
			fprintf( stderr, "unable to initalize cusparse\n" );
			return cu_status;
		}

		//create a cusparse matrix object
		cusparseMatDescr_t cu_descr = 0;
		cu_status = cusparseCreateMatDescr(&cu_descr); 

		//throw an error if the cusparse object was not created
		if (cu_status != 0)
		{
			fprintf( stderr, "unable to create cusparse matrix descriptor\n" );
			return cu_status;
		}

		cusparseSetMatType(cu_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(cu_descr,CUSPARSE_INDEX_BASE_ONE);

		//allocate GPU space
		cudaMalloc(&cuda_iac_gpu, (niac1+1) * sizeof(int));
		cudaMalloc(&cuda_jac_gpu, nnzc1 * sizeof(int));
		cudaMalloc(&cuda_Ac_gpu, nnzc1 * sizeof(double)) ;
		cudaMalloc(&cuda_xc_gpu, niac1 * sizeof(double)) ;  
		cudaMalloc(&cuda_dc_gpu, niac1 * sizeof(double)) ;
		cudaMalloc(&cuda_zc_gpu, niac1 * sizeof(double)) ;
		cudaMalloc(&cuda_pc_gpu, niac1 * sizeof(double)) ;
		cudaMalloc(&cuda_qc_gpu, niac1 * sizeof(double)) ;
		cudaMalloc(&cuda_scl_gpu, niac1 * sizeof(double)) ;
		cudaMalloc(&cuda_scli_gpu, niac1 * sizeof(double)) ;

		//pinned memory testing
		//cudaHostAlloc(&pl_dc_gpu,niac1 * sizeof(double),cudaHostAllocMapped);
		//cudaHostAlloc(&pl_zc_gpu,niac1 * sizeof(double),cudaHostAllocMapped);

		//cast GPU pointers to pointers that can be accessed by Fortran 
		*handle_ptr = (long long)cu_handle;
		*descr_ptr = (long long)cu_descr;

		*cuda_iac  = (long long)cuda_iac_gpu;	
		*cuda_jac  = (long long)cuda_jac_gpu;
		*cuda_Ac   = (long long)cuda_Ac_gpu;
		*cuda_xc   = (long long)cuda_xc_gpu;
		*cuda_dc   = (long long)cuda_dc_gpu;
		*cuda_zc   = (long long)cuda_zc_gpu;
		*cuda_pc   = (long long)cuda_pc_gpu;
		*cuda_qc   = (long long)cuda_qc_gpu;
		*cuda_scl  = (long long)cuda_scl_gpu;
		*cuda_scli = (long long)cuda_scli_gpu;

		//preconditioner is used..
		if ( *npc == 1 )
		{
			cudaMalloc(&cuda_Apc_gpu, niac1 * sizeof(double)) ;
		}
		else
		{
			cudaMalloc(&cuda_Apc_gpu, 1 * sizeof(double)) ;
		}
		*cuda_Apc = (long long)cuda_Apc_gpu;
		
		//if ILU0 and MILU0 preconditioner is used..
		if ( (*npc == 2) || (*npc == 3) )
		{
			cudaHostAlloc(&pl_dc_gpu,niac1 * sizeof(double),cudaHostAllocDefault);
			cudaHostAlloc(&pl_zc_gpu,niac1 * sizeof(double),cudaHostAllocDefault);
			//cast GPU pointers to pointers that can be accessed by Fortran 
			*pl_dc = (long long)pl_dc_gpu;
			*pl_zc = (long long)pl_zc_gpu;
		}
		//if polynomial preconditioner is used..
		if ( *npc == 4 )
		{
			cudaMalloc(&cuda_v_gpu, niac1 * sizeof(double)) ;
			cudaMalloc(&cuda_v0_gpu, niac1 * sizeof(double)) ;
			cudaMalloc(&cuda_v1_gpu, niac1 * sizeof(double)) ;

			//cast GPU pointers to pointers that can be accessed by Fortran 
			*cuda_v  = (long long)cuda_v_gpu;
			*cuda_v0  = (long long)cuda_v0_gpu;
			*cuda_v1  = (long long)cuda_v1_gpu;
		}

		//copy iac and jac to card
		cudaMemcpy((void**)*cuda_iac,iac,(*niac+1) * sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_jac,jac,*nnzc * sizeof(int),cudaMemcpyHostToDevice);

		return (0);
	}

	int getprintmaxdouble(const long long *var, const int *length)
	{
		double *temp = NULL;
		double max = -1.0e+32;
		temp = (double*)malloc(sizeof(double) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(double),cudaMemcpyDeviceToHost);

		for (int i=0;i<*length;i++){
			//printf("%d,%3.5e\n",i,temp[i]);
			if (abs(temp[i]) > max) max = abs(temp[i]);
		}

		printf("%3.5e\n\n",max);
		free(temp);
		return (0);
	}

	double upcg_timer(double value)
	{
		double dt = 0.0;
		double  r = 0.0;

		clock_t t = clock();
		if ( value == 0.0 )
		{
			r =  double(t);
		}
		else
		{
			dt = ( double(t) - value ) / CLOCKS_PER_SEC;
			r  = dt;
		}
		return r;
	}

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
		int *ierr)
	{      	 	  
		//locals
		int iiter = 0;
		int idx_deltax = -1, idx_rmax = -1;
		double *zero1 = (double*)malloc(sizeof(double));
		*zero1 = 0.0;
		double *one1 = (double*)malloc(sizeof(double));
		*one1 = 1.0;
		int *izero1 = (int*)malloc(sizeof(int));
		*izero1 = 0;
		int *ione1 = (int*)malloc(sizeof(int));
		*ione1 = 1;
		int *ione2 = (int*)malloc(sizeof(int));
		*ione2 = 1;
		double alpha = 0.0, beta = 0.0;
		double rho0s = double(*rho0), rhos = double(*rho);
		double rmax1 = 0.0, deltax1 = 0.0;

		double pcat1  = 0.0;
		double dpt1   = 0.0;
		double mvt1   = 0.0;
		double axpyt1 = 0.0;
		double vvpt1  = 0.0;
		double misct1 = 0.0;
		double gputt1 = 0.0;

		double pcastart = 0.0;
		double start    = 0.0;

		int k;
		int iicnvg = 0;

		start = upcg_timer(0.0);
		cudaMemcpy((void**)*cuda_Ac,Ac,*nnzc * sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_xc,xc,*niac * sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_dc,dc,*niac * sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_zc,zc,*niac * sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_scl,scl,*niac * sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy((void**)*cuda_scli,scli,*niac * sizeof(double),cudaMemcpyHostToDevice);

		//jacobi, ilu0, or milu0 preconditiones
		//if ( (*npc > 0) && (*npc < 4) )
		//if jacobi preconditioner is used..
		if ( *npc == 1 )
		{
			cudaMemcpy((void**)*cuda_Apc,Apc, *niac * sizeof(double),cudaMemcpyHostToDevice); 
		}
		gputt1 += upcg_timer(start);

		while (iiter < *niter)
		{
			iiter++;

			//preconditioning
			pcastart = 0.0;
			pcastart = upcg_timer(0.0);
			//no preconditioning
			if (*npc == 0)
			{
				start = upcg_timer(0.0);
				cublasDcopy(*niac,(double*)*cuda_dc,1,(double*)*cuda_zc,1);
				misct1 += upcg_timer(start);
			}
			//jacobi
			else if (*npc == 1)
			{
				start = upcg_timer(0.0);
				cuda_Dvxv(niac,(double*) *cuda_Apc,(double*) *cuda_dc,(double*) *cuda_zc);
				vvpt1 += upcg_timer(start);
			}
			//ILU0 and MILU0
			else if ( (*npc == 2) || (*npc == 3) )
			{
				//pageable memory = very slow
				//cudaMemcpy(dc,(void**)(long long)*dc,*niac *sizeof(double), cudaMemcpyDeviceToHost);				
				//SPCGCILU0A(one,nnzc,niac,niapc,Apc,iac,jac,iuc,dc,zc);
				//cudaMemcpy((void**)*zc,zc,*niac * sizeof(double),cudaMemcpyHostToDevice);				

				//page_locked memory
				start = upcg_timer(0.0);
				cudaMemcpy((void**)*pl_dc,(void**)(long long)*cuda_dc,*niac *sizeof(double), cudaMemcpyDeviceToHost);
				gputt1 += upcg_timer(start);

				SUPCGILU0A(nnzc,niac,niapc,Apc,iac,jac,iuc,(double*)*pl_dc,(double*)*pl_zc);

				start = upcg_timer(0.0);
				cudaMemcpy((void**)*cuda_zc,(void**)*pl_zc,*niac * sizeof(double),cudaMemcpyHostToDevice);				
				gputt1 += upcg_timer(start);

				//pinned
				//SPCGCILU0A(one,nnzc,niac,niapc,Apc,iac,jac,iuc,(double*)*pl_dc,(double*)*pl_zc);
			}
			//polynomial preconditioner
			else if ( *npc == 4 )
			{
				double bet = 0.0;

				// zc = M^{-1} * dc
				// assume dc0 = 0
				// zc0 = dc - A * dc0 = dc
				// v1 = zc0 / beta(1) = dc / beta(1)
				start = upcg_timer(0.0);
				cublasDcopy(*niac, (double*)*cuda_dc, 1, (double*)*cuda_v1, 1);
				cublasDscal(*niac, 1.0/pbeta[0], (double*)*cuda_v1, 1);
				misct1 += upcg_timer(start);

				// zc  = zc0 + pgamma(1)*v1;
				start = upcg_timer(0.0);
				cudaMemset((double*)*cuda_zc, 0, *niac*sizeof(double));
				misct1 += upcg_timer(start);
				start = upcg_timer(0.0);
				cublasDaxpy(*niac, pgamma[0], (double*)*cuda_v1, 1, (double*)*cuda_zc, 1);
				axpyt1 += upcg_timer(start);

				// v0 = zeros(n,1);
				start = upcg_timer(0.0);
				cudaMemset((double*)*cuda_v0, 0, *niac*sizeof(double));
				misct1 += upcg_timer(start);

				for (k=0; k<*ndegree; k++) {
					// v = A*v1;
					start = upcg_timer(0.0);
					cusparseDcsrmv((cusparseHandle_t)*handle_ptr,CUSPARSE_OPERATION_NON_TRANSPOSE,
						*niac, *niac, 1.0, (cusparseMatDescr_t)*descr_ptr, (double*)*cuda_Ac, 
						(int*)*cuda_iac,(int*)*cuda_jac, (double*)*cuda_v1, 0.0, (double*)*cuda_v);		
					mvt1 += upcg_timer(start);

					// v  = v - alpha(k)*v1 - bet*v0; 
					start = upcg_timer(0.0);
					cublasDaxpy(*niac, -palpha[k], (double*)*cuda_v1, 1, (double*)*cuda_v, 1);
					axpyt1 += upcg_timer(start);

					if ( k > 0 )
					{
						start = upcg_timer(0.0);
						cublasDaxpy(*niac, -bet, (double*)*cuda_v0, 1, (double*)*cuda_v, 1);
						axpyt1 += upcg_timer(start);
					}

					// v0 = v1;
					start = upcg_timer(0.0);
					cublasDcopy(*niac, (double*)*cuda_v1, 1, (double*)*cuda_v0, 1);
					misct1 += upcg_timer(start);

					// v1 = v/beta(k+1);
					start = upcg_timer(0.0);
					cudaMemset((double*)*cuda_v1, 0, *niac*sizeof(double));
					misct1 += upcg_timer(start);
					start = upcg_timer(0.0);
					cublasDaxpy(*niac, 1.0/pbeta[k+1], (double*)*cuda_v, 1, (double*)*cuda_v1, 1);
					axpyt1 += upcg_timer(start);

					// zc = zc + gamma(k+1)*v1;
					start = upcg_timer(0.0);
					cublasDaxpy(*niac, pgamma[k+1], (double*)*cuda_v1, 1, (double*)*cuda_zc, 1);
					axpyt1 += upcg_timer(start);

					bet = pbeta[k+1];
				}
			}
			pcat1 += upcg_timer(pcastart);

			start = upcg_timer(0.0);
			rhos = cublasDdot(*niac, (double*)*cuda_dc, 1, (double*)*cuda_zc, 1);
			dpt1 += upcg_timer(start);

			if ( iiter > 1 )
			{
				beta = rhos / rho0s;

				start = upcg_timer(0.0);
				cublasDscal(*niac, beta, (double*)*cuda_pc, 1);
				misct1 += upcg_timer(start);

				start = upcg_timer(0.0);
				cublasDaxpy(*niac, 1.0, (double*)*cuda_zc, 1, (double*)*cuda_pc, 1);
				axpyt1 += upcg_timer(start);
			}
			else
			{
				start = upcg_timer(0.0);
				cublasDcopy(*niac, (double*)*cuda_zc, 1, (double*)*cuda_pc, 1);
				misct1 += upcg_timer(start);
			}

			start = upcg_timer(0.0);
			cusparseDcsrmv((cusparseHandle_t)*handle_ptr,CUSPARSE_OPERATION_NON_TRANSPOSE,
				*niac, *niac, 1.0, (cusparseMatDescr_t)*descr_ptr, (double*)*cuda_Ac, 
				(int*)*cuda_iac,(int*)*cuda_jac, (double*)*cuda_pc, 0.0, (double*)*cuda_qc);		
			mvt1 += upcg_timer(start);

			start = upcg_timer(0.0);
			alpha = rhos / cublasDdot(*niac, (double*)*cuda_pc, 1, (double*)*cuda_qc, 1);		
			dpt1 += upcg_timer(start);

			//update x with delta x
			start = upcg_timer(0.0);
			cublasDaxpy(*niac, alpha, (double*)*cuda_pc, 1, (double*)*cuda_xc, 1);
			axpyt1 += upcg_timer(start);

			//update residual
			start = upcg_timer(0.0);
			cublasDaxpy(*niac, -alpha, (double*)*cuda_qc, 1, (double*)*cuda_dc, 1);
			axpyt1 += upcg_timer(start);

			//calculate maximum head change
			//calculate deltax
			start = upcg_timer(0.0);
			cudaMemset((double*)*cuda_qc, 0, *niac*sizeof(double));
			misct1 += upcg_timer(start);
			start = upcg_timer(0.0);
			cublasDaxpy(*niac, alpha, (double*)*cuda_pc, 1, (double*)*cuda_qc, 1);
			axpyt1 += upcg_timer(start);
			start = upcg_timer(0.0);
			cuda_Dvxv(niac,(double*) *cuda_qc,(double*) *cuda_scl,(double*) *cuda_qc);
			vvpt1 += upcg_timer(start);
			start = upcg_timer(0.0);
			idx_deltax = cublasIdamax(*niac, (double*)*cuda_qc, 1);
			misct1 += upcg_timer(start);
			//find maximum unscaled head change
			start = upcg_timer(0.0);
			cudaMemcpy(deltax,(void**)(long long)*cuda_qc+idx_deltax-1,sizeof(double), cudaMemcpyDeviceToHost);
			gputt1 += upcg_timer(start);

			//calculate residual error
			//unscale residual
			start = upcg_timer(0.0);
			cuda_Dvxv(niac,(double*) *cuda_dc,(double*) *cuda_scli,(double*) *cuda_qc);
			vvpt1 += upcg_timer(start);
			start = upcg_timer(0.0);
			idx_rmax   = cublasIdamax(*niac, (double*)*cuda_qc, 1);
			misct1 += upcg_timer(start);
			//find maximum unscaled residual
			start = upcg_timer(0.0);
			cudaMemcpy(rmax,(void**)(long long)*cuda_qc+idx_rmax-1,sizeof(double), cudaMemcpyDeviceToHost);
			gputt1 += upcg_timer(start);

			*deltax = abs(*deltax);
			*rmax = abs(*rmax);

			//convergence check
			if ((*deltax <= *hclose) && (*rmax <= *rclose)) iicnvg = 1;
			if (*mxiter == 1)
			{
				if ( iicnvg == 1 ) *icnvg = 1;
			}
			else 
			{
				if (iiter == 1 && iicnvg == 1) *icnvg = 1;
			}
			if ( iicnvg == 1 ) break;

			rho0s = rhos;

		}

		*rho0 = rho0s;
		*rho  =  rhos;
		*iter = iiter;

		//update timers
		*pcat  += pcat1;
		*dpt   += dpt1;
		*mvt   += mvt1;
		*axpyt += axpyt1;
		*vvpt  += vvpt1;
		*misct += misct1;
		*gputt += gputt1;

		//copy cuda_xc back to xc
		start = upcg_timer(0.0);
		cudaMemcpy(xc,(void**)(long long)*cuda_xc,*niac *sizeof(double), cudaMemcpyDeviceToHost);
		gputt1 += upcg_timer(start);

		return (0);
	}

	//free up all of the GPU resources
	int UPCGC7_FINAL(const int *npc, long long *handle_ptr, 
		long long *cuda_jac, long long *cuda_iac, 
		long long *cuda_Ac, long long *cuda_Apc, long long *cuda_xc, 
		long long *cuda_dc, long long *cuda_zc, long long *cuda_pc, long long *cuda_qc,
		long long *cuda_scl, long long *cuda_scli, long long *cuda_v, long long *cuda_v0, long long *cuda_v1, 
		long long *pl_dc, long long *pl_zc, 
		int *ierr)
	{
		cusparseDestroy((cusparseHandle_t)*handle_ptr);
		cudaFree((void**)*cuda_jac);
		cudaFree((void**)*cuda_iac);
		cudaFree((void**)*cuda_Ac);
		cudaFree((void**)*cuda_Apc);
		cudaFree((void**)*cuda_xc);
		cudaFree((void**)*cuda_zc);
		cudaFree((void**)*cuda_dc);
		cudaFree((void**)*cuda_pc);
		cudaFree((void**)*cuda_qc);
		cudaFree((void**)*cuda_scl);
		cudaFree((void**)*cuda_scli);
		if ( *npc == 4 )
		{
			cudaFree((void**)*cuda_v);
			cudaFree((void**)*cuda_v0);
			cudaFree((void**)*cuda_v1);
		}
		if ( (*npc == 2) || (*npc == 3) )
		{
			cudaFreeHost((void**)*pl_dc);
			cudaFreeHost((void**)*pl_zc);
		}
		cudaThreadExit();
		return (0);
	}

	//****************************************
	//The Bone Yard...enter at your own risk
	int compprintfloat(const int *var, const float *var2, const int *length)
	{
		float *temp = NULL;
		temp = (float*)malloc(sizeof(float) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(float),cudaMemcpyDeviceToHost);
		for (int i=0;i<*length;i++){
			printf("%d,%3.5f,%3.5f\n",i,temp[i],var2[i]);
		}
		free(temp);
		return (0);
	}
	int getprintdouble(const long long *var, const int *length)
	{
		double *temp = NULL;
		temp = (double*)malloc(sizeof(double) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(double),cudaMemcpyDeviceToHost);
		for (int i=0;i<*length;i++){
			printf("%d,%3.5e\n",i,temp[i]);
		}
		free(temp);
		return (0);
	}
	int getprintmaxfloat(const int *var, const int *length)
	{
		float *temp = NULL;
		float max   = -FLT_MAX; //-1.0e+32;
		temp = (float*)malloc(sizeof(float) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(float),cudaMemcpyDeviceToHost);

		for (int i=0;i<*length;i++){
			//printf("%d,%3.5e\n",i,temp[i]);
			if (abs(temp[i]) > max) max = abs(temp[i]);
		}

		printf("%3.5e\n",max);
		free(temp);
		return (0);
	}

	int compprintint(const int *var, const int *var2, const int *length)
	{
		int *temp = NULL;
		temp = (int*)malloc(sizeof(int) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(int),cudaMemcpyDeviceToHost);
		for (int i=0;i<*length;i++){
			printf("%d,%d,%d\n",i,temp[i],var2[i]);
		}
		free(temp);
		return (0);
	}

	int getprintint(const int *var, const int *length)
	{
		int *temp = NULL;
		temp = (int*)malloc(sizeof(int) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(int),cudaMemcpyDeviceToHost);
		for (int i=0;i<*length;i++){
			printf("%d,%d\n",i,temp[i]);
		}
		free(temp);
		return (0);
	}

	//casting single to double and back
	int d2s(const double *doub, float *sing, const int *length)
	{
		for (int i=0;i<*length;i++){ 
			sing[i] = (float)doub[i];
			//printf("%f,%f\n",sing[i],doub[i]);
		}
		return (0);
	}

	int s2d(double *doub,const float *sing, const int *length)
	{
		for (int i=0;i<*length;i++){
			doub[i] = (double)sing[i];
		}
		return (0);
	}

	//extract diagonal of A for Jacobi
	void diag_d(const int *niac,const int *iac, const double *A, double *a_diag)
	{
		for (int n=0;n<*niac;n++)
		{
			a_diag[n] = A[iac[n]-1];
			//printf("%d,%d\n",iac[n],A[iac[n]-1]);
		}
		return;
	}

	double getmaxdouble(const long long *var, const int *length)
	{
		double *temp = NULL;
		double max = -1.0e+32;
		temp = (double*)malloc(sizeof(double) * (*length));
		cudaMemcpy(temp,(void**)*var,*length * sizeof(double),cudaMemcpyDeviceToHost);

		for (int i=0;i<*length;i++){
			//printf("%d,%3.5e\n",i,temp[i]);
			if (abs(temp[i]) > max) max = abs(temp[i]);
		}

		free(temp);
		return (max);
	}

	void write_time(char* location, clock_t start, FILE* f)
	{
		clock_t end = clock();
		float total = ((float)end - (float)start)/CLOCKS_PER_SEC;
		//FILE *f_out = fopen(fName,"w");
		//printf("%s - %f seconds\n",location,total);
		fprintf(f,"%s - %f seconds\n",location,total);
		//fclose(f_out);
		return;
	}

}