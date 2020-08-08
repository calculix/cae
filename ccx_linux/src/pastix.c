/* CalculiX - A 3-dimensional finite element program */
/*Copyright (C) 1998-2020 Guido Dhondt*/

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License as*/
/* published by the Free Software Foundation(version 2);*/
/**/

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of*/ 
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the*/
/* GNU General Public License for more details.*/

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifdef PASTIX

#include <pastix.h>
#include <spm.h>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "pastix.h"

// Time structs for benchmarking
extern struct timespec totalCalculixTimeStart, totalCalculixTimeEnd;
double totalPastixTime;

// Variables for storing reuse history
int totalIterations = 0;
int totalReused = 0;

// Current sparse matrix in STI format
double* auPtr = NULL;
double* adPtr = NULL;
ITG *icolTotal = NULL, *irowTotal = NULL;

// Matrix data from previous iteration
ITG neqPrev=0, nzsPrev=0;
ITG *icolPrev=NULL,*irowPrev=NULL,*jqPrev=NULL;
ITG inputformatPrev=-1;
ITG basePrev=1;

// Current sparse matrix in CSC
double *aupastix=NULL;
ITG *icolpastix=NULL,*irowpastix=NULL;

// Global variable that indicates whether we are currently reusing or not
char redo = 1;

// Global variable which data set was previously used (basic/radiation) and now
#define BASIC 1
#define AS 2
char modePrev = BASIC;
char mode = BASIC;
char needCleanup=0;

// PaStiX configuration
spm_int_t iparm_basic[IPARM_SIZE];
spm_int_t iparm_as[IPARM_SIZE];
double dparm_basic[DPARM_SIZE];
double dparm_as[DPARM_SIZE];
spm_int_t *iparm = iparm_basic;
double *dparm = dparm_basic;
pastix_data_t* pastix_data = NULL;
spmatrix_t *spm = NULL;

// GPU active or not
char gpu = 0;


// Store how many nzs the merged Matrix has
ITG nzsTotal = 0;

// Size of allocated space for sparse matrix
ITG pastix_nnzBound = 0;

// Number of iterations that failed with mixed precision
ITG mixedFailed = 0;

// indicates whether this is the first invocation of PaStiX or not
char firstIter = 1;

// When this flag is activated, PaStiX will not reuse in the next iteration
char forceRedo = 0;

// Use double or mixed precision
char mixed = 1;
char globDoublePrecision = 0;

// This is set to one, when to many iterations with mixed precision did not converge
char stickToDouble = 0;


// Pointers for faster matrix transpose
ITG *irowacc = NULL;
ITG *irowPrediction = NULL;

// Number of threads
ITG nthread_mkl=0;

struct pastix_data_s {
    int totalIterations;
    int totalReused;
    ITG *icolTotal;
    ITG *irowTotal;
    ITG neqPrev;
    ITG nzsPrev;
    ITG *icolPrev;
    ITG *irowPrev;
    ITG *jqPrev;
    ITG inputformatPrev;
    ITG basePrev;
    double *aupastix;
    ITG *icolpastix;
    ITG *irowpastix;
    char redo;
    ITG nzsTotal;
    ITG pastix_nnzBound;
    ITG mixedFailed;
    char firstIter;
    char forceRedo;
    char globDoublePrecision;
    char stickToDouble;
    ITG *irowacc;
    ITG *irowPrediction;
    spm_int_t *iparm;
    double *dparm;
    char gpu;
    char mixed;
    pastix_data_t* pastix_data;
    spmatrix_t *spm;
};

typedef struct pastix_data_s pastix_data_object;

pastix_data_object pastix_mode_basic = {
    0,0,NULL,NULL,0,0,NULL,NULL,NULL,-1,0,
    NULL,NULL,NULL,1,0,0,0,1,0,0,0,NULL,NULL,iparm_basic,dparm_basic,
    0,1,NULL,NULL
    };

pastix_data_object pastix_mode_as = {
    0,0,NULL,NULL,0,0,NULL,NULL,NULL,-1,0,
    NULL,NULL,NULL,1,0,0,0,1,0,0,0,NULL,NULL,iparm_as,dparm_as,
    0,1,NULL,NULL
    };


// Initializes and configurates PaStiX environment. Also forwards the sparse matrix pointers
void pastix_init(double *ad, double *au, double *adb, double *aub, 
	        double *sigma,ITG *icol, ITG *irow, 
			ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
			ITG *jq, ITG *nzs3){
	// if reusing, only update the value pointer of the sparse matrix
	if(!redo){
	    pastixResetSteps(pastix_data);
		spm->values = aupastix;

		printf("\n");
		spmPrintInfo( spm, stdout );
		printf("\n");
		
		return;
	}

	ITG nthread, nthread_v;
	char *env;

	 /*set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
	 must be done once*/
	if (nthread_mkl == 0) {
	nthread=1;
	env=getenv("MKL_NUM_THREADS");
	if(env) {
	nthread=atoi(env);}
	else {
	env=getenv("OMP_NUM_THREADS");
	if(env) {nthread=atoi(env);}
	}
	env=getenv("CCX_NPROC_EQUATION_SOLVER");
	if(env) {
	nthread_v=atoi(env);
	if (nthread_v <= nthread) {nthread=nthread_v;}
	}
	if (nthread < 1) {nthread=1;}
	nthread_mkl=nthread;
	}
	
	// Init integer and double parameters with default values
	pastixInitParam( iparm, dparm );
	
	// Set best PaStiX parameters for CalculiX usage
    iparm[IPARM_ORDERING]  				= PastixOrderScotch;
    if( mode == AS ){
	    iparm[IPARM_SCHEDULER] 			= PastixSchedStatic;
    }
    else{
	    iparm[IPARM_SCHEDULER] 			= PastixSchedParsec;
    }
	iparm[IPARM_THREAD_NBR]				= nthread_mkl;
	iparm[IPARM_GPU_NBR]   				= (int) gpu;
	iparm[IPARM_FLOAT]			 		= globDoublePrecision ? 3 : 2;
	iparm[IPARM_MIN_BLOCKSIZE] 			= 1024;
	iparm[IPARM_MAX_BLOCKSIZE] 			= 2048;
	iparm[IPARM_FACTORIZATION] 			= PastixFactLU;
	iparm[IPARM_TASKS2D_WIDTH] 			= globDoublePrecision ? 256 : 128;

	iparm[IPARM_REUSE_LU] 				= firstIter ? 0 : 1;
	iparm[IPARM_REUSE_LU] 				= forceRedo ? 2 : 1;
	
    iparm[IPARM_GPU_MEMORY_PERCENTAGE] 	= 95;
    iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] 	= 64 * 1024;

    
    dparm[DPARM_EPSILON_REFINEMENT] 	= 1e-12;
    iparm[IPARM_ITERMAX]            	= 50;
    iparm[IPARM_GMRES_IM]            	= 50;

	// Initialize sparse matrix
	spm = malloc( sizeof( spmatrix_t ) );
	spmInit(spm);
	spm->flttype = globDoublePrecision ? SpmDouble : SpmFloat;
	spm->values = aupastix;
	spm->fmttype = SpmCSC;
	spm->nexp = spm->gNexp = spm->gN = spm->n = *neq;
    spm->mtxtype = SpmGeneral;
    if( *inputformat == 3 ){
	    spm->nnzexp = spm->gnnzexp = spm->gnnz = spm->nnz = nzsTotal + *neq;
        spm->flttype = SpmDouble;
    } else{
	    spm->nnzexp = spm->gnnzexp = spm->gnnz = spm->nnz = nzsTotal * 2 + *neq;
    }
	spm->colptr = (spm_int_t*) icolpastix;
	spm->rowptr = (spm_int_t*) irowpastix;

    if(*inputformat == 3){
        spmatrix_t *spm_tmp=NULL;
        spm_tmp = spmCopy(spm);
        spmCheckAndCorrect( spm_tmp, spm);
    	spm->flttype = globDoublePrecision ? SpmDouble : SpmFloat;
    }

	// initialize pastix
	pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

	printf("\n");
	spmPrintInfo( spm, stdout );
	printf("\n");
	
	// perform reordering, analysis and symbolic factorization
	pastix_task_analyze( pastix_data, spm );
}

void pastix_csc_conversion(double *ad, double *au, double *adb, double *aub, 
double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){
	
	ITG i,j;
    char merged=0;
	
	// jq for the merged matrix
	ITG* jqTotal = NULL;
	
	redo = forceRedo ? 1 : 0;
	
	if(*neq != neqPrev || *inputformat != inputformatPrev)
		redo = 1;
	
	if(!redo){
		nzsTotal = 0;
		NNEW(icolTotal, ITG, *neq);
        ITG base = jq[0];
		
		// Compute the number of entries in the merged matrix
		#pragma omp parallel for reduction(+:nzsTotal)
		for(i=0;i<*neq;i++){
			ITG kCur = jq[i] - base;
			ITG kPrev = jqPrev[i] - basePrev;
			ITG curColTotal = 0;
			while(kCur < jq[i+1] - base && kPrev < jqPrev[i+1] - basePrev) {
				if(irowPrev[kPrev] == irow[kCur]){
					kCur++;
					kPrev++;
				}
				else{
					if(irowPrev[kPrev] < irow[kCur])
						kPrev++;
					else // irowPrev[kPrev] > irow[k]
						kCur++;
				}
				curColTotal++;
			}
			while(kCur < jq[i+1] - base){
				kCur++;
				curColTotal++;
			}
			while(kPrev < jqPrev[i+1] - basePrev){
				kPrev++;
				curColTotal++;
			}
			icolTotal[i] = curColTotal;
			nzsTotal += curColTotal;
		}
		
		
		// compute jq for the merged matrix
		NNEW(jqTotal, ITG, (*neq+1));
		jqTotal[0] = base;
		for(i = 0; i < *neq; i++){
			jqTotal[i+1] = jqTotal[i] + icolTotal[i];
		}

		
		// If the number of entries in the merged matrix is the same as in the last iteration, we can reuse
		if(nzsTotal == nzsPrev){
			printf("Reusing csc.!\n");
		}
		else{
			redo = 1;
			printf("Not reusing csc, merging patterns!\n");
		}

		// allocate space for the sparse matrix
        if(*symmetryflag || *inputformat != 3)
			NNEW(auPtr,double,2 * nzsTotal);
		else
			NNEW(auPtr,double,nzsTotal);

		NNEW(adPtr,double,neqPrev);
		NNEW(irowTotal, ITG, nzsTotal);

        if(*symmetryflag || *inputformat != 3){
            j=2*nzsTotal;
        }
        else{
            j=nzsTotal;
        }
        #pragma omp parallel for
        for(i=0;i<j;i++){
            auPtr[i] = 0.0;
        }

        #pragma omp parallel for
        for(i=0;i<*neq;i++){
            adPtr[i] = 0.0;
        }
		
		// merge the old and the new sparsity pattern
		#pragma omp parallel for
		for(i=0;i<*neq;i++){
			ITG kCur = jq[i] - base;
			ITG kPrev = jqPrev[i] - base;
			ITG kTotal = jqTotal[i] - base;
			adPtr[i] = ad[i] - (*sigma == 0 ? 0.0 : (*sigma)*adb[i]);
			while(kCur < jq[i+1] - base && kPrev < jqPrev[i+1] - base) {
				if(irowPrev[kPrev] == irow[kCur]){
					auPtr[kTotal] = au[kCur] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur]);

                    if(*symmetryflag && *inputformat != 3 )
						auPtr[kTotal + nzsTotal] = au[kCur + *nzs3] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur + *nzs3]);
					
					irowTotal[kTotal] = irow[kCur];
					kCur++;
					kPrev++;
				}
				else{
					if(irowPrev[kPrev] < irow[kCur]){
//						auPtr[kTotal] = 0.0; 
//						if(*symmetryflag)
//							auPtr[kTotal + nzsTotal] = 0.0;
							
						irowTotal[kTotal] = irowPrev[kPrev];
						
						kPrev++;
					}
					else // irowPrev[kPrev] > irow[k]
					{
						auPtr[kTotal] = au[kCur] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur]);

                        if(*symmetryflag && *inputformat != 3)
							auPtr[kTotal + nzsTotal] = au[kCur + *nzs3] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur + *nzs3]);
						
						irowTotal[kTotal] = irow[kCur];
						kCur++;
					}
				}
				kTotal++;
			}
			while(kCur < jq[i+1] - base){
				auPtr[kTotal] = au[kCur] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur]);

                if(*symmetryflag && *inputformat != 3)
					auPtr[kTotal + nzsTotal] = au[kCur + *nzs3] - (*sigma == 0 ? 0.0 : (*sigma)*aub[kCur + *nzs3]);
						
				irowTotal[kTotal] = irow[kCur];
					
				kCur++;
				kTotal++;
			}
			while(kPrev < jqPrev[i+1] - base){
//				auPtr[kTotal] = 0.0;
//				if(*symmetryflag)
//					auPtr[kTotal + nzsTotal] = 0.0;
						
				irowTotal[kTotal] = irowPrev[kPrev];
					
				kPrev++;
				kTotal++;
			}
		}
		
		SFREE(irowPrev);
		SFREE(icolPrev);
        SFREE(jqPrev);
		
		
		// update pointers to the merged matrix
		icol = icolTotal;
		icolPrev = icolTotal;
		irow = irowTotal;
		irowPrev = irowTotal;
        jqPrev = jqTotal;
		nzsPrev = nzsTotal;
        basePrev = base;
		
		au = auPtr;
		ad = adPtr;
        merged = 1;
		
			
	}
	else
	{
		// This is executed in either the first iteration, or when the number of equations changed
	
		printf("Not reusing csc.\n");
		if(icolPrev != NULL){
			SFREE(icolPrev);
			SFREE(irowPrev);
            SFREE(jqPrev);
		}

		NNEW(icolPrev,ITG,*neq);
		NNEW(irowPrev,ITG,*nzs);
		NNEW(jqPrev,ITG,*neq+1);
		
		memcpy(icolPrev, icol, sizeof(ITG) * *neq);
		memcpy(irowPrev, irow, sizeof(ITG) * *nzs);
		memcpy(jqPrev,   jq,   sizeof(ITG) * (*neq+1));
		
		nzsTotal = *nzs;
		nzsPrev = *nzs;
		neqPrev = *neq;
        jqTotal = jqPrev;
        inputformatPrev = *inputformat;
	}


	// Convert Matrix Format
	if(*inputformat==1 || *symmetryflag==0){
	
	/* lower triangular matrix is stored column by column in
	 au, followed by the upper triangular matrix row by row;
	 the diagonal terms are stored in ad */

		// allocate space for the matrix and the utility arrays
		if(redo){
			// We allocate 10% more space for the values than required so that we have to perform the expensive cudaMallocHost only once, even when the size of the matrix increases slightly
			if((nzsTotal * 2 + *neq) > pastix_nnzBound){
				// perform the call with PaStiX because pinned memory allocation via CUDA is performed if gpu is activated
				pastixAllocMemory((void**)&aupastix, sizeof(double) * 1.1 * (nzsTotal * 2 + *neq), gpu);
				pastix_nnzBound = 1.1 * (nzsTotal * 2 + *neq);
			}
			NNEW(irowpastix,ITG,nzsTotal*2+*neq);
			NNEW(icolpastix,ITG,*neq+1);
			
			if(irowacc != NULL){
				SFREE(irowacc);		
			}
			NNEW(irowacc,ITG,*neq);		
			if(irowPrediction != NULL){
				SFREE(irowPrediction);
			}
			NNEW(irowPrediction,ITG,nzsTotal);
		}


		// Compute utility pointers for parallelization
		// irowPrediction stores the offset to the first entry in it's column of each entry 
		// irowacc stores the number of elements in each row
		if(redo){
			for(i=0;i<nzsTotal;i++){
				irowPrediction[i] = irowacc[irow[i]-1]++;
			}
	
			icolpastix[0] = 1;
			for(i=0;i<*neq;i++){
				icolpastix[i+1] = icolpastix[i] + icol[i] + irowacc[i] + 1;
			}
		}
		
		// copy lower triangular values to the right position in the CSC
		#pragma omp parallel for private(i, j)
		for(i=0;i<*neq;i++){
			ITG k_pastix = icolpastix[i] + irowacc[i];
			ITG k = jqTotal[i] - 1;
			aupastix[k_pastix-1] = ad[i] - (merged != 0 ? 0 : (*sigma == 0 ? 0.0 : (*sigma)*adb[i]));
   			memcpy(aupastix + k_pastix, au + k, sizeof(double) * icol[i]);
            if(*sigma != 0.0 && !merged ){
                for(j=0;j<icol[i];j++){
                    aupastix[k_pastix+j] -= (*sigma)*aub[k+j];
                }
            }
                
		}
		
		// copy the upper triangular values to the right position in the CSC
		#pragma omp parallel for private(i, j)
		for(i=0;i<*neq;i++){
			ITG k = jqTotal[i] - 1;
			for(j=0;j<icol[i];j++){
    			aupastix[irowPrediction[k] + icolpastix[irow[k]-1] - 1] = au[k+(*symmetryflag == 0 ? 0 : (*nzs == *nzs3 ? nzsTotal : *nzs3))] - (merged != 0 ? 0 : (*sigma == 0 ? 0.0 : (*sigma *aub[k+(*symmetryflag == 0 ? 0 : (*nzs == *nzs3 ? nzsTotal : *nzs3))])));
				k++;
			}
		}
		
		// do the same for the rowptr (does not change when reusing)
		if(redo){
		
			#pragma omp parallel for
			for(i=0;i<*neq;i++){
				ITG k_pastix = icolpastix[i] + irowacc[i];
				ITG k = jqTotal[i] - 1;
				irowpastix[k_pastix-1] = i+1;
				memcpy(irowpastix + k_pastix, irow + k, sizeof(ITG) * icol[i]);
			}
			
			#pragma omp parallel for private(i, j)
			for(i=0;i<*neq;i++){
				ITG k = jqTotal[i] - 1;
				for(j=0;j<icol[i];j++){
					irowpastix[irowPrediction[k] + icolpastix[irow[k]-1] - 1] = i+1;
					k++;
				}
			}
		}
	}
	else if(*inputformat==3){

        if(redo){
       		NNEW(irowpastix,ITG,nzsTotal+*neq);
    		NNEW(icolpastix,ITG,*neq+1);
       		if((nzsTotal + *neq) > pastix_nnzBound){
            	pastixAllocMemory((void**)&aupastix, sizeof(double) * 1.1 * (nzsTotal + *neq), gpu);
   				pastix_nnzBound = 1.1 * (nzsTotal + *neq);
            }
            memcpy(icolpastix, jqTotal, sizeof(ITG) * (*neq+1));
        }


        ITG kflag,ilength;
        kflag = 2;

        #pragma omp parallel for private(i, j)
        for(i=0;i<*neq;i++){
            if(redo){
                memcpy(irowpastix + i + (jqTotal[i] - jqTotal[0]), irow + (jqTotal[i]-jqTotal[0]), sizeof(ITG) * (jqTotal[i+1] - jqTotal[i]));
                irowpastix[i+jqTotal[i+1]-1] = i+1;
                icolpastix[i+1] += i+1;
            }
            memcpy(aupastix + i + (jqTotal[i] - jqTotal[0]), au + (jqTotal[i]-jqTotal[0]), sizeof(double) * (jqTotal[i+1] - jqTotal[i]));
            if(*sigma != 0 && merged == 0){
                for(j=0;j<jqTotal[i+1]-jqTotal[i];j++){
                    aupastix[i+(jqTotal[i]-jqTotal[0])+j] -= (*sigma)*aub[jqTotal[i]-jqTotal[0]+j];
                }
            }
            aupastix[i+jq[i+1]-1] = ad[i] - (merged != 0 ? 0.0 : (*sigma == 0 ? 0.0 : (*sigma)*adb[i]));

        }

        if(redo) icolpastix[*neq] = nzsTotal + *neq + 1;
	} 


	// Free the merged array in STI format
	if(auPtr){
		SFREE(auPtr);
		auPtr = NULL;
	}
	if(adPtr){
		SFREE(adPtr);
		adPtr = NULL;
	}
	
}



// PaStiX invocation when the factorization function is called individually
void pastix_factor_main_generic(double *ad, double *au, double *adb, double *aub, 
double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){
			
    pastix_set_globals(mode);	
    gpu = 0;
	// Set GPU flag from environment
//	const char* pastix_gpu = getenv("PASTIX_GPU");
//	if(pastix_gpu)
//		//gpu = (*pastix_gpu == '1') ? 1 : 0;
//		gpu = ( mode == AS ) ? 0 : (*pastix_gpu == '1') ? 1 : 0;
	
	// Perform individual invocations always in double precision. If previous iterations were in single precision, do not reuse.
	forceRedo=1;
	globDoublePrecision = 1;

	// invoke PaStiX
	pastix_csc_conversion(ad, au, adb, aub, sigma, icol, irow, neq, nzs, symmetryflag, inputformat, jq, nzs3);
	pastix_cleanup(neq,symmetryflag);
	pastix_init(ad, au, adb, aub, sigma, icol, irow, neq, nzs, symmetryflag, inputformat, jq, nzs3);	
	pastix_factor(ad, au, adb, aub, sigma, icol, irow, neq, nzs, symmetryflag, inputformat, jq, nzs3);
}

// invokes the factorization routine of PaStiX
void pastix_factor(double *ad, double *au, double *adb, double *aub, 
double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){
			
	if(spm->n == 1)
		return;
		
	pastix_task_numfact( pastix_data, spm );
}


// invokes the solve and iterative refinement routines of PaStiX
ITG pastix_solve_generic(double *x, ITG *neq,ITG *symmetryflag,ITG *nrhs){
	ITG i;
	double* b;
	float* buffer;
	
	// dont call pastix with only one equation, might lead to segfault
	if(spm->n == 1)
	{
		x[0] = x[0] / aupastix[0];
		return 0;
	}
	
	
	// check whether the RHS consists of only Zeroes and return in that case
	char allZero = 1;
	for(i = 0; i < *neq; i++){
		if(x[i] != 0){
			allZero = 0;
			break;
		}
	}
	if(allZero){
		printf("RHS only consists of 0.0\n");
		return 0;
	}

	//Copy the b so that we can modify x without losing b
	NNEW(b,double,*nrhs**neq);
	memcpy(b, x, sizeof(double) * (*nrhs) * (*neq));
	

	// If we are in mixed precision mode, cast double x to float x and call solve. Afterwards upcast the solution.
	if(!globDoublePrecision){
		NNEW(buffer,float,*nrhs**neq);
		#pragma omp parallel for
		for(i = 0; i < (*nrhs) * (*neq); i++){ 
			buffer[i] = (float) x[i];
		}
		
		pastix_task_solve( pastix_data, *nrhs, buffer, spm->n );

		#pragma omp parallel for
		for(i = 0; i < (*nrhs) * (*neq); i++){
			x[i] = (double) buffer[i];
		}
		SFREE(buffer);
	}
	else{
		pastix_task_solve( pastix_data, *nrhs, x, spm->n );
	}
	
	// invoke iterative refinement in double precision
	ITG rc = pastix_task_refine( pastix_data, spm->n, *nrhs, (void*)b, spm->n, (void*)x, spm->n );

//    FILE *f=fopen("spm.out","a");
//    fprintf(f,"\n\nMatrix\n");
//    spmConvert(SpmCSR, spm);
//    spmPrint( spm, f);
//    spmConvert(SpmCSC, spm);
//
//    fprintf(f,"Solution vector b:\n");
//    for(i=0;i<*neq;i++){
//        fprintf(f,"b[%d] = %.17f\n",i,b[i]);
//    }
//
//    fprintf(f,"Solution vector x:\n");
//    for(i=0;i<*neq;i++){
//        fprintf(f,"x[%d] = %.17f\n",i,x[i]);
//    }
//
//    fclose(f);

	SFREE(b);

    modePrev = mode;

	return rc;
}

ITG pastix_solve(double *x, ITG *neq,ITG *symmetryflag,ITG *nrhs){
    mode = BASIC;
    if( modePrev != mode ) pastix_set_globals(mode);	
    ITG rc = pastix_solve_generic(x,neq,symmetryflag,nrhs);
    return rc;
}

ITG pastix_solve_as(double *x, ITG *neq,ITG *symmetryflag,ITG *nrhs){
    mode = AS;
    if( modePrev != mode ) pastix_set_globals(mode);	
    ITG rc = pastix_solve_generic(x,neq,symmetryflag,nrhs);
    return rc;
}

// Invokes pastixFinalize and spmExit which frees everything but the dense LU array and parsec pointer
void pastix_cleanup(ITG *neq,ITG *symmetryflag){
	if( redo && !firstIter ){
		spmExit( spm );
		if(spm != NULL) free( spm );
			
		pastixFinalize( &pastix_data );
	}
	
	firstIter = 0;
	return;
}

void pastix_cleanup_as(ITG *neq,ITG *symmetryflag){
    pastix_cleanup(neq,symmetryflag);
    return;
}

// main method for executing PaStiX
void pastix_main_generic(double *ad, double *au, double *adb, double *aub, 
     double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	 ITG *jq, ITG *nzs3,ITG *nrhs){
		 
		 
	if(*neq==0) return;

    pastix_set_globals( mode );	
	
	const char* pastix_gpu = getenv("PASTIX_GPU");
	if(pastix_gpu)
		gpu = (*pastix_gpu == '1') ? 1 : 0;
	
	// check mixed precision environment variable
    const char* pastix_mixed = getenv("PASTIX_MIXED_PRECISION");	
    if( pastix_mixed && globDoublePrecision == 0 ){
        mixed = (*pastix_mixed == '1') ? 1 : 0;
        globDoublePrecision = (mixed == 0) ? 1 : 0;
    }
    else{
        mixed = 1;
    }

	if(!stickToDouble && mixed == '1'){
		if(globDoublePrecision == 1)
			forceRedo = 1;

		globDoublePrecision = 0;
	}


    if(*inputformat==3){
        forceRedo=1;
    }

	// backup b in case mixed precision solve corrupts the original array
    double* b_backup = NULL;
    NNEW(b_backup, double, *nrhs * *neq);
    memcpy(b_backup, b, sizeof(double) * (*nrhs)*(*neq));

	// benchmarking structs
	struct timespec start, end; 
	struct timespec stepCscConversionStart, stepCscConversionEnd; 
	struct timespec stepInitStart, stepInitEnd; 
	struct timespec stepFactorizeStart, stepFactorizeEnd; 
	struct timespec stepSolveStart, stepSolveEnd; 
	struct timespec stepCleanUpStart, stepCleanUpEnd; 
	double pastixTime, stepCscConversion, stepInit, stepFactorize, stepSolve, stepCleanUp, totalCCXTime, CCXwithoutPastix; 

	clock_gettime(CLOCK_MONOTONIC, &start); 
	clock_gettime(CLOCK_MONOTONIC, &stepCscConversionStart); 
	
	// invoke csc conversion
	pastix_csc_conversion(ad,au,adb,aub,sigma,icol,irow, 
			 neq,nzs,symmetryflag,inputformat,jq,nzs3);
			 
	clock_gettime(CLOCK_MONOTONIC, &stepCscConversionEnd); 
	clock_gettime(CLOCK_MONOTONIC, &stepCleanUpStart); 

	// invoke cleanup
	pastix_cleanup(neq,symmetryflag);
			 
	clock_gettime(CLOCK_MONOTONIC, &stepCleanUpEnd); 
	clock_gettime(CLOCK_MONOTONIC, &stepInitStart); 
	
	//invoke init
	pastix_init(ad,au,adb,aub,sigma,icol,irow, 
		 neq,nzs,symmetryflag,inputformat,jq,nzs3);
		 
	clock_gettime(CLOCK_MONOTONIC, &stepInitEnd); 
	clock_gettime(CLOCK_MONOTONIC, &stepFactorizeStart); 
 
	// invoke factor
	pastix_factor(ad,au,adb,aub,sigma,icol,irow, 
		 neq,nzs,symmetryflag,inputformat,jq,nzs3);
 
	clock_gettime(CLOCK_MONOTONIC, &stepFactorizeEnd); 
	clock_gettime(CLOCK_MONOTONIC, &stepSolveStart); 
	
	// if solve does not converge
	if(pastix_solve_generic(b,neq,symmetryflag,nrhs) == -1){
		
		// Give up, if we tried it with double precision, use backup b otherwise
		if(globDoublePrecision == 1){
            printf("PaStiX could not converge to a valid result\n");
			exit(1);
		}
        else{
		    memcpy(b, b_backup, sizeof(double) * (*nrhs)*(*neq));
        }
        
        // make sure that we switch to double and do not reuse in the next iteration
        globDoublePrecision = 1;
		forceRedo = 1;
		stickToDouble = 1;
		mixedFailed++;
		
	 	// call pastix_main recursively, but now in double precision
		pastix_main_generic(ad, au, adb, aub, sigma, b, icol, irow, neq, nzs, symmetryflag, inputformat, jq, nzs3, nrhs);
		
		// if we do not converge with mixed precision for the third time, permanently switch to double precision
		if(mixedFailed <= 2){
			stickToDouble = 0;
			forceRedo = 1;
		}
			
		return;
	}
	else{
		forceRedo = 0;
	}
	
	clock_gettime(CLOCK_MONOTONIC, &stepSolveEnd); 
	clock_gettime(CLOCK_MONOTONIC, &end);

	// compute benchmark times
	
	pastixTime = (end.tv_sec - start.tv_sec) * 1e9; 
	pastixTime = (pastixTime + (end.tv_nsec - start.tv_nsec)) * 1e-9;

	totalPastixTime += pastixTime;

	clock_gettime(CLOCK_MONOTONIC, &totalCalculixTimeEnd); 

	totalCCXTime = (totalCalculixTimeEnd.tv_sec - totalCalculixTimeStart.tv_sec) * 1e9; 
	totalCCXTime = (totalCCXTime + (totalCalculixTimeEnd.tv_nsec - totalCalculixTimeStart.tv_nsec)) * 1e-9;

	CCXwithoutPastix = totalCCXTime - totalPastixTime;

	stepCscConversion = (stepCscConversionEnd.tv_sec - stepCscConversionStart.tv_sec) * 1e9; 
	stepCscConversion = (stepCscConversion + (stepCscConversionEnd.tv_nsec - stepCscConversionStart.tv_nsec)) * 1e-9;

	stepInit = (stepInitEnd.tv_sec - stepInitStart.tv_sec) * 1e9; 
	stepInit = (stepInit + (stepInitEnd.tv_nsec - stepInitStart.tv_nsec)) * 1e-9;

	stepFactorize = (stepFactorizeEnd.tv_sec - stepFactorizeStart.tv_sec) * 1e9; 
	stepFactorize = (stepFactorize + (stepFactorizeEnd.tv_nsec - stepFactorizeStart.tv_nsec)) * 1e-9;

	stepSolve = (stepSolveEnd.tv_sec - stepSolveStart.tv_sec) * 1e9; 
	stepSolve = (stepSolve + (stepSolveEnd.tv_nsec - stepSolveStart.tv_nsec)) * 1e-9;

	stepCleanUp = (stepCleanUpEnd.tv_sec - stepCleanUpStart.tv_sec) * 1e9; 
	stepCleanUp = (stepCleanUp + (stepCleanUpEnd.tv_nsec - stepCleanUpStart.tv_nsec)) * 1e-9;
	
	
	// upate iteration timer
	totalIterations++;
	if(!redo)
		totalReused++;

	// benchmark output

	printf("________________________________________\n\n");

	printf("CSC Conversion Time: %lf\n", stepCscConversion);
	printf("Init Time: %lf\n", stepInit);
	printf("Factorize Time: %lf\n", stepFactorize);
	printf("Solve Time: %lf\n", stepSolve);
	printf("Clean up Time: %lf\n", stepCleanUp);
	printf("---------------------------------\n");
	printf("Sum: %lf\n", pastixTime);
	printf("\n");
	printf("Total PaStiX Time: %lf\n", totalPastixTime);
	printf("CCX without PaStiX Time: %lf\n", CCXwithoutPastix);
	printf("Share of PaStiX Time: %lf\n", totalPastixTime / totalCCXTime );
	printf("Total Time: %lf\n", totalCCXTime);
	printf("Reusability: %d : %d \n", totalReused, totalIterations);

	printf("________________________________________\n\n");
	
	
    SFREE(b_backup);
	return;
	
	
}


void pastix_factor_main(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){

        mode = BASIC;
        pastix_factor_main_generic(ad, au, adb, aub, sigma, icol, irow, neq,
                                    nzs, symmetryflag, inputformat, jq, nzs3);

        return;
}

void pastix_factor_main_as(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){

        mode = AS;
        pastix_factor_main_generic(ad, au, adb, aub, sigma, icol, irow, neq,
                                    nzs, symmetryflag, inputformat, jq, nzs3);

        return;
}

void pastix_main(double *ad, double *au, double *adb, double *aub, 
     double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	 ITG *jq, ITG *nzs3,ITG *nrhs){
	 
     mode = BASIC;
     pastix_main_generic(ad, au, adb, aub, sigma, b, icol, irow, neq,
                          nzs, symmetryflag, inputformat, jq, nzs3, nrhs);

     return;
}

void pastix_main_as(double *ad, double *au, double *adb, double *aub, 
     double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	 ITG *jq, ITG *nzs3,ITG *nrhs){
	 
     mode = AS;
     pastix_main_generic(ad, au, adb, aub, sigma, b, icol, irow, neq,
                          nzs, symmetryflag, inputformat, jq, nzs3, nrhs);

     return;
}

void pastix_set_globals(char mode){

    if( modePrev != mode ){
        pastix_data_object *temp,*temp2;
        switch(mode){
            case BASIC:
                temp = &pastix_mode_basic;
                temp2 = &pastix_mode_as;
                break;

            case AS:
                temp = &pastix_mode_as;
                temp2 = &pastix_mode_basic;
                break;
        }

        // saving old data set
        temp2->totalIterations = totalIterations;
        temp2->totalReused = totalReused;
        
        // Current sparse matrix in STI format
        temp2->icolTotal = icolTotal;
        temp2->irowTotal = irowTotal;
        
        // Matrix data from previous iteration
        temp2->neqPrev = neqPrev;
        temp2->nzsPrev = nzsPrev;
        temp2->icolPrev = icolPrev;
        temp2->irowPrev = irowPrev;
        temp2->jqPrev = jqPrev;
        temp2->inputformatPrev = inputformatPrev;
        temp2->basePrev = basePrev;
        
        // Current sparse matrix in CSC
        temp2->aupastix = aupastix;
        temp2->icolpastix = icolpastix;
        temp2->irowpastix = irowpastix;
        
        // Global variable that indicates whether we are currently reusing or not
        temp2->redo = redo;
        
        // PaStiX configuration
        temp2->iparm = iparm;
        temp2->dparm = dparm;
        temp2->pastix_data = pastix_data;
        temp2->spm = spm;
        
        // GPU active or not
        temp2->gpu = gpu;

        // Store how many nzs the merged Matrix has
        temp2->nzsTotal = nzsTotal;
        
        // Size of allocated space for sparse matrix
        temp2->pastix_nnzBound = pastix_nnzBound;
        
        // Number of iterations that failed with mixed precision
        temp2->mixedFailed = mixedFailed;
        
        // indicates whether this is the first invocation of PaStiX or not
        temp2->firstIter = firstIter;
        
        // When this flag is activated, PaStiX will not reuse in the next iteration
        temp2->forceRedo = forceRedo;
        
        // Use double or mixed precision
        temp2->mixed = mixed;
        temp2->globDoublePrecision = globDoublePrecision;
        
        // This is set to one, when to many iterations with mixed precision did not converge
        temp2->stickToDouble = stickToDouble;
        
        // Pointers for faster matrix transpose
        temp2->irowacc = irowacc;
        temp2->irowPrediction = irowPrediction;
        
        // Number of threads
        // nthread_mkl=0;


        // setting new data set
        totalIterations = temp->totalIterations;
        totalReused = temp->totalReused;
        
        // Current sparse matrix in STI format
        icolTotal = temp->icolTotal;
        irowTotal = temp->irowTotal;
        
        // Matrix data from previous iteration
        neqPrev = temp->neqPrev;
        nzsPrev = temp->nzsPrev;
        icolPrev = temp->icolPrev;
        irowPrev = temp->irowPrev;
        jqPrev = temp->jqPrev;
        inputformatPrev = temp->inputformatPrev;
        basePrev = temp->basePrev;
        
        // Current sparse matrix in CSC
        aupastix = temp->aupastix;
        icolpastix = temp->icolpastix;
        irowpastix = temp->irowpastix;
        
        // Global variable that indicates whether we are currently reusing or not
        redo = temp->redo;
        
        // PaStiX configuration
        iparm = temp->iparm;
        dparm = temp->dparm;
        pastix_data = temp->pastix_data;
        spm = temp->spm;
        
        // GPU active or not
        gpu = temp->gpu;
        
        // Store how many nzs the merged Matrix has
        nzsTotal = temp->nzsTotal;
        
        // Size of allocated space for sparse matrix
        pastix_nnzBound = temp->pastix_nnzBound;
        
        // Number of iterations that failed with mixed precision
        mixedFailed = temp->mixedFailed;
        
        // indicates whether this is the first invocation of PaStiX or not
        firstIter = temp->firstIter;
        
        // When this flag is activated, PaStiX will not reuse in the next iteration
        forceRedo = temp->forceRedo;
        
        // Use double or mixed precision
        mixed = temp->mixed;
        globDoublePrecision = temp->globDoublePrecision;
        
        // This is set to one, when to many iterations with mixed precision did not converge
        stickToDouble = temp->stickToDouble;
        
        // Pointers for faster matrix transpose
        irowacc = temp->irowacc;
        irowPrediction = temp->irowPrediction;
        
        // Number of threads
        // nthread_mkl=0;
        if( firstIter ){
            forceRedo=1;
        }

    }

    modePrev = mode;

    return;
}

#endif

