#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_XYZ(){
    _apply_mask = &WF::apply_mask_XYZ;
}

void WF::apply_mask_XYZ(cdouble *imask, cdouble *jmask, cdouble *kmask){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf[i][j][k] *= imask[i]*jmask[j]*kmask[k];
            }
        }
    }
}

void WF::apply_mask_buf_XYZ(cdouble *imask, cdouble *jmask, cdouble *kmask, const int idx){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk; k++){
                _wf_buf[idx][i][j][k] *= imask[i]*jmask[j]*kmask[k];
            }
        }
    }

}

cdouble WF::_norm2_XYZ(){
	cdouble integral=0.0;
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for(int i=0;i<_ni/_nproc_i;i++){
		for(int j=0;j<_nj/_nproc_j;j++){
			for(int k=0;k<_nk;k++){
				sum += (_wf[i][j][k]*conj(_wf[i][j][k])).real()*_di*_dk*_dj;
			}
		}
	}

	#ifndef MPI
	integral = sum;
	return integral;
	#endif

	#ifdef MPI
	if(_mpi_grid->rank != 0)
		MPI_Send(&sum,1,MPI_DOUBLE,0,42,_mpi_grid->comm);

	if(_mpi_grid->rank == 0){
		MPI_Request request_arr[_mpi_grid->size];
		MPI_Status  status_arr[_mpi_grid->size];
		double sum_array[_mpi_grid->size];
		sum_array[0] = sum;
		for(int i=1;i<_mpi_grid->size;i++){
			MPI_Irecv(&sum_array[i],1,MPI_DOUBLE,i,42,_mpi_grid->comm,&request_arr[i]);
		}
		MPI_Waitall(_mpi_grid->size-1,request_arr,status_arr);
		for(int i=0;i<_mpi_grid->size;i++){
			integral += sum_array[i];
		}
		
		for(int i= 1; i<_mpi_grid->size;i++)
			MPI_Isend(&integral, 1, MPI_DOUBLE_COMPLEX,i,43,_mpi_grid->comm,&request_arr[i]);
		MPI_Waitall(_mpi_grid->size-1,request_arr,status_arr);
	}
	if(_mpi_grid->rank != 0)
		MPI_Recv(&integral,1,MPI_DOUBLE_COMPLEX,0,43,_mpi_grid->comm,MPI_STATUS_IGNORE);
	return integral;
	#endif
}
