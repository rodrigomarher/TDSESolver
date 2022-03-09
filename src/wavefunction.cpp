#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "utils.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters *param){
	_param = param;
}

#ifdef MPI
void WF::set_mpi(mpi_grid *grid){
	_mpi_grid = grid;	
    _nproc_i = _mpi_grid->dims[0];
	_nproc_j = _mpi_grid->dims[1];
	//if(_mpi_grid->rank == 0){
	//	std::cout<<"Grid: ("<<_nproc_i<<","<<_nproc_j<<")"<<std::endl;
	//}
    //MPI_Barrier(_mpi_grid->comm);	
	//std::cout<<"Node number: "<<_mpi_grid->rank<<" Coords: "<<_mpi_grid->coords[0]<<" "<<_mpi_grid->coords[1]<<std::endl; 
    //MPI_Barrier(_mpi_grid->comm);	
}
#endif

void WF::set_geometry( double *i, double *j, double *k, const double di, const double dj, const double dk){    
    _ni = _param->ni;
    _nj = _param->nj;
    _nk = _param->nk;
    _wf = alloc3d<cdouble>(_ni/_nproc_i, _nj/_nproc_j, _nk);
    if(_param->geometry == XYZ)
        _wf_buf = alloc4d<cdouble>(_ni/_nproc_i, _nj/_nproc_j, _nk, 1);
    else
        _wf_buf = alloc4d<cdouble>(_ni/_nproc_i, _nj/_nproc_j, _nk, _param->nt_diag);
    _i_row = new cdouble[_ni];
    _j_row = new cdouble[_nj];
    _k_row = new cdouble[_nk];

    _diag_buf = new cdouble[_param->nt_diag];

    for(int i=0; i<_ni/_nproc_i; i++){
        for(int j=0; j<_nj/_nproc_j;j++){
            for(int k=0;k<_nk;k++){
                _wf[i][j][k] = 0.0;
            }
        }
    }

    _i = i; _j = j; _k = k; _di = di; _dj = dj; _dk = dk;

    switch(_param->geometry){
        case X:
            _geom_X();
            break;

        case XZ:
            _geom_XZ();
            break;
        case RZ:
            _geom_RZ();
            break;
        case XYZ:
            _geom_XYZ();
            break;
    }
}

void WF::gaussian(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni/_nproc_i; i++){
        for(int j=0; j<_nj/_nproc_j;j++){
			int idx_i = _mpi_grid->coords[0]*_ni/_nproc_i + i;
			int idx_j = _mpi_grid->coords[1]*_nj/_nproc_j + j;
            for(int k=0; k<_nk; k++){
                _wf[i][j][k] = exp(-(_i[idx_i]-i0)*(_i[idx_i]-i0)/sigma- (_j[idx_j]-j0)*(_j[idx_j]-j0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
            }
        }
    }

}

void WF::exponential(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni/_nproc_i; i++){
        for(int j=0; j<_nj/_nproc_j;j++){
			int idx_i = _mpi_grid->coords[0]*_ni/_nproc_i + i;
			int idx_j = _mpi_grid->coords[1]*_nj/_nproc_j + j;
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] = exp(-sqrt((_i[idx_i]-i0)*(_i[idx_i]-i0) + (_j[idx_j]-j0)*(_j[idx_j]-j0) + (_k[k]-k0)*(_k[k]-k0)));
            }
        }
    }
}

cdouble WF::norm(){
    cdouble norm2 = 0.0;
    switch(_param->geometry){
        case X:
			norm2 = _norm2_X();
            break;
        case XZ:
			norm2 = _norm2_XZ();
            break;
        case RZ:
			norm2 = _norm2_RZ();
            break;
        case XYZ:
            norm2 = _norm2_XYZ(); 
            break;
    }
    return sqrt(norm2);
}

void WF::apply_mask(cdouble *imask, cdouble *jmask, cdouble *kmask){
    (this->*(this->_apply_mask))(imask, jmask, kmask);
}


cdouble*** WF::get(){
    return _wf;
}

cdouble* WF::i_row(int j, int k){
	for(int i=0;i<_ni;i++){
        _i_row[i] = _wf[i][j][k];
    }
    return _i_row;
}

cdouble* WF::j_row(int i, int k){
    for(int j=0; j<_nj; j++){
        _j_row[j] = _wf[i][j][k];
    }
    return _j_row;
}

cdouble* WF::k_row(int i, int j){
    for(int k=0; k<_nk;k++){
        _k_row[k] = _wf[i][j][k];
    }
    return _k_row;
}

void WF::set_i_row(cdouble* i_row, int j, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][j][k] = i_row[i];
    }
}

void WF::set_i_row_mask(cdouble* i_row, cdouble* imask, int j, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][j][k] = i_row[i]*imask[i];
    }
}

void WF::set_j_row(cdouble* j_row, int i, int k){
    for(int j=0; j<_nj;j++){
        _wf[i][j][k] = j_row[j];
    }
}

void WF::set_j_row_mask(cdouble* j_row, cdouble *jmask, int i, int k){
    for(int j=0; j<_nj;j++){
        _wf[i][j][k] = j_row[j]*jmask[j];
    }
}

void WF::set_k_row(cdouble* k_row, int i, int j){
    for(int k=0; k<_nk; k++){
        _wf[i][j][k] = k_row[k];
    }
}

void WF::set_k_row_mask(cdouble* k_row, cdouble *kmask, int i, int j){
    for(int k=0; k<_nk; k++){
        _wf[i][j][k] = k_row[k]*kmask[k];
    }
}

void WF::get_i_row(cdouble* i_row, int j, int k){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf[i][j][k];
    }
}

void WF::get_j_row(cdouble* j_row, int i, int k){
    for(int j=0; j<_nj;j++){
        j_row[j] = _wf[i][j][k];
    }
}

void WF::get_k_row(cdouble* k_row, int i, int j){
    for(int k=0; k<_nk; k++){
        k_row[k] = _wf[i][j][k];
    }
}
void WF::set(cdouble*** arr){
    for(int i=0; i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] = arr[i][j][k];
            }
        }
    }
}

void WF::set_to_buf(const int idx){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf_buf[idx][i][j][k] = _wf[i][j][k];
            }
        }
    }
    //std::memcpy(_wf_buf[idx],_wf,_ni*_nk*sizeof(cdouble));
}

void WF::set_i_row_buf(cdouble* i_row, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        _wf_buf[idx][i][j][k] = i_row[i];
    }
}

void WF::set_j_row_buf(cdouble *j_row, const int i, const int k, const int idx){
    for(int j=0;j<_nj;j++)
        _wf_buf[idx][i][j][k] = j_row[j];
}

void WF::set_k_row_buf(cdouble *k_row, const int i, const int j, const int idx){
    for(int k=0; k<_nk;k++)
        _wf_buf[idx][i][j][k] = k_row[k];
}

void WF::get_i_row_buf(cdouble* i_row, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf_buf[idx][i][j][k];
    }
}

void WF::get_j_row_buf(cdouble *j_row, const int i, const int k, const int idx){
    for(int j=0;j<_nj;j++)
        j_row[j] = _wf_buf[idx][i][j][k];
}

void WF::get_k_row_buf(cdouble *k_row, const int i, const int j, const int idx){
    for(int k=0; k<_nk;k++)
        k_row[k] = _wf_buf[idx][i][j][k];
}


void WF::set_i_row_buf_mask(cdouble* i_row, cdouble* imask, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++)
        _wf_buf[idx][i][j][k] = i_row[i]*imask[i];
}

void WF::set_j_row_buf_mask(cdouble* j_row, cdouble* jmask, const int i, const int k, const int idx){
    for(int j=0; j<_nj; j++)
        _wf_buf[idx][i][j][k] = j_row[j]*jmask[j];
}

void WF::set_k_row_buf_mask(cdouble* k_row, cdouble* kmask, const int i, const int j, const int idx){
    for(int k=0; k<_nk; k++)
        _wf_buf[idx][i][j][k] = k_row[k]*kmask[k];
}



void WF::get_from_buf(cdouble*** arr, const int idx){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0;k<_nk; k++){
                arr[i][j][k] = _wf_buf[idx][i][j][k];
            }
        }
    }
}

cdouble**** WF::get_buf(){
    return _wf_buf;
}

cdouble* WF::get_diag_buf(){
    return _diag_buf;
}

cdouble WF::operator()(int i, int j, int k){
    return _wf[i][j][k];
}

void WF::operator/=(cdouble val){
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0; i<_ni/_nproc_i;i++){
        for(int j=0; j<_nj/_nproc_j;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] /= val;
            }
        }
    }
}

#ifdef MPI
void get_i_row_MPI(cdouble *i_row, int j, int k, int rank){
	
} 
#endif

WF::~WF(){
    free3d(&_wf,_ni/_nproc_i,_nj/_nproc_j,_nk);
    if(_param->geometry == XYZ)
    	free4d(&_wf_buf,_ni/_nproc_i,_nj/_nproc_j, _nk, 1);
    else
    	free4d(&_wf_buf, _ni/_nproc_i, _nj/_nproc_j, _nk, _param->nt_diag);
    delete[] _diag_buf;
    delete[] _i_row;
    delete[] _j_row;
    delete[] _k_row;
}
