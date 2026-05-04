#include "probe.CUSTOM.h"
#include <iostream>

ProbeCUSTOM::ProbeCUSTOM(): Probe(){}

ProbeCUSTOM::ProbeCUSTOM(std::string def): Probe::Probe(def){}

void ProbeCUSTOM::_acc_i(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=1;i<_ni-1;i++){
        for(int j=1; j<_nj-1; j++){
            for(int k=1; k<_nk-1;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                cdouble dV_i = _ham->dpotential_i(_i[i],_j[j],_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*dV_i)*wf[i][j][k]*dr).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_acc_j(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=0;i<_ni;i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                cdouble dV_j = _ham->dpotential_j(_i[i],_j[j],_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*dV_j)*wf[i][j][k]*dr).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_acc_k(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=0;i<_ni;i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                cdouble dV_k = _ham->dpotential_k(_i[i],_j[j],_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*dV_k)*wf[i][j][k]*dr).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_dip_i(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=0;i<_ni;i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*_i[i])*wf[i][j][k]).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_dip_j(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=0;i<_ni;i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*_j[j])*wf[i][j][k]).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_dip_k(const int idx){
    //In 3D is not using buffer for memory economy
    cdouble *** wf = _wf->get();
    cdouble result;
    double sum = 0.0;
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(+: sum)
    for(int i=0;i<_ni;i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk;k++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                sum += (conj(wf[i][j][k])*(-1.0*_k[k])*wf[i][j][k]).real();
            }
        }
    }
    result = sum*_tempmask[idx];
    _data[idx] = result;
}

void ProbeCUSTOM::_pop(const int idx){}
void ProbeCUSTOM::_dens(const int idx){}
void ProbeCUSTOM::_wf_snap(const int idx){}

ProbeCUSTOM::~ProbeCUSTOM(){
    delete[] _data;
}
