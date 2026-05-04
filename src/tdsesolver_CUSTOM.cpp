#include <iostream>
#include <tuple>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

void TDSESolver::_geom_CUSTOM(){
    _di = 1.0;
    _dj = 1.0;
    _dk = 1.0;

    const int ni = calc_n_elem(_param->file_geom_i);
    if(ni != _param->ni){debug0("[TDSESolver::_geom_CUSTOM] Size of file_geom_i differs from _ni.\n");exit(1);}
    const int nj = calc_n_elem(_param->file_geom_j);
    if(nj != _param->nj){debug0("[TDSESolver::_geom_CUSTOM] Size of file_geom_j differs from _nj.\n");exit(1);}
    const int nk = calc_n_elem(_param->file_geom_k);
    if(nk != _param->nk){debug0("[TDSESolver::_geom_CUSTOM] Size of file_geom_k differs from _nk.\n");exit(1);}
    
    _i = new double[_param->ni];
    _j = new double[_param->nj];
    _k = new double[_param->nk];

    std::ifstream file;
    file.open(_param->file_geom_i);
    std::string line;
    int i = 0;
    if(file.is_open()){
        while(getline(file, line)){
            std::stringstream ss(line);
            ss >> _i[i];
            i++;
        }
        file.close();
    }
    else{debug0("[TDSESolver_geom_CUSTOM] Unable to open file. \n"); exit(1);}
    
    file.open(_param->file_geom_j);
    i = 0;
    if(file.is_open()){
        while(getline(file, line)){
            std::stringstream ss(line);
            ss >> _j[i];
            i++;
        }
        file.close();
    }
    else{debug0("[TDSESolver_geom_CUSTOM] Unable to open file. \n"); exit(1);}

    file.open(_param->file_geom_k);
    i = 0;
    if(file.is_open()){
        while(getline(file, line)){
            std::stringstream ss(line);
            ss >> _k[i];
            i++;
        }
        file.close();
    }
    else{debug0("[TDSESolver_geom_CUSTOM] Unable to open file. \n"); exit(1);}

    _param->imax = _i[_param->ni-1];
    _param->imin = _i[0];
    _param->jmax = _j[_param->nj-1];
    _param->jmin = _j[0];
    _param->kmax = _k[_param->nk-1];
    _param->kmin = _k[0];
  
    _propagate = &TDSESolver::_propagate_CUSTOM;
    _ipropagate = &TDSESolver::_ipropagate_CUSTOM;
}

void TDSESolver::_fields_CUSTOM(){
    if (_param->use_field_file == 0){
        Afield_i = new Field_TDSESolver(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
        Bfield_i = new Field_TDSESolver(_param->B0i, _param->w0Bi, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
        Afield_j = new Field_TDSESolver(_param->E0j, _param->w0Ej, _param->phiEj, _param->env, _param->tmax_ev, _t, _param->nt);
        Bfield_j = new Field_TDSESolver(_param->B0j, _param->w0Bj, _param->phiBj, _param->env, _param->tmax_ev, _t, _param->nt);
        Afield_k = new Field_TDSESolver(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
        Bfield_k =new Field_TDSESolver(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
    }

    else{
        Afield_i = new Field_TDSESolver(_param->file_fieldx, _param->tmax_ev, _t, _param->nt);
        Bfield_i = new Field_TDSESolver(_param->B0i, _param->w0Bi, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
        Afield_j = new Field_TDSESolver(_param->file_fieldy, _param->tmax_ev, _t, _param->nt);
        Bfield_j = new Field_TDSESolver(_param->B0j, _param->w0Bj, _param->phiBj, _param->env, _param->tmax_ev, _t, _param->nt);
        Afield_k = new Field_TDSESolver(_param->file_fieldz, _param->tmax_ev, _t, _param->nt);
        Bfield_k = new Field_TDSESolver(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
    }

    std::string path;
    path = _param->path_results + "/Efield_i.dat";
    write_array(Afield_i->get(),_param->nt,path);
    Afield_i->calc_pot();
    path = _param->path_results + "/Afield_i.dat";
    write_array(Afield_i->get(), _param->nt, path);
    path = _param->path_results + "/Bfield_i.dat";
    write_array(Bfield_i->get(), _param->nt, path);

    path = _param->path_results + "/Efield_j.dat";
    write_array(Afield_j->get(), _param->nt, path);
    Afield_j->calc_pot();
    path = _param->path_results + "/Afield_j.dat";
    write_array(Afield_j->get(), _param->nt, path);
    path = _param->path_results + "/Bfield_j.dat";
    write_array(Bfield_j->get(), _param->nt, path);

    path = _param->path_results + "/Efield_k.dat";
    write_array(Afield_k->get(),_param->nt, path);
    Afield_k->calc_pot();
    path = _param->path_results + "/Afield_k.dat";
    write_array(Afield_k->get(),_param->nt, path);
    path = _param->path_results + "/Bfield_k.dat";
    write_array(Bfield_k->get(),_param->nt, path);
}

void TDSESolver::_masks_CUSTOM(){
    double ib = _param->imax/5.0;
    double jb = _param->jmax/5.0;
    double kb = _param->kmax/5.0;
    double gamma = 0.95;
    for(int i=0; i<_param->ni; i++){
        if(_i[i] < _i[0]+ib){
            _imask[i] = pow(cos(M_PI*(_i[i]-(_i[0]+ib))*gamma/(2.0*ib)),1.0/8.0);
            if(_imask[i].real() < 0.0)
                _imask[i] = 0.0;
        }
        else if(_i[i]>(_i[_param->ni-1] - ib)){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param->ni - 1]-ib))*gamma/(2.0*ib)),1.0/8.0);
        }
        else{
            _imask[i] = 1.0;
        }

    }

    for(int j=0; j<_param->nj; j++){
        if(_j[j] < _j[0]+jb){
            _jmask[j] = pow(cos(M_PI*(_j[j]-(_j[0]+jb))*gamma/(2.0*jb)),1.0/8.0);
            if(_jmask[j].real() < 0.0)
                _jmask[j] = 0.0;
        }
        else if(_j[j]>(_j[_param->nj-1] - jb)){
            _jmask[j] = pow(cos(M_PI*(_j[j] - (_j[_param->nj - 1]-jb))*gamma/(2.0*jb)),1.0/8.0);
        }
        else{
            _jmask[j] = 1.0;
        }
    }

    for(int k=0; k<_param->nk; k++){
        if(_k[k] < _k[0]+kb){
            _kmask[k] = pow(cos(M_PI*(_k[k]-(_k[0]+kb))*gamma/(2.0*kb)),1.0/8.0);
            if(_kmask[k].real() < 0.0)
                _kmask[k] = 0.0;
        }
        else if(_k[k]>(_k[_param->nk-1] - kb)){
            _kmask[k] = pow(cos(M_PI*(_k[k] - (_k[_param->nk - 1]-kb))*gamma/(2.0*kb)),1.0/8.0);
        }
        else{
            _kmask[k] = 1.0;
        }
    }
}

void TDSESolver::_ipropagate_CUSTOM(){
    cdouble ener = 0.0;
    cdouble ener_old=0.0;
    double eps=10000.0;
    int counter =0;
    cdouble norm;
    cdouble proj=0.0;
    cdouble **psi_i_row, **psi_j_row, **psi_k_row;
    const int ni = _param->ni;
    const int nj = _param->nj;
    const int nk = _param->nk;
    psi_i_row = alloc2d<cdouble>(_param->n_threads,ni);
    psi_j_row = alloc2d<cdouble>(_param->n_threads,nj);
    psi_k_row = alloc2d<cdouble>(_param->n_threads,nk);

    for(int m=0; m<1; m++){
        ener =0.0;
        ener_old=0.0;
        eps = 1000;
        
        //_wf->gaussian_anti(0.0,0.0,0.0,0.1);
        
        cdouble norm = _wf->norm();
        (*_wf) /= norm;

        for(int n=0; n<_param->nt_ITP;n++){
        //while((eps > 1e-4)){ // || (counter<_param->nt_ITP)){
            double tstart, tend;
            tstart = omp_get_wtime();
            #pragma omp parallel for collapse(1) schedule(dynamic)
            for(int j=0;j<nj;j++){
                for(int k=0;k<nk;k++){
                    int id = omp_get_thread_num();
                    _wf->get_i_row(psi_i_row[id],j,k);
                    (_ham->*(_ham->step_i))(psi_i_row[id],j,k,0,1,id);
                    _wf->set_i_row(psi_i_row[id],j,k);
                }
            }

            #pragma omp parallel for collapse(1) schedule(dynamic)
            for(int i=0;i<ni;i++){
                for(int k=0;k<nk;k++){
                    int id = omp_get_thread_num();
                    _wf->get_j_row(psi_j_row[id],i,k);
                    (_ham->*(_ham->step_j))(psi_j_row[id],i,k,0,1,id);
                    _wf->set_j_row(psi_j_row[id],i,k);
                }
            }

            #pragma omp parallel for collapse(1) schedule(dynamic)
            for(int i=0;i<ni;i++){
                for(int j=0;j<nj;j++){
                    int id = omp_get_thread_num();
                    _wf->get_k_row(psi_k_row[id],i,j);
                    (_ham->*(_ham->step_k))(psi_k_row[id],i,j,0,1,id);
                    _wf->set_k_row(psi_k_row[id],i,j);
                }
            }
            //_wf->anti_sym_k();
            norm = _wf->norm();
            (*_wf) /= norm;
            tend = omp_get_wtime();
            std::cout<<"n: "<<counter<<" timestep: "<<tend-tstart<<"\n";

            if(counter%5==0){
                //_wf->grand_schmidt(); 
                ener = (_ham->*(_ham->ener))(_wf->get());
	            eps = std::abs((std::real(ener)-std::real(ener_old))/std::real(ener));
                std::cout<<"State: "<<m<<" Norm: "<< norm<<" Ener: "<<ener<<" Eps: "<<eps<<"\n";
                ener_old = ener;
            }
            counter++;

        }
        std::cout<<"Ener: "<<ener<<"\n";
        //_wf->set_to_eigen(m);
	//std::string name = "argon_"+ std::to_string(m);
	//_wf->save_wf2(name);
    }
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_j_row,_param->n_threads,nj);
    free2d(&psi_k_row,_param->n_threads,nk);
}

void TDSESolver::_propagate_CUSTOM(){
    cdouble ener = 0.0;
    cdouble ener_old=0.0;
    double eps=0.0;
    cdouble **psi_i_row, **psi_j_row, **psi_k_row;
    const int ni = _param->ni;
    const int nj = _param->nj;
    const int nk = _param->nk;
    cdouble norm;
    psi_i_row = alloc2d<cdouble>(_param->n_threads,ni);
    psi_j_row = alloc2d<cdouble>(_param->n_threads,nj);
    psi_k_row = alloc2d<cdouble>(_param->n_threads,nk);

    double tstart, tend;
    tstart = omp_get_wtime();
    ener = (_ham->*(_ham->ener))(_wf->get());
    std::cout<<"Initial energy:  "<<ener<<std::endl;
    for(int n=0; n<_param->nt;n++){
        //norm = _wf->norm();
        //std::cout<<"norm: "<<norm<<std::endl;
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int j=0;j<nj;j++){
            for(int k=0;k<nk;k++){
                int id = omp_get_thread_num();
                _wf->get_i_row(psi_i_row[id],j,k);
                (_ham->*(_ham->step_i))(psi_i_row[id],j,k,n,0,id);
                _wf->set_i_row_mask(psi_i_row[id],_imask,j,k);
            }
        }

        //ener = (_ham->*(_ham->ener))(_wf->get());
        //norm = _wf->norm();
        //std::cout<<"After x – norm: "<<norm<<" ener: "<<ener<<std::endl;
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int i=0;i<ni;i++){
            for(int k=0;k<nk;k++){
                int id = omp_get_thread_num();
                _wf->get_j_row(psi_j_row[id],i,k);
                (_ham->*(_ham->step_j))(psi_j_row[id],i,k,n,0,id);
                _wf->set_j_row_mask(psi_j_row[id],_jmask,i,k);
            }
        }

        //ener = (_ham->*(_ham->ener))(_wf->get());
        //norm = _wf->norm();
        //std::cout<<"After y – norm: "<<norm<<" ener: "<<ener<<std::endl;
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int i=0;i<ni;i++){
            for(int j=0;j<nj;j++){
                int id = omp_get_thread_num();
                _wf->get_k_row(psi_k_row[id],i,j);
                (_ham->*(_ham->step_k))(psi_k_row[id],i,j,n,0,id);
                _wf->set_k_row_mask(psi_k_row[id],_kmask,i,j);
            }
        }

        //ener = (_ham->*(_ham->ener))(_wf->get());
        //norm = _wf->norm();
        //std::cout<<"After z – norm: "<<norm<<" ener: "<<ener<<std::endl;
        //_diag->run_diagnostics(n);
        if(n%1==0){
            //_wf->grand_schmidt(); 
            norm = _wf->norm();
            ener = (_ham->*(_ham->ener))(_wf->get());
            eps = std::abs((std::real(ener)-std::real(ener_old))/std::real(ener));
            std::cout<<"n: "<<n<<" Norm: "<< norm<<" Ener: "<<ener<<" Eps: "<<eps<<std::endl;
            ener_old = ener;
        }
    }
    tend = omp_get_wtime();
    norm = _wf->norm();
    std::cout<<"Norm: "<<norm<<" Elapsed time: "<< tend-tstart<<std::endl;
    //_diag->write_diagnostics();
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_j_row,_param->n_threads,nj);
    free2d(&psi_k_row,_param->n_threads,nk);
}
