#include <iostream>
#include "hamiltonian.h"
#include <omp.h>
void Hamiltonian::_allocate_CUSTOM(){
    step_i = &Hamiltonian::step_i_CUSTOM;
    step_k = &Hamiltonian::step_k_CUSTOM;
    step_j = &Hamiltonian::step_j_CUSTOM;
    ener = &Hamiltonian::ener_CUSTOM;

    _Mk_du = new cdouble[_param->n_threads*_nk];
    _Mk_d  = new cdouble[_param->n_threads*_nk];
    _Mk_dl = new cdouble[_param->n_threads*_nk];
    _Mpk_du = new cdouble[_param->n_threads*_nk];
    _Mpk_d  = new cdouble[_param->n_threads*_nk];
    _Mpk_dl = new cdouble[_param->n_threads*_nk];
    _lhs_k = new cdouble[_param->n_threads*_nk];
    _res_k = new cdouble[_param->n_threads*_nk];

    _Mj_du = new cdouble[_param->n_threads*_nj];
    _Mj_d  = new cdouble[_param->n_threads*_nj];
    _Mj_dl = new cdouble[_param->n_threads*_nj];
    _Mpj_du = new cdouble[_param->n_threads*_nj];
    _Mpj_d  = new cdouble[_param->n_threads*_nj];
    _Mpj_dl = new cdouble[_param->n_threads*_nj];
    _lhs_j = new cdouble[_param->n_threads*_nj];
    _res_j = new cdouble[_param->n_threads*_nj];

    _Mi_du = new cdouble[_param->n_threads*_ni];
    _Mi_d  = new cdouble[_param->n_threads*_ni];
    _Mi_dl = new cdouble[_param->n_threads*_ni];
    _Mpi_du = new cdouble[_param->n_threads*_ni];
    _Mpi_d  = new cdouble[_param->n_threads*_ni];
    _Mpi_dl = new cdouble[_param->n_threads*_ni];
    _lhs_i = new cdouble[_param->n_threads*_ni];
    _res_i = new cdouble[_param->n_threads*_ni];
}

void Hamiltonian::step_i_CUSTOM(cdouble *psi_i_row, const int j, const int k, const int ti, const int imag, const int id_thread){
    cdouble Hx_du;
    cdouble Hx_d;
    cdouble Hx_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;

    cdouble afield_i = (*Afield_i)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a, b, c;
    for(int i=1;i<_ni-1;i++){
        a =  1.0/((_i[i+1]-_i[i])*(_i[i]-_i[i-1])) + 0.0*0.5/(C*C)*afield_i*afield_i;
        b =  -1.0/((_i[i+1]-_i[i-1])*(_i[i+1]-_i[i])) + 0.0*I/(C*(_i[i+1]-_i[i-1]))*afield_i;
        c =  -1.0/((_i[i+1]-_i[i-1])*(_i[i]-_i[i-1])) - 0.0*I/(C*(_i[i+1]-_i[i-1]))*afield_i;
        Hx_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k], 0);
        Hx_du =  b;
        Hx_dl =  c;

        _Mi_du[id_thread*_ni + i] = I*Hx_du*dt/2.0;
        _Mi_d[id_thread*_ni + i]  = 1.0 + I*Hx_d*dt/2.0;
        _Mi_dl[id_thread*_ni + i] = I*Hx_dl*dt/2.0;
        _Mpi_du[id_thread*_ni + i] = -I*Hx_du*dt/2.0;
        _Mpi_d[id_thread*_ni + i]  = 1.0 - I*Hx_d*dt/2.0;
        _Mpi_dl[id_thread*_ni + i] = -I*Hx_dl*dt/2.0;
    }
    _Mi_du[id_thread*_ni] = _Mi_du[id_thread*_ni + 1];
    _Mi_d[id_thread*_ni] = _Mi_d[id_thread*_ni + 1];
    _Mi_dl[id_thread*_ni] = _Mi_dl[id_thread*_ni + 1];
    _Mpi_du[id_thread*_ni] = _Mpi_du[id_thread*_ni + 1];
    _Mpi_d[id_thread*_ni] = _Mpi_d[id_thread*_ni + 1];
    _Mpi_dl[id_thread*_ni] = _Mpi_dl[id_thread*_ni + 1];

    _Mi_du[id_thread*_ni + _ni-1] = _Mi_du[id_thread*_ni + _ni-2];
    _Mi_d[id_thread*_ni + _ni-1] = _Mi_d[id_thread*_ni + _ni-2];
    _Mi_dl[id_thread*_ni + _ni-1] = _Mi_dl[id_thread*_ni + _ni-2];
    _Mpi_du[id_thread*_ni + _ni-1] = _Mpi_du[id_thread*_ni + _ni-2];
    _Mpi_d[id_thread*_ni + _ni-1] = _Mpi_d[id_thread*_ni + _ni-2];
    _Mpi_dl[id_thread*_ni + _ni-1] = _Mpi_dl[id_thread*_ni + _ni-2];

    tridot(&_Mpi_dl[id_thread*_ni],
           &_Mpi_d[id_thread*_ni],
           &_Mpi_du[id_thread*_ni],
           psi_i_row,&_lhs_i[id_thread*_ni],_ni); 
    tdma(&_Mi_dl[id_thread*_ni],
         &_Mi_d[id_thread*_ni],
         &_Mi_du[id_thread*_ni],
         &_lhs_i[id_thread*_ni],&_res_i[id_thread*_ni],_ni);
    for(int i=0; i<_ni; i++){
        psi_i_row[i] = _res_i[id_thread*_ni + i];
    }
}

void Hamiltonian::step_j_CUSTOM(cdouble *psi_j_row, const int i, const int k, const int ti, const int imag, const int id_thread){
    cdouble Hy_du;
    cdouble Hy_d;
    cdouble Hy_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;
    
    cdouble afield_j = (*Afield_j)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a,b,c;
    for(int j=1;j<_nj-1;j++){
        a =  1.0/((_j[j+1]-_j[j])*(_j[j]-_j[j-1])) + 0.0*0.5/(C*C)*afield_j*afield_j;
        b = -1.0/((_j[j+1]-_j[j-1])*(_j[j+1]-_j[j])) + 0.0*I/(C*(_j[j+1]-_j[j-1]))*afield_j;
        c = -1.0/((_j[j+1]-_j[j-1])*(_j[j]-_j[j-1])) - 0.0*I/(C*(_j[j+1]-_j[j-1]))*afield_j;
        Hy_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
        Hy_du =  b;
        Hy_dl =  c;

        _Mj_du[id_thread*_nj + j] = I*Hy_du*dt/2.0;
        _Mj_d[id_thread*_nj + j]  = 1.0 + I*Hy_d*dt/2.0;
        _Mj_dl[id_thread*_nj + j] = I*Hy_dl*dt/2.0;
        _Mpj_du[id_thread*_nj + j] = -I*Hy_du*dt/2.0;
        _Mpj_d[id_thread*_nj + j]  = 1.0 - I*Hy_d*dt/2.0;
        _Mpj_dl[id_thread*_nj + j] = -I*Hy_dl*dt/2.0;
    }
    _Mj_du[id_thread*_nj] = _Mj_du[id_thread*_nj + 1];
    _Mj_d[id_thread*_nj] = _Mj_d[id_thread*_nj + 1];
    _Mj_dl[id_thread*_nj] = _Mj_dl[id_thread*_nj + 1];
    _Mpj_du[id_thread*_nj] = _Mpj_du[id_thread*_nj + 1];
    _Mpj_d[id_thread*_nj] = _Mpj_d[id_thread*_nj + 1];
    _Mpj_dl[id_thread*_nj] = _Mpj_dl[id_thread*_nj + 1];

    _Mj_du[id_thread*_nj + _nj-1] = _Mj_du[id_thread*_nj + _nj-2];
    _Mj_d[id_thread*_nj + _nj-1] = _Mj_d[id_thread*_nj + _nj-2];
    _Mj_dl[id_thread*_nj + _nj-1] = _Mj_dl[id_thread*_nj + _nj-2];
    _Mpj_du[id_thread*_nj + _nj-1] = _Mpj_du[id_thread*_nj + _nj-2];
    _Mpj_d[id_thread*_nj + _nj-1] = _Mpj_d[id_thread*_nj + _nj-2];
    _Mpj_dl[id_thread*_nj + _nj-1] = _Mpj_dl[id_thread*_nj + _nj-2];
    tridot(&_Mpj_dl[id_thread*_nj],
           &_Mpj_d[id_thread*_nj],
           &_Mpj_du[id_thread*_nj],
           psi_j_row, &_lhs_j[id_thread*_nj], _nj); 
    tdma(&_Mj_dl[id_thread*_nj],
         &_Mj_d[id_thread*_nj],
         &_Mj_du[id_thread*_nj],
         &_lhs_j[id_thread*_nj],&_res_j[id_thread*_nj],_nj);
    for(int j=0; j<_nj; j++){
        psi_j_row[j] = _res_j[id_thread*_nj + j];
    }
}
void Hamiltonian::step_k_CUSTOM(cdouble *psi_k_row, const int i, const int j, const int ti, const int imag, const int id_thread){
    cdouble Hz_du;
    cdouble Hz_d;
    cdouble Hz_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;
    
    cdouble afield_k = (*Afield_k)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a,b,c;
    for(int k=1;k<_nk-1;k++){
        a =  1.0/((_k[k+1]-_k[k])*(_k[k]-_k[k-1])) + 0.0*0.5/(C*C)*afield_k*afield_k;
        b = -1.0/((_k[k+1]-_k[k-1])*(_k[k+1]-_k[k])) + 0.0*I/(C*(_k[k+1]-_k[k-1]))*afield_k;
        c = -1.0/((_k[k+1]-_k[k-1])*(_k[k]-_k[k-1])) - 0.0*I/(C*(_k[k+1]-_k[k-1]))*afield_k;
        Hz_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
        Hz_du =  b;
        Hz_dl =  c;

        _Mk_du[id_thread*_nk + k] = I*Hz_du*dt/2.0;
        _Mk_d[id_thread*_nk + k]  = 1.0 + I*Hz_d*dt/2.0;
        _Mk_dl[id_thread*_nk + k] = I*Hz_dl*dt/2.0;
        _Mpk_du[id_thread*_nk + k] = -I*Hz_du*dt/2.0;
        _Mpk_d[id_thread*_nk + k]  = 1.0 - I*Hz_d*dt/2.0;
        _Mpk_dl[id_thread*_nk + k] = -I*Hz_dl*dt/2.0;
    }
    _Mk_du[id_thread*_nk] = _Mk_du[id_thread*_nk + 1];
    _Mk_d[id_thread*_nk] = _Mk_d[id_thread*_nk + 1];
    _Mk_dl[id_thread*_nk] = _Mk_dl[id_thread*_nk + 1];
    _Mpk_du[id_thread*_nk] = _Mpk_du[id_thread*_nk + 1];
    _Mpk_d[id_thread*_nk] = _Mpk_d[id_thread*_nk + 1];
    _Mpk_dl[id_thread*_nk] = _Mpk_dl[id_thread*_nk + 1];

    _Mk_du[id_thread*_nk + _nk-1] = _Mk_du[id_thread*_nk + _nk-2];
    _Mk_d[id_thread*_nk + _nk-1] = _Mk_d[id_thread*_nk + _nk-2];
    _Mk_dl[id_thread*_nk + _nk-1] = _Mk_dl[id_thread*_nk + _nk-2];
    _Mpk_du[id_thread*_nk + _nk-1] = _Mpk_du[id_thread*_nk + _nk-2];
    _Mpk_d[id_thread*_nk + _nk-1] = _Mpk_d[id_thread*_nk + _nk-2];
    _Mpk_dl[id_thread*_nk + _nk-1] = _Mpk_dl[id_thread*_nk + _nk-2];
    tridot(&_Mpk_dl[id_thread*_nk],
           &_Mpk_d[id_thread*_nk],
           &_Mpk_du[id_thread*_nk],
           psi_k_row, &_lhs_k[id_thread*_nk], _nk); 
    tdma(&_Mk_dl[id_thread*_nk],
         &_Mk_d[id_thread*_nk],
         &_Mk_du[id_thread*_nk],
         &_lhs_k[id_thread*_nk],&_res_k[id_thread*_nk],_nk);
    for(int k=0; k<_nk; k++){
        psi_k_row[k] = _res_k[id_thread*_nk + k];
    }
}

cdouble Hamiltonian::ener_CUSTOM(cdouble ***psi){
    cdouble integral=0.0;

    cdouble *temp_x, *temp_y, *temp_z;
    cdouble *i_row, *j_row, *k_row;
    cdouble *Hx_du, *Hx_dl, *Hx_d;
    cdouble *Hz_du, *Hz_dl, *Hz_d; 
    cdouble *Hy_du, *Hy_dl, *Hy_d;

    double integral_z, integral_y, integral_x;
    // Apply Hz to wavefunction
    temp_z = new cdouble[_param->n_threads*_nk];
    k_row = new cdouble[_param->n_threads*_nk];
    Hz_dl = new cdouble[_param->n_threads*_nk];
    Hz_d = new cdouble[_param->n_threads*_nk];
    Hz_du = new cdouble[_param->n_threads*_nk];
    integral_z = 0.0;
    {    
    #pragma omp parallel for schedule(dynamic) collapse(1)  reduction(+: integral_z)
    for(int i=1;i<_ni-1;i++){
        for(int j=1;j<_nj-1;j++){
            int id = omp_get_thread_num();
            for(int k=1;k<_nk-1;k++){
                cdouble a,b,c;
                a =  1.0/((_k[k+1]-_k[k])*(_k[k]-_k[k-1]));
                b = -1.0/((_k[k+1]-_k[k-1])*(_k[k+1]-_k[k]));
                c = -1.0/((_k[k+1]-_k[k-1])*(_k[k]-_k[k-1]));
                Hz_du[id*_nk + k] = b;
                Hz_d[id*_nk + k]  = a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hz_dl[id*_nk + k] = c;
                k_row[id*_nk + k] = psi[i][j][k];
            }
            Hz_du[id*_nk] = Hz_du[id*_nk+1];
            Hz_d[id*_nk] = Hz_d[id*_nk+1];
            Hz_dl[id*_nk] = Hz_dl[id*_nk+1];
            k_row[id*_nk] = psi[i][j][0];

            Hz_du[id*_nk+_nk-1] = Hz_du[id*_nk+_nk-2];
            Hz_d[id*_nk+_nk-1] = Hz_d[id*_nk+_nk-2];
            Hz_dl[id*_nk+_nk-1] = Hz_dl[id*_nk+_nk-2];
            k_row[id*_nk+_nk-1] = psi[i][j][_nk-1];

            tridot(&Hz_dl[id*_nk], &Hz_d[id*_nk], &Hz_du[id*_nk], &k_row[id*_nk], &temp_z[id*_nk], _nk);
            for(int  k=1;k<_nk-1;k++){
                double dr =(_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                integral_z += (conj(psi[i][j][k])*temp_z[id*_nk + k]).real()*dr;
            }
        }
    }
    }
    delete[] temp_z;
    delete[] k_row;
    delete[] Hz_du;
    delete[] Hz_d;
    delete[] Hz_dl; 

    // Apply Hx to wavefunction
    temp_x = new cdouble[_param->n_threads*_ni];
    i_row = new cdouble[_param->n_threads*_ni];
    Hx_dl = new cdouble[_param->n_threads*_ni];
    Hx_d = new cdouble[_param->n_threads*_ni];
    Hx_du = new cdouble[_param->n_threads*_ni];
    integral_x = 0.0;
    {
    #pragma omp parallel for schedule(dynamic) collapse(1)  reduction(+: integral_x)
    for(int j=1;j<_nj-1;j++){
        for(int k=1;k<_nk-1;k++){
            int id = omp_get_thread_num();
            for(int i=1;i<_ni-1;i++){
                cdouble a,b,c;
                a =  1/((_i[i+1]-_i[i])*(_i[i]-_i[i-1]));
                b = -1/((_i[i+1]-_i[i-1])*(_i[i+1]-_i[i]));
                c = -1/((_i[i+1]-_i[i-1])*(_i[i]-_i[i-1]));

                Hx_du[id*_ni + i] = b;
                Hx_d[id*_ni + i]  = a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hx_dl[id*_ni + i] = c;
                i_row[id*_ni + i] = psi[i][j][k];
            }
            Hx_du[id*_ni] = Hx_du[id*_ni+1];
            Hx_d[id*_ni] = Hx_d[id*_ni+1];
            Hx_dl[id*_ni] = Hx_dl[id*_ni+1];
            i_row[id*_ni] = psi[0][j][k];

            Hx_du[id*_ni+_ni-1] = Hx_du[id*_ni+_ni-2];
            Hx_d[id*_ni+_ni-1] = Hx_d[id*_ni+_ni-2];
            Hx_dl[id*_ni+_ni-1] = Hx_dl[id*_ni+_ni-2];
            i_row[id*_ni+_ni-1] = psi[_ni-1][j][k];

            tridot(&Hx_dl[id*_ni], &Hx_d[id*_ni], &Hx_du[id*_ni], &i_row[id*_ni], &temp_x[id*_ni], _ni);
            for(int i=1; i<_ni-1; i++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                integral_x += (conj(psi[i][j][k])*temp_x[id*_ni + i]).real()*dr;
            }
        }
    }
    }
    delete[] temp_x;
    delete[] i_row;
    delete[] Hx_du;
    delete[] Hx_d;
    delete[] Hx_dl;

    // Apply Hy to wavefunction
    temp_y = new cdouble[_param->n_threads*_nj];
    j_row = new cdouble[_param->n_threads*_nj];
    Hy_dl = new cdouble[_param->n_threads*_nj];
    Hy_d = new cdouble[_param->n_threads*_nj];
    Hy_du = new cdouble[_param->n_threads*_nj];
    integral_y = 0.0;
    {
    #pragma omp parallel for schedule(dynamic) collapse(1) reduction(+: integral_y)
    for(int i=1;i<_ni-1;i++){
        for(int k=1;k<_nk-1;k++){
            cdouble a,b,c;
            int id = omp_get_thread_num();
            for(int j=1;j<_nj-1;j++){
                a =  1.0/((_j[j+1]-_j[j])*(_j[j]-_j[j-1]));
                b = -1.0/((_j[j+1]-_j[j-1])*(_j[j+1]-_j[j]));
                c = -1.0/((_j[j+1]-_j[j-1])*(_j[j]-_j[j-1]));

                Hy_du[id*_nj + j] = b;
                Hy_d[id*_nj + j]  = a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hy_dl[id*_nj + j] = c;
                j_row[id*_nj + j] = psi[i][j][k];
            }
            Hy_du[id*_nj] = Hy_du[id*_nj+1];
            Hy_d[id*_nj] = Hy_d[id*_nj+1];
            Hy_dl[id*_nj] = Hy_dl[id*_nj+1];
            j_row[id*_nj] = psi[i][0][k];

            Hy_du[id*_nj+_nj-1] = Hy_du[id*_nj+_nj-2];
            Hy_d[id*_nj+_nj-1] = Hy_d[id*_nj+_nj-2];
            Hy_dl[id*_nj+_nj-1] = Hy_dl[id*_nj+_nj-2];
            j_row[id*_nj+_nj-1] = psi[i][_nj-1][k];

            tridot(&Hy_dl[id*_nj], &Hy_d[id*_nj], &Hy_du[id*_nj], &j_row[id*_nj], &temp_y[id*_nj], _nj);
            for(int j=1; j<_nj-1; j++){
                double dr = (_i[i+1]-_i[i])*(_j[j+1]-_j[j])*(_k[k+1]-_k[k]);
                integral_y += (conj(psi[i][j][k])*temp_y[id*_nj + j]).real()*dr;
            }
        }
    }
    }
    delete[] temp_y;
    delete[] j_row;
    delete[] Hy_dl;
    delete[] Hy_d;
    delete[] Hy_du;
    // Integrate 
    integral = integral_x + integral_y + integral_z;
    return integral; 
}

