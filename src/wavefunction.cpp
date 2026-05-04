#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include "debug.h"
#include "utils.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters *param){
    _param = param;
}

void WF::set_geometry( double *i, double *j, double *k, const double di, const double dj, const double dk){    
    _ni = _param->ni;
    _nj = _param->nj;
    _nk = _param->nk;
    
    _wf = alloc3d<cdouble>(_ni, _nj, _nk);
    if(_param->geometry == XYZ || _param->geometry == CUSTOM)
        _wf_buf = alloc4d<cdouble>(_ni, _nj, _nk, 1);
    else
        _wf_buf = alloc4d<cdouble>(_ni, _nj, _nk, _param->nt_diag);

    _eigen_wf = alloc4d<cdouble>(_ni,_nj,_nk,1);

    _i_row = new cdouble[_ni];
    _j_row = new cdouble[_nj];
    _k_row = new cdouble[_nk];

    _diag_buf = new cdouble[_param->nt_diag];

    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf[i][j][k] = 0.0;
            }
        }
    }


    if(_param->geometry != XYZ && _param->geometry != CUSTOM){
        for(int n=0;n<2;n++){
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        //std::cout<<i<<" "<<k<<std::endl;
                        _eigen_wf[n][i][j][k] = 0.0;
                    }
                }
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
        case CUSTOM:
            _geom_CUSTOM();
            break;
    }
}

void WF::gaussian(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk; k++){
                _wf[i][j][k] = exp(-(_i[i]-i0)*(_i[i]-i0)/sigma- (_j[j]-j0)*(_j[j]-j0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
            }
        }
    }
}

void WF::gaussian_anti(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk; k++){
                _wf[i][j][k] = _k[k]*exp(-(_i[i]-i0)*(_i[i]-i0)/sigma- (_j[j]-j0)*(_j[j]-j0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
            }
        }
    }
}

void WF::exponential(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] = exp(-sqrt((_i[i]-i0)*(_i[i]-i0) + (_j[j]-j0)*(_j[j]-j0) + (_k[k]-k0)*(_k[k]-k0)));
            }
        }
    }
}

void WF::cubic_spline(std::string path){

    std::cout<<"Starting cubic spline data loading..."<<std::endl; 
    int nr_interp = calc_n_elem(path);
    std::vector<double> rho_interp(nr_interp);
    std::vector<std::vector<double>> data_cs(4, std::vector<double>(nr_interp));
    
    std::ifstream file;
    std::string line;
    int idx = 0;
    file.open(path);
    if(file.is_open()){
        while(getline(file, line)){
            std::stringstream s(line);
            double rho, c0, c1, c2, c3;
            std::string tmp;
            getline(s, tmp, ' ');
            rho = stod(tmp);
            getline(s, tmp, ' ');
            c0 = stod(tmp);
            getline(s, tmp, ' ');
            c1 = stod(tmp);
            getline(s, tmp, ' ');
            c2 = stod(tmp);
            getline(s, tmp, ' ');
            c3 = stod(tmp);

            rho_interp[idx] = rho;
            data_cs[0][idx] = c0;
            data_cs[1][idx] = c1;
            data_cs[2][idx] = c2;
            data_cs[3][idx] = c3;
            idx++;
        }
    }
    std::cout<<"Cubic spline data loaded. Number of points: "<<idx<<std::endl; 
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk; k++){
                 
                int idx_rho_node;
                std::vector<double>::iterator iterator_pos;
                double rho = sqrt(_i[i]*_i[i] + _j[j]*_j[j] + _k[k]*_k[k]);
                double theta = atan2(sqrt(_i[i]*_i[i] + _j[j]*_j[j]),_k[k]);
                double phi = atan2(_j[j],_i[i]);
                double interpolation_result;
                cdouble radial_component;
                cdouble angular_component;

                rho = rho>rho_interp[nr_interp-1] ? rho_interp[nr_interp-1] : rho;
                rho = rho<rho_interp[0] ? rho_interp[0] : rho;
                iterator_pos = std::lower_bound(rho_interp.begin(), rho_interp.end(), rho);
                idx_rho_node = iterator_pos - rho_interp.begin() - 1; 
                
                double k1 = rho-rho_interp[idx_rho_node];
                double k2 = pow(rho-rho_interp[idx_rho_node], 2);
                double k3 = pow(rho-rho_interp[idx_rho_node], 3);
                interpolation_result = data_cs[3][idx_rho_node] + 
                                       data_cs[2][idx_rho_node]*k1 + 
                                       data_cs[1][idx_rho_node]*k2 + 
                                       data_cs[0][idx_rho_node]*k3;
                radial_component = cdouble(interpolation_result, 0.0); 
                angular_component = 0.5*(_spherical_harmonics(theta, phi, 1, -1) - _spherical_harmonics(theta,phi,1,1));

                _wf[i][j][k] = radial_component*angular_component;
            }
        }
    }

}

cdouble  WF::_spherical_harmonics(double theta, double phi, int l, int m){
    switch(l){
        case 0:
            return cdouble(0.5*sqrt(1.0/M_PI),0.0);
            break;

        case 1:
            switch(m){
                case -1:
                    return 0.5*sqrt(3.0/(2.0*M_PI))*sin(theta)*cdouble(cos(phi), -sin(phi));
                    break;
                case 0:
                    return 0.5*sqrt(3.0/(M_PI))*cos(theta)*cdouble(1.0, 0.0);
                    break;
                case 1:
                    return -0.5*sqrt(3.0/(2.0*M_PI))*sin(theta)*cdouble(cos(phi), sin(phi));
                    break;
                default:
                    debug0("[WF::_spherical_harmonics] m value out of bounds (|m|>l)\n");
                    exit(1);
                    break;
                }
         break;

         default:
            debug0("[WF::_spherical_harmonics] l value not implemented\n");
            exit(1);
            break;
    }
}

cdouble WF::norm(){
    cdouble integral = 0.0;
    double sum = 0.0;
    double dr = 0.0;
    switch(_param->geometry){
        case X:
            for(int i=0; i<_ni;i++){
                integral += _wf[i][0][0]*conj(_wf[i][0][0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf[i][0][k]*conj(_wf[i][0][k])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += 2*M_PI*_i[i]*_wf[i][0][k]*conj(_wf[i][0][k])*_di*_dk;
                }
            }
            break;
        case XYZ:
            //#pragma omp parallel for reduction(+:sum)
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        sum += (_wf[i][j][k]*conj(_wf[i][j][k])).real()*_di*_dk*_dj;
                    }
                }
            }
            integral = sum;
            break;
        case CUSTOM:
            for (int i=0; i<_ni-1; i++){
                for (int j=0; j<_nj-1; j++){
                    for (int k=0; k<_nk-1; k++){
                        dr = (_i[i+1] - _i[i])*(_j[j+1] - _j[j])*(_k[k+1] - _k[k]);
                        sum += (_wf[i][j][k]*conj(_wf[i][j][k])).real()*dr;
                    }
                }
            } 
            integral = sum;
            break;  
    }
    return sqrt(integral);
}

cdouble WF::norm_buf(int idx){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case X:{
            for(int i=0; i<_ni;i++){
                integral += _wf_buf[idx][i][0][0]*conj(_wf_buf[idx][i][0][0])*_di;
            }
            }
            break;
        case XZ:{
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf_buf[idx][i][0][k]*conj(_wf_buf[idx][i][0][k])*_di*_dk;
                }
            }
            }
            break;
        case RZ:{
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += 2*M_PI*_i[i]*_wf_buf[idx][i][0][k]*conj(_wf_buf[idx][i][0][k])*_di*_dk;
                }
            }
            }
            break;
        case XYZ:{
            double sum = 0.0;
            //#pragma omp parallel for reduction(+:sum)
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        sum += (_wf_buf[idx][i][j][k]*conj(_wf_buf[idx][i][j][k])).real()*_di*_dk*_dj;
                    }
                }
            }
            integral = sum;
            }
            break;
        case CUSTOM:{
            double sum = 0.0;
            double dr = 0.0;
            for (int i=0; i<_ni-1; i++){
                for (int j=0; j<_nj-1; j++){
                    for (int k=0; k<_nk-1; k++){
                        dr = (_i[i+1] - _i[i])*(_j[j+1] - _j[j])*(_k[k+1] - _k[k]);
                        sum += (_wf_buf[idx][i][j][k]*conj(_wf_buf[idx][i][j][k])).real()*dr;
                    }
                }
            } 
            integral = sum;
            }
            break;
    }
    return sqrt(integral);
}



void WF::apply_mask(cdouble *imask, cdouble *jmask, cdouble *kmask){
    (this->*(this->_apply_mask))(imask, jmask, kmask);
}


cdouble*** WF::get(){
    return _wf;
}

cdouble* WF::i_row(int j, int k){
    for(int i=0; i<_ni; i++){
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

void WF::anti_sym_k(){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0;k<_nk/2;k++){
                _wf[i][j][k] = (_wf[i][j][k] - _wf[i][j][_nk-k-1])/2.0;           
                _wf[i][j][_nk-k-1] = -_wf[i][j][k];
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

void WF::set_to_eigen(int n){
    for(int i=0; i<_ni; i++ ){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk; k++){
                _eigen_wf[n][i][j][k] = _wf[i][j][k];
            }
        }
    }
}

cdouble*** WF::get_eigen(int n){
    return _eigen_wf[n];
}


void WF::set_to_ground(){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _eigen_wf[0][i][j][k] = _wf[i][j][k];
            }
        }
    }
}

cdouble*** WF::get_ground(){
    return _eigen_wf[0];
}

cdouble WF::operator()(int i, int j, int k){
    return _wf[i][j][k];
}

cdouble WF::project(cdouble ***phi){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case X:
            for(int i=0; i<_ni;i++){
                integral += _wf[i][0][0]*conj(phi[i][0][0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf[i][0][k]*conj(phi[i][0][k])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += 2*M_PI*_i[i]*_wf[i][0][k]*conj(phi[i][0][k])*_di*_dk;
                }
            }
            break;
        case XYZ:
            //#pragma omp parallel for reduction(+:sum)
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        integral += (_wf[i][j][k]*conj(phi[i][j][k]))*_di*_dk*_dj;
                    }
                }
            }
            break;
    }
    return integral;
}

void WF::grand_schmidt(){
    std::vector<cdouble> proj;
    for(int i=0;i<5;i++){
        proj.push_back(project(_eigen_wf[i]));
    }


    std::cout<<"Projection: "<<proj[0]<<" "<<proj[1]<<" "<<proj[2]<<" "<<proj[3]<<" "<<proj[4]<<std::endl;
    for (int i=0; i<_ni; i++){
        for (int j=0; j<_nj; j++){
            for (int k=0; k<_nk; k++){
                for(int n=0; n<5; n++){
                    _wf[i][j][k] -= proj[n]*_eigen_wf[n][i][j][k];
                } 
            }
        }
    }
    proj.clear();
}


void WF::operator/=(cdouble val){
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0; i<_ni;i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] /= val;
            }
        }
    }
}

void WF::save_wf2(std::string name){

    switch(_param->geometry){
        case XYZ: case CUSTOM:
            cdouble value;
            std::ofstream outfile_X;
            std::ofstream outfile_Y;
            std::ofstream outfile_Z;
            outfile_X.open(_param->path_results + "/" + name + "_X.dat");
            outfile_Y.open(_param->path_results + "/" + name + "_Y.dat");
            outfile_Z.open(_param->path_results + "/" + name + "_Z.dat");

            for(int j=0;j<_nj;j++){
                for(int k=0; k<_nk;k++){
                    value = _wf[_ni/2][j][k]*conj(_wf[_ni/2][j][k]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_X<<doubleStr.str()<<std::endl;
                }
            }

            for(int i=0;i<_ni;i++){
                for(int k=0; k<_nk;k++){
                    value = _wf[i][_nj/2][k]*conj(_wf[i][_nj/2][k]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_Y<<doubleStr.str()<<std::endl;
                }
            }

            for(int i=0;i<_ni;i++){
                for(int j=0; j<_nj;j++){
                    value = _wf[i][j][_nk/2]*conj(_wf[i][j][_nk/2]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_Z<<doubleStr.str()<<std::endl;
                }
            }
            outfile_X.close();
            outfile_Y.close();
            outfile_Z.close();                 
            break;
    }
}

WF::~WF(){
    free3d(&_wf,_ni,_nj,_nk);
    if(_param->geometry == XYZ)
        free4d(&_wf_buf,_ni,_nj,_nk, 1);
    else
        free4d(&_wf_buf, _ni, _nj, _nk, _param->nt_diag);
    free4d(&_eigen_wf, _ni, _nj, _nk, 5);
    //if(_param->geometry != XYZ)
    //    free3d(&_wf_0,_ni,_nk,_nk);
    //else
    //    free3d(&_wf_0,1,1,1);
    delete[] _diag_buf;
    delete[] _i_row;
    delete[] _j_row;
    delete[] _k_row;
}
