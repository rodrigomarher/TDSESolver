#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "fields.h"
#include "debug.h"
#include "parameters.h"
#include "utils.h"

Field_TDSESolver::Field_TDSESolver(){
    _field = NULL;
}

Field_TDSESolver::Field_TDSESolver(double amp, double w, double phi, int env, double tmax,  double *t, int nt){
    _nt = nt; _tmax = tmax; _w = w; _t = t;
    _dt = _t[_nt-1]/(double)(_nt);
    _field = new double[_nt];

    for(int i=0; i<_nt; i++){
        if(env == 0){
            _field[i] = amp*env_sin2(_t[i])*sin(w*_t[i]+phi);
        }        
        else if(env == 1){
            _field[i] = amp*env_trap(_t[i])*sin(w*_t[i]+phi);
        }        
    }
    _flag = true;
}

Field_TDSESolver::Field_TDSESolver(std::string &path, double tmax, double *t, const int nt){
   const int n_file = calc_n_elem(path);
   if(n_file != nt){
        debug0("[Field::Field] File length not equal to number of temporal points.\n"); exit(1);
   }

   _nt  = nt; _tmax = tmax; _w = -1, _t = t;
   _dt = _t[_nt-1]/(double)(_nt);
   _field = new double[_nt];

   std::ifstream file;
   file.open(path);
   std::string line;
   int i = 0;
   if(file.is_open()){
        while(getline(file, line)){
            std::stringstream ss(line);
            ss >> _field[i]; 
            i++;
        }
        file.close();
   }
   else{debug0("[Field::Field] Unable to open file.\n"); exit(1);}
}

void Field_TDSESolver::calc_pot(){
    double *temp;
    temp = new double[_nt];
    
    for(int i=0; i<_nt;i++){
        temp[i] = 0.0;
        for (int j=0; j<=i; j++){
            temp[i] += _field[j];
        }
        temp[i] *= -C*_dt;
    }
    for(int i=0; i<_nt;i++){
        _field[i] = temp[i];
    }
    delete[] temp;
}

double Field_TDSESolver::operator[](int i){
    return _field[i];
}

double Field_TDSESolver::get(int i){
    return _field[i];
}

double* Field_TDSESolver::get(){
    return _field;
}
double Field_TDSESolver::env_sin2(double ti){
    if (ti<_tmax){
        return pow(sin(M_PI*ti/_tmax),2);
    }
    else {
        return 0.0;
    }

}

double Field_TDSESolver::env_trap(double ti){
    double T = 2*M_PI/_w;
    if (ti<_tmax+2.0*T){
         if (ti<T){
             return pow(sin(M_PI*ti/(2.0*T)),2);
         }
         else if (T<ti && ti<_tmax+T){
             return 1.0;
         }
         else{
             return pow(sin(M_PI*(ti-_tmax-T)/(2.0*T)+M_PI/2.0),2);
         }
    }
    else{
        return 0.0;
    }
}

Field_TDSESolver::~Field_TDSESolver(){
    delete[] _field;
}
