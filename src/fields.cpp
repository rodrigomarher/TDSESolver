#include <iostream>
#include <cmath>
#include <vector>
#include "fields.h"
#include "debug.h"
#include "parameters.h"

Field::Field(){
    _field = NULL;
}

Field::Field(std::string def, double *t, int nt){
    _def = def;
    _t =t;
    _nt = nt;
    _dt = _t[_nt-1]/(double)(_nt);
    _parse_def();
    _init_field();
}

Field::Field(double amp, double w, double phi, t_envelope env, double tmax, double t0,  double *t, int nt){
    _nt = nt; _tmax = tmax; _w = w; _t = t;
    _dt = _t[_nt-1]/(double)(_nt);
    _amp = amp; _phi = phi; _env = env; _t0 = t0; 

    _init_field();
}

void Field::_init_field(){
    _field = new double[_nt];

    for(int i=0; i<_nt; i++){
        if(_env == SIN2){
            _field[i] = _amp*env_sin2(_t[i],_t0)*sin(_w*(_t[i]-_t0)+_phi);
        }        
        else if(_env == TRAP){
            _field[i] = _amp*env_trap(_t[i],_t0)*sin(_w*(_t[i]-_t0)+_phi);
        }        
        else if(_env == CTE){
            _field[i] = _amp*env_cte(_t[i],_t0)*sin(_w*(_t[i]-_t0)+_phi);
        }        
    }
}

void Field::_parse_def(){
    char delim = ',';
    std::stringstream ss(_def);
    std::string item;
    std::vector<std:string> parsed;
    while(getline(ss,item,delim)){
        parsef.push_back(item);
    }
    
    if(parsed[0] == "efield_i"){
        _tfield = EFIELD_I; 
    }
    else if (parsed[0] == "efield_j"){
        _tfield = EFIELD_J;
    }
    else if (parsed[0] == "efield_k"){
        _tfield = EFIELD_K;
    }
    else if(parsed[0] == "bfield_i"){
        _tfield = BFIELD_I;
    }
    else if(parsed[0] == "bield_j"){
        _tfield = BFIELD_J;
    }
    else if(parsed[0] == "bfield_k"){
        _tfield = BFIELD_K;
    }
    
    if (parsed[1] == "sin2"){
        _env = SIN2;
    }
    else if (parsed[1] == "trap"){
        _env = TRAP;
    }
    else{
        _env = CTE,
    }
    

    _amp  = parsed[2];
    _w    = parsed[3];
    _t0   = parsed[4];
    _tmax = parsed[5];
    _phi  = parsed[6];
}

void Field::calc_pot(){
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

double Field::operator[](int i){
    return _field[i];
}

double Field::get(int i){
    return _field[i];
}

double* Field::get(){
    return _field;
}
double Field::env_sin2(double ti, double t0){
    if (abs(ti-t0)<_tmax){
        return pow(sin(M_PI*(ti-t0)/_tmax + M_PI/2),2);
    }
    else {
        return 0.0;
    }

}

double Field::env_trap(double ti,double t0){
    double T = 2*M_PI/_w;
    if (ti>(t0-_tmax/2-T) && abs(ti-t0+_tmax/2+T)<(_tmax+2.0*T)){
         if ((ti-t0+_tmax/2.0 + T)<T){
             return pow(sin(M_PI*(ti-t0+T+_tmax/2.0)/(2.0*T)),2);
         }
         else if (T<(ti-t0+_tmax/2.0+T) && (ti-t0+_tmax/2.0+T)<_(tmax+T)){
             return 1.0;
         }
         else{
             return pow(sin(M_PI*(ti-t0-_tmax/2.0-T)/(2.0*T)),2);
         }
    }
    else{
        return 0.0;
    }
}

double Field::env_cte(double ti, double t0){
    return 1.0;
    }

std::string Field::type(){
    return _t_field;
}

Field::~Field(){
    delete[] _field;
}
