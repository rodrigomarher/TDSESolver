#ifndef FIELDS_H
#define FIELDS_H

#define SIN2 0
#define TRAP 1

#include <string>
class Field_TDSESolver{
    private:
        int _nt;
        bool _flag = false;
        double _dt, _tmax, _w;
        double *_t;
        double *_field;
        double env_sin2(double ti);
        double env_trap(double ti);
        
    public:
        Field_TDSESolver();
        Field_TDSESolver(double amp, double w, double phi, int env, double tmax, double *t, int nt);
        Field_TDSESolver(std::string &path, double tmax, double *t, const int nt);
        void calc_pot();
        double get(int i);
        double* get();
        double operator[](int i);
        ~Field_TDSESolver();
};

#endif
