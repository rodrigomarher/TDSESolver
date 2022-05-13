#ifndef FIELDS_H
#define FIELDS_H

#define SIN2 0
#define TRAP 1

#include <string>

enum t_envelope { SIN2,
                  TRAP,
                  CTE,                 
		        };
enum t_field { EFIELD_I,
               EFIELD_J,
               EFIELD_K,
               BFIELD_I,
               BFIELD_J,
               BFIELD_K
             }

class Field{
    private:
        int _nt;
        bool _flag = false;
        double _dt, _tmax, _amp, _w, _t0, _phi;
        t_envelope _env;
        t_field _tfield;
        double *_t;
        double *_field;
        std::string _def;
        double env_sin2(double ti,t0);
        double env_trap(double ti,t0);
        double env_cte(double ti,t0);
        void _parse_def();
        void _init_field();
        
    public:
        Field();
        Field(std::string def, double *t, int nt);
        Field(double amp, double w, double phi, int env, double tmax, double t0, double *t, int nt);
        void calc_pot();
        double get(int i);
        double* get();
        double operator[](int i);
        void save();
        int type();
        ~Field();
};

#endif
