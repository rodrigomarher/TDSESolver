#ifndef FIELD_CONTAINER_H
#define FIELD_CONTAINER_H

#include "fields.h"
#include "parameters.h"
#include <vectoreters>
#include <string>

class FieldContainer{
    private:
        int _n_fields;
        std::vector<std::string> _def_vec;
        std::vector<Field*> _efield_i_vec;
        std::vector<Field*> _efield_j_vec;
        std::vector<Field*> _efield_k_vec;
        std::vector<Field*> _bfield_i_vec;
        std::vector<Field*> _bfield_j_vec;
        std::vector<Field*> _bfield_k_vec;
        Parameters* _param;
        int _nt;
    public:
        Field Efield_i;
        Field Efield_j;
        Field Efield_k;
        Field Bfield_i;
        Field Bfield_j;
        Field Bfield_k;

        FieldContainer(const int n_fields, std::string defs);
        void set_parameters(Parameters *param);
        void create_fields();

};


#endif
