#include "fields.h"
#include "field_container.h"

FieldContainer::FieldContainer(const int n_fields, std::string defs){
    _n_fields = n_fields;

    char delim = ';';
    int counter = 0;
    std::strinstream ss(defs);
    std::string item;

    while(getline(ss,item,delim)){
        _def_vec.push_back(item);
        counter++;
    }
    if(counter != _n_fields){debug0("[Diagnostics] n_fields != number of definitions.");}
}

void FieldContainer::set_parameters(Parameters *param){
    _param = param;
    _nt = _param->nt;
}

void FieldContainer::create_fields(){
    for(int i=0; i<_n_fields;i++){
        Field temp_field = new Field(_def_vec[i]);    
    }     
} 
