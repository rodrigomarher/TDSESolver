#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_X(){
    _apply_mask = &WF::apply_mask_X;
}

void WF::apply_mask_X(cdouble *imask, cdouble *jmask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        _wf[i][0][0] = _wf[i][0][0]*imask[i];
    }
}

cdouble WF::_norm2_X(){
	cdouble integral = 0.0;
	for(int i=0; i<_ni;i++){
		integral += _wf[i][0][0]*conj(_wf[i][0][0])*_di;
	}
	return integral;
}
