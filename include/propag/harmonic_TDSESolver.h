#ifndef HARMONIC_TDSESOLVER_H
#define HARMONIC_TDSESOLVER_H

#include "propag/harmonic.h"
#include "tdsesolver.h"

class HarmonicTDSESolver : public Harmonic{
    public:
        HarmonicTDSESolver(std::shared_ptr<Settings> settings);
        virtual ~HarmonicTDSESolver();

        virtual void calculateAcceleration(Field *field);
    
    protected:
        Parameters *param;
        TDSESolver *tdsesolver;
};

#endif
