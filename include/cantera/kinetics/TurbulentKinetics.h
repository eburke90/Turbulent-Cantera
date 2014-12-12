/**
 * @file TurbulentKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_TURBULENTKINETICS_H
#define CT_TURBULENTKINETICS_H

#include "cantera/thermo/mix_defs.h"
#include "Kinetics.h"
#include "cantera/base/utilities.h"
#include "cantera/kinetics/GasKinetics.h"
#include "ReactionStoichMgr.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"

namespace Cantera
{

class TurbulentKinetics : public GasKinetics {
public:
    virtual int type() const {
        return cTurbulentKinetics;
    }

    void setTprime(double Tprime) {
        m_Tprime = Tprime;
    }
	
	doublereal Tprime() const {
		return m_Tprime;
   }
   virtual void update_rates_T();

protected:
    double m_Tprime;
};
}
#endif