/**
 *  @file GasKinetics.cpp
 *
 * Homogeneous kinetics in ideal gases
 */

// Copyright 2001  California Institute of Technology
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Enhanced3BConc.h"
#include "cantera/kinetics/ThirdBodyMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/TurbulentKinetics.h"

using namespace std;

namespace Cantera
{

void TurbulentKinetics::update_rates_T()
{
    doublereal T = thermo().temperature();
    doublereal P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());
    doublereal logT = log(T);
	doublereal TempFluc = Tprime();

    if (T != m_temp) {
        if (!m_rfn.empty()) {
            m_rates.updateTurb(T, logT, &m_rfn[0],TempFluc);
        }

        if (!m_rfn_low.empty()) {
            m_falloff_low_rates.updateTurb(T, logT, &m_rfn_low[0],TempFluc);
            m_falloff_high_rates.updateTurb(T, logT, &m_rfn_high[0],TempFluc);
        }
        if (!falloff_work.empty()) {
            m_falloffn.updateTemp(T, &falloff_work[0]);
        }
        updateKc();
        m_ROP_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.updateTurb(T, logT, &m_rfn[0],TempFluc);
            m_ROP_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.updateTurb(T, logT, &m_rfn[0],TempFluc);
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}


}
