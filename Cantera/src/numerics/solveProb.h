/**
 * @file solveProb.h
 *       Header file for implicit nonlinear solver with the option of a pseudotransient
 *  (see \ref numerics and class \link Cantera::solveProb solveProb\endlink).
 */
/*
 * $Id: solveSP.h 381 2010-01-15 21:20:41Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef SOLVEPROB_H
#define SOLVEPROB_H
/**
 * @defgroup solverGroup Solvers for Equation Systems
 */


#include <vector>
#include "Array.h"
#include "ResidEval.h"

//! Solution Methods 
/*!
 * Flag to specify the solution method
 *
 *  1: SOLVEPROB_INITIALIZE   = This assumes that the initial guess supplied to the
 *                          routine is far from the correct one. Substantial
 *                          work plus transient time-stepping is to be expected
 *                          to find a solution.
 *  2:  SOLVEPROB_RESIDUAL    = Need to solve the surface problem in order to
 *                          calculate the surface fluxes of gas-phase species.
 *                          (Can expect a moderate change in the solution
 *                           vector -> try to solve the system by direct
 *                            methods
 *                           with no damping first -> then, try time-stepping
 *                           if the first method fails)
 *                          A "time_scale" supplied here is used in the
 *                          algorithm to determine when to shut off
 *                          time-stepping.
 *  3:  SOLVEPROB_JACOBIAN    = Calculation of the surface problem is due to the
 *                          need for a numerical jacobian for the gas-problem.
 *                          The solution is expected to be very close to the
 *                          initial guess, and accuracy is needed.
 *  4:  SOLVEPROB_TRANSIENT   = The transient calculation is performed here for an
 *                          amount of time specified by "time_scale".  It is
 *                          not garraunted to be time-accurate - just stable
 *                          and fairly fast. The solution after del_t time is
 *                          returned, whether it's converged to a steady
 *                          state or not.
 */
const int SOLVEPROB_INITIALIZE = 1;
const int SOLVEPROB_RESIDUAL   = 2;
const int SOLVEPROB_JACOBIAN   = 3;
const int SOLVEPROB_TRANSIENT  = 4;




namespace Cantera {


  //! Method to solve a pseudo steady state of a nonlinear problem
  /*!
   *   The following class handles solving nonlinear problem.s
   *   
   *
   *   Note there are a couple of different types of species indecices
   *   floating around in the formulation of this object.
   *
   *
   *
   *  Solution Method
   *
   *  This routine is typically used within a residual calculation in a large code.
   *  It's typically invoked millions of times for large calculations, and it must
   *  work every time. Therefore, requirements demand that it be robust but also
   *  efficient.
   * 
   *  The solution methodology is largely determined by the <TT>ifunc<\TT> parameter,
   *  that is input to the solution object. This parameter may have the following
   *  4 values:
   *
   *
   *  1: SFLUX_INITIALIZE   = This assumes that the initial guess supplied to the
   *                          routine is far from the correct one. Substantial
   *                          work plus transient time-stepping is to be expected
   *                          to find a solution.
   *
   *  2:  SFLUX_RESIDUAL    = Need to solve the nonlinear problem in order to
   *                          calculate quantities for a residual calculation
   *                          (Can expect a moderate change in the solution
   *                           vector -> try to solve the system by direct methods
   *                           with no damping first -> then, try time-stepping
   *                           if the first method fails)
   *                          A "time_scale" supplied here is used in the
   *                          algorithm to determine when to shut off
   *                          time-stepping.
   *
   *  3:  SFLUX_JACOBIAN    = Calculation of the surface problem is due to the
   *                          need for a numerical jacobian for the gas-problem.
   *                          The solution is expected to be very close to the
   *                          initial guess, and extra accuracy is needed because
   *                          solution variables have been delta'd from
   *                          nominal values to create jacobian entries.
   *
   *  4:  SFLUX_TRANSIENT   = The transient calculation is performed here for an
   *                          amount of time specified by "time_scale".  It is
   *                          not garraunted to be time-accurate - just stable
   *                          and fairly fast. The solution after del_t time is
   *                          returned, whether it's converged to a steady
   *                          state or not. This is a poor man's time stepping
   *                          algorithm.
   *
   * Psuedo time stepping algorithm:
   *  The time step is determined from sdot[],  so that the time step
   *   doesn't ever change the value of a variable by more than 100%.
   *
   *  This algorithm does use a damped Newton's method to relax the equations.
   *  Damping is based on a "delta damping" technique. The solution unknowns
   *  are not allowed to vary too much between iterations.
   *
   *
   *   EXTRA_ACCURACY:A constant that is the ratio of the required update norm in
   *    this Newton iteration compared to that in the nonlinear solver.
   *     A value of 0.1 is used so surface species are safely  overconverged.
   *
   *  Functions called:
   *----------------------------------------------------------------------------
   *
   * ct_dgetrf    -- First half of LAPACK direct solve of a full Matrix
   *
   * ct_dgetrs    -- Second half of LAPACK direct solve of a full matrix. Returns
   *                 solution vector in the right-hand-side vector, resid.
   *
   *----------------------------------------------------------------------------
   *
   *  @ingroup solverGroup
   */
  class solveProb {

  public:

    //! Constructor for the object
    /*!
     *  @param surfChemPtr  Pointer to the ImplicitSurfChem object that
     *                      defines the surface problem to be solved.
     *
     *  @param bulkFunc     Integer representing how the bulk phases 
     *                      should be handled. Currently, only the
     *                      default value of BULK_ETCH is supported.
     */ 
    solveProb(ResidEval* resid);
   
    //! Destructor. Deletes the integrator.
    ~solveProb();

  private:

    //! Unimplemented private copy constructor
    solveProb(const solveProb &right);

    //! Unimplemented private assignment operator
    solveProb& operator=(const solveProb &right);

  public:
  
    //! Main routine that actually calculates the pseudo steady state
    //! of the surface problem
    /*!
     *   The actual converged solution is returned as part of the
     *   internal state of the InterfaceKinetics objects.
     *
     * @param ifunc Determines the type of solution algorithm to be
     *                  used.  Possible values are  SFLUX_INITIALIZE  ,
     *                  SFLUX_RESIDUAL SFLUX_JACOBIAN  SFLUX_TRANSIENT   .
     *
     * @param time_scale  Time over which to integrate the surface equations,
     *                    where applicable
     * 
     * @param reltol      Relative tolerance to use
     * @param abstol      absolute tolerance.
     *
     * @return  Returns 1 if the surface problem is successfully solved.
     *          Returns -1 if the surface problem wasn't solved successfully.
     *          Note the actual converged solution is returned as part of the
     *          internal state of the InterfaceKinetics objects.
     */
    int solve(int ifunc, doublereal time_scale,
	      doublereal reltol, doublereal abstol);

  private:

    //! Printing routine that gets called at the start of every
    //! invocation
    virtual void print_header(int ioflag, int ifunc, doublereal time_scale,
			      int damping, doublereal reltol, doublereal abstol,
			      doublereal netProdRate[]);

#ifdef DEBUG_SOLVEPROB

    virtual void printResJac(int ioflag, int neq, const Array2D &Jac,
			     doublereal resid[], doublereal wtResid[], doublereal norm);
#endif

    //! Printing routine that gets called after every iteration
    virtual void printIteration(int ioflag, doublereal damp, int label_d, int label_t,
				doublereal inv_t, doublereal t_real, int iter,
				doublereal update_norm, doublereal resid_norm,
				doublereal netProdRate[], doublereal CSolnSP[],
				doublereal resid[], 
				doublereal wtSpecies[], int dim, bool do_time);


    //! Print a summary of the solution
    /*!
     *
     */
    virtual void printFinal(int ioflag, doublereal damp, int label_d, int label_t,
			    doublereal inv_t, doublereal t_real, int iter,
			    doublereal update_norm, doublereal resid_norm,
			    doublereal netProdRateKinSpecies[], const doublereal CSolnSP[],
			    const doublereal resid[],  
			    const doublereal wtSpecies[], const doublereal wtRes[],
			    int dim, bool do_time);
    
    //! Calculate a conservative delta T to use in a pseudo-steady state
    //! algorithm
    /*!
     *    This routine calculates a pretty conservative 1/del_t based
     *    on  MAX_i(sdot_i/(X_i*SDen0)).  This probably guarantees
     *    diagonal dominance.
     *
     *     Small surface fractions are allowed to intervene in the del_t
     *     determination, no matter how small.  This may be changed.
     *     Now minimum changed to 1.0e-12,
     *
     *     Maximum time step set to time_scale.
     *
     *    @param netProdRateSolnSP  Output variable. Net production rate 
     *             of all of the species in the solution vector.
     *    @param XMolSolnSP output variable. 
     *            Mole fraction of all of the species in the  solution vector
     *    @param label Output variable. Pointer to the value of the
     *                 species index (kindexSP) that is controlling
     *                 the time step
     *    @param label_old Output variable. Pointer to the value of the
     *                 species index (kindexSP) that controlled
     *                 the time step at the previous iteration
     *    @param label_factor Output variable. Pointer to the current
     *                 factor that is used to indicate the same species
     *                 is controlling the time step.
     *
     *    @param ioflag Level of the output requested.
     *
     *    @return  Returns the 1. /  delta T to be used on the next step
     */
    virtual doublereal calc_t(doublereal netProdRateSolnSP[], doublereal Csoln[],
			      int *label, int *label_old,
			      doublereal *label_factor, int ioflag);

    //! Calculate the solution and residual weights
    /*!
     *   @param wtSpecies Weights to use for the soln unknowns. These
     *                    are in concentration units
     *  @param wtResid    Weights to sue for the residual unknowns.
     * 
     *  @param CSolnSP    Solution vector for the surface problem
     */
    virtual void calcWeights(doublereal wtSpecies[], doublereal wtResid[],
			     const doublereal CSolnSP[]);

#ifdef DEBUG_SOLVEPROB
    //! Utility routine to print a header for high lvls of debugging
    /*!
     *  @param ioflag Lvl of debugging
     *  @param damp   lvl of damping
     *  @param inv_t  Inverse of the value of delta T
     *  @param t_real Value of the time
     *  @param iter   Interation number
     *  @param do_time boolean indicating whether time stepping is taking
     *                 place
     */
    virtual void printIterationHeader(int ioflag, doublereal damp,
				      doublereal inv_t, doublereal t_real, int iter,
				      bool do_time);
#endif

    /**
     * Update the surface states of the surface phases.
     */
    virtual void updateState(const doublereal *cSurfSpec);


   
    //! Main Function evalulation
    /*!
     *
     *  @param resid output Vector of residuals, length = m_neq     
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param CSolnSPOld Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem. 
     */
    virtual void  fun_eval(doublereal* resid, const doublereal *CSolnSP, 
			   const doublereal *CSolnOldSP,  const bool do_time, const doublereal deltaT);

    //! Main routine that calculates the current residual and Jacobian 
    /*!
     *  @param JacCol  Vector of pointers to the tops of columns of the
     *                 Jacobian to be evalulated.
     *  @param resid   output Vector of residuals, length = m_neq     
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq. These are tweaked in order
     *                  to derive the columns of the jacobian.
     *  @param CSolnSPOld Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem. 
     */
    virtual void resjac_eval(std::vector<doublereal *>& JacCol, doublereal * resid, 
			     doublereal *CSolnSP, 
			     const doublereal *CSolnSPOld,  const bool do_time, 
			     const doublereal deltaT);

    virtual doublereal calc_damping(doublereal x[], doublereal dxneg[], int dim, int *label);

    ResidEval *m_residFunc;

    //! Total number of equations to solve in the implicit problem.
    /*!
     * Note, this can be zero, and frequently is
     */
    int  m_neq;

    //! m_atol is the absolute tolerance in real units.
    vector_fp m_atol;
    
    //! m_rtol is the relative error tolerance.
    doublereal m_rtol;
  
    //! maximum value of the time step
    /*!
     * units = seconds
     */
    doublereal m_maxstep;        
    
    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_netProductionRatesSave;

    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_numEqn1;

    //! Temporary vector with  length MAX(1, m_neq)
    vector_fp m_numEqn2;

    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_CSolnSave;

    //! Solution vector
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSP;

    //! Saved inital solution vector
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSPInit;

    //! Saved  solution vector at the old time step
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSPOld;

    //!  Weights for the residual norm calculation
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_wtResid;

    //!  Weights for the species concentrations norm calculation
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_wtSpecies;

    //!  Residual for the surface problem
    /*!
     *  The residual vector of length "dim" that, that has the value
     *  of "sdot" for surface species.  The residuals for the bulk
     *  species are a function of the sdots for all species in the bulk
     *  phase. The last residual of each phase enforces {Sum(fractions)
     *  = 1}. After linear solve (dgetrf_ & dgetrs_), resid holds the
     *  update vector.
     *
     * length MAX(1, m_neq)
     */
    vector_fp m_resid;

    //!  pivots
    /*!
     * length MAX(1, m_neq)
     */
    vector_int m_ipiv;

    //! Vector of pointers to the top of the columns of the
    //! jacobians
    /*!
     *   The "dim" by "dim" computed Jacobian matrix for the
     *   local Newton's method.
     */
    std::vector<doublereal *> m_JacCol;

    //! Jacobian
    /*!
     *   m_neq by m_neq computed Jacobian matrix for the
     *   local Newton's method.
     */
    Array2D m_Jac;


  public:
    int m_ioflag;
  };
}
#endif
