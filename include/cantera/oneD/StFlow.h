//! @file StFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef ct_stflow_h
#define ct_stflow_h

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/ct_defs.h"

#include <algorithm> // fill,min
#include <cmath> // exp,log, sqrt, pow
#include <numeric> // Quick sum
#include <string> // Spray flame prototype


namespace Cantera
{

//------------------------------------------
//   constants
//------------------------------------------

// Offsets of solution components in the solution array.
const size_t c_offset_U = 0; // axial velocity
const size_t c_offset_V = 1; // strain rate
const size_t c_offset_T = 2; // temperature
const size_t c_offset_L = 3; // (1/r)dP/dr
const size_t c_offset_Y = 4; // mass fractions

// Parameters used in spray flame
const size_t c_offset_Ul = 0; // liquid radial velocity 
const size_t c_offset_vl = 1; // liquid axial velocity
const size_t c_offset_Tl = 2; // liquid temperature
const size_t c_offset_nl = 3; // number density
const size_t c_offset_ml = 4; // droplet species masses
const doublereal mmHg2Pa = 133.322365;
const doublereal bar2Pa  = 1.0e+05;

class Transport;

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric flows.
 *  @ingroup onedim
 */
class StFlow : public Domain1D
{

friend class ModelGeneric;

public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! @name Problem Specification
    //! @{

    virtual void setupGrid(size_t n, const doublereal* z);

    virtual void resetBadValues(double* xg);


    thermo_t& phase() {
        return *m_thermo;
    }
    Kinetics& kinetics() {
        return *m_kin;
    }

    /**
     * Set the thermo manager. Note that the flow equations assume
     * the ideal gas equation.
     */
    void setThermo(IdealGasPhase& th) {
        m_thermo = &th;
    }

    //! Set the kinetics manager. The kinetics manager must
    void setKinetics(Kinetics& kin) {
        m_kin = &kin;
    }

    //! set the transport manager
    void setTransport(Transport& trans);

    void enableSoret(bool withSoret);
    bool withSoret() const {
        return m_do_soret;
    }

    //! Set the pressure. Since the flow equations are for the limit of small
    //! Mach number, the pressure is very nearly constant throughout the flow.
    void setPressure(doublereal p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    doublereal pressure() const {
        return m_press;
    }

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(double* x);

    virtual void _finalize(const doublereal* x);

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /*!
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, doublereal t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    //! The fixed temperature value at point j.
    doublereal T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    // @}

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    //! Print the solution.
    virtual void showSolution(const doublereal* x);

    //! Save the current solution for this domain into an XML_Node
    /*!
     *  @param o    XML_Node to save the solution to.
     *  @param sol  Current value of the solution vector. The object will Pick
     *              out which part of the solution vector pertains to this
     *              object.
     */
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    virtual void restore(const XML_Node& dom, doublereal* soln,
                         int loglevel);

    // overloaded in subclasses
    virtual std::string flowType() {
        return "<none>";
    }

    void solveEnergyEqn(size_t j=npos);

    //! Turn radiation on / off.
    /*!
     *  The simple radiation model used was established by Y. Liu and B. Rogg
     *  [Y. Liu and B. Rogg, Modelling of thermally radiating diffusion flames
     *  with detailed chemistry and transport, EUROTHERM Seminars, 17:114-127,
     *  1991]. This model considers the radiation of CO2 and H2O.
     */
    void enableRadiation(bool doRadiation) {
        m_do_radiation = doRadiation;
    }

    //! Returns `true` if the radiation term in the energy equation is enabled
    bool radiationEnabled() const {
        return m_do_radiation;
    }

    //! Set the emissivities for the boundary values
    /*!
     * Reads the emissivities for the left and right boundary values in the
     * radiative term and writes them into the variables, which are used for the
     * calculation.
     */
    void setBoundaryEmissivities(doublereal e_left, doublereal e_right);

    void fixTemperature(size_t j=npos);

    bool doEnergy(size_t j) {
        return m_do_energy[j];
    }

    //! Change the grid size. Called after grid refinement.
    void resize(size_t components, size_t points);

    virtual void setFixedPoint(int j0, doublereal t0) {}

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j);

    doublereal density(size_t j) const {
        return m_rho[j];
    }

    virtual bool fixed_mdot() {
        return true;
    }
    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  j == npos, the residual function is evaluated at all grid points.
     *  Otherwise, the residual function is only evaluated at grid points
     *  j-1, j, and j+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt) = 0;

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt) = 0;

    //! Index of the species on the left boundary with the largest mass fraction
    size_t leftExcessSpecies() const {
        return m_kExcessLeft;
    }

    //! Index of the species on the right boundary with the largest mass fraction
    size_t rightExcessSpecies() const {
        return m_kExcessRight;
    }

protected:
    doublereal wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(doublereal* x, size_t j) {
        setGas(x,j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(const doublereal* x, size_t j0, size_t j1) {
        for (size_t j = j0; j <= j1; j++) {
            setGas(x,j);
            m_rho[j] = m_thermo->density();
            m_wtm[j] = m_thermo->meanMolecularWeight();
            m_cp[j] = m_thermo->cp_mass();
        }
    }

    //! @name Solution components
    //! @{

    doublereal T(const doublereal* x, size_t j) const {
        return x[index(c_offset_T, j)];
    }
    doublereal& T(doublereal* x, size_t j) {
        return x[index(c_offset_T, j)];
    }
    doublereal T_prev(size_t j) const {
        return prevSoln(c_offset_T, j);
    }

    doublereal rho_u(const doublereal* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U, j)];
    }

    doublereal u(const doublereal* x, size_t j) const {
        return x[index(c_offset_U, j)];
    }

    doublereal V(const doublereal* x, size_t j) const {
        return x[index(c_offset_V, j)];
    }
    doublereal V_prev(size_t j) const {
        return prevSoln(c_offset_V, j);
    }

    doublereal lambda(const doublereal* x, size_t j) const {
        return x[index(c_offset_L, j)];
    }

    doublereal Y(const doublereal* x, size_t k, size_t j) const {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal& Y(doublereal* x, size_t k, size_t j) {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y + k, j);
    }

    doublereal X(const doublereal* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    doublereal flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }
    //! @}

    //! @name convective spatial derivatives.
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    doublereal dVdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dYdz(const doublereal* x, size_t k, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
    }

    doublereal dTdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    doublereal shear(const doublereal* x, size_t j) const {
        doublereal c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
        doublereal c2 = m_visc[j]*(V(x,j+1) - V(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal divHeatFlux(const doublereal* x, size_t j) const {
        doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
        doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //! Update the diffusive mass fluxes.
    void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    doublereal m_press; // pressure

    // grid parameters
    vector_fp m_dz;

    // mixture thermo properties
    vector_fp m_rho;
    vector_fp m_wtm;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_cp;

    // transport properties
    vector_fp m_visc;
    vector_fp m_tcon;
    vector_fp m_diff;
    vector_fp m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;

    // production rates
    Array2D m_wdot;

    size_t m_nsp;

    IdealGasPhase* m_thermo;
    Kinetics* m_kin;
    Transport* m_trans;

    // boundary emissivities for the radiation calculations
    doublereal m_epsilon_left;
    doublereal m_epsilon_right;

    //! Indices within the ThermoPhase of the radiating species. First index is
    //! for CO2, second is for H2O.
    std::vector<size_t> m_kRadiating;

    // flags
    std::vector<bool> m_do_energy;
    bool m_do_soret;
    std::vector<bool> m_do_species;
    bool m_do_multicomponent;

    //! flag for the radiative heat loss
    bool m_do_radiation;

    //! radiative heat loss vector
    vector_fp m_qdotRadiation;

    // fixed T and Y values
    vector_fp m_fixedtemp;
    vector_fp m_zfix;
    vector_fp m_tfix;

    //! Index of species with a large mass fraction at each boundary, for which
    //! the mass fraction may be calculated as 1 minus the sum of the other mass
    //! fractions
    size_t m_kExcessLeft;
    size_t m_kExcessRight;

    bool m_dovisc;

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    void updateTransport(const doublereal* x, size_t j0, size_t j1);

    // Update transport properties in range. Inclusive. Use gas at location and used in droplet model
    void updateTransportLocal(const doublereal* x, size_t j0, size_t j1);

private:
    vector_fp m_ybar;

};

/**
 * A class for axisymmetric stagnation flows.
 * @ingroup onedim
 */
class AxiStagnFlow : public StFlow
{
public:
    AxiStagnFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1) :
        StFlow(ph, nsp, points) {
        m_dovisc = true;
    }

    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Axisymmetric Stagnation";
    }
};

/**
 * A class for freely-propagating premixed flames.
 * @ingroup onedim
 */
class FreeFlame : public StFlow
{
public:
    FreeFlame(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);
    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Free Flame";
    }
    virtual bool fixed_mdot() {
        return false;
    }
    virtual void _finalize(const doublereal* x);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    //! Location of the point where temperature is fixed
    doublereal m_zfixed;

    //! Temperature at the point used to fix the flame location
    doublereal m_tfixed;
};

/**
 * A class for spray flames.
 * @ingroup onedim
 */
class SprayFlame : public AxiStagnFlow
{

friend class SprayInlet1D;
friend class SprayOutlet1D;
friend class SprayPhase;

public:
    SprayFlame(IdealGasPhase* ph = 0, std::string fuel = "", std::vector<std::string> palette = {""}, std::string evapModel = "", size_t nsp = 1, size_t points = 1);

    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    virtual void evalNumberDensity(size_t j, doublereal* x, doublereal* rsd,
                      integer* diag, doublereal rdt);

    virtual void evalRightBoundaryLiquid(doublereal* x, doublereal* rsd,
                      integer* diag, doublereal rdt);

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    virtual std::string flowType() {
        return "Axisymmetric Spray Stagnation";
    }

    void updateFuelSpecies(const std::vector<std::string> palette) {
        numFuelSpecies = palette.size();
	m_palette = palette;
	c_offset_fuel.resize(numFuelSpecies);
        for (size_t i = 0; i < numFuelSpecies; i++)
          c_offset_fuel[i] = componentIndex(palette[i]) - c_offset_Y;
    }

    size_t getNumFuelSpecies() {
         return numFuelSpecies;
    }

   
    std::vector<doublereal> getMW();

    std::vector<doublereal> getRhoL(double T); 


protected:

    //! @name Solution components
    //! @{
    doublereal Tl(size_t j) const {
         return m_sprayPhase->Tl(j); 
    }

    doublereal Tl(const doublereal* x, size_t j) const {
        return x[index(c_offset_Y+m_nsp+c_offset_Tl,j)];
    }

    doublereal Tl_prev(size_t j) const {
        return prevSoln(c_offset_Y+m_nsp+c_offset_Tl, j);
    }

    doublereal vl(const doublereal* x, size_t j) const {
        return x[index(c_offset_Y+m_nsp+c_offset_vl,j)];
    }

    doublereal vl_prev(size_t j) const {
        return prevSoln(c_offset_Y+m_nsp+c_offset_vl, j);
    }

    doublereal Ul(const doublereal* x, size_t j) const {
        return x[index(c_offset_Y+m_nsp+c_offset_Ul,j)];
    }

    doublereal Ul_prev(size_t j) const {
        return prevSoln(c_offset_Y+m_nsp+c_offset_Ul, j);
    }

    doublereal ml(const doublereal* x, size_t j) const {
        doublereal ret = 0.0;
        for (size_t i = 0; i < numFuelSpecies; i++)
             ret += x[index(c_offset_Y+m_nsp+c_offset_ml+i,j)];
        return ret;
    }

    doublereal mlk(const doublereal* x, size_t k, size_t j) const {
        if (k >= numFuelSpecies) return 0.0;
        return x[index(c_offset_Y+m_nsp+c_offset_ml+k,j)];
    }

    doublereal ml_prev(size_t j) const {
        doublereal ret = 0.0;
        for (size_t i = 0; i < numFuelSpecies; i++)
             ret += prevSoln(c_offset_Y+m_nsp+c_offset_ml+i,j);
        return ret;
    }

    doublereal mlk_prev(size_t i, size_t j) const {
         return prevSoln(c_offset_Y+m_nsp+c_offset_ml+i,j); 
    }

    doublereal nl(size_t j) const {
         return m_sprayPhase->nl(j); 
    }

    std::vector<doublereal> mkdot(size_t k, size_t j) const {
         return m_sprayPhase->mkdot(k,j); 
    }
    
    doublereal mdot(size_t j) const {
         return m_sprayPhase->mdot(j); 
    }

    doublereal q(size_t j) const {
         return m_sprayPhase->q(j); 
    }

    doublereal nl(const doublereal* x, size_t j) const {
        return x[index(c_offset_Y+m_nsp+c_offset_nl,j)];
        // return 0.0;
    }

    doublereal nl_prev(size_t j) const {
        return prevSoln(c_offset_Y+m_nsp+c_offset_nl, j);
    }

    doublereal rhol(const doublereal* x, size_t j) const;

    doublereal ml_vl(const doublereal* x, size_t j) const {
        return ml(x,j)*vl(x,j);
    }

    doublereal ml_Ul(const doublereal* x, size_t j) const {
        return ml(x,j)*Ul(x,j);
    }

    doublereal nl_Ul(const doublereal* x, size_t j) const {
        return nl(x,j)*Ul(x,j);
    }

    doublereal nl_vl(const doublereal* x, size_t j) const {
        return nl(x,j)*vl(x,j);
    }

    doublereal dl(const doublereal* x, size_t j) const {
         return m_sprayPhase->dl(x,j); 
    }

    doublereal dl(const doublereal* x, size_t j) const {
        if (ml(x,j) <= 0.0 || rhol(x,j) <= 0.0) 
            return 0.0;
        else
            return std::pow(6.0*ml(x,j)/Pi/rhol(x,j),1.0/3.0);
    }

    doublereal cpl(const doublereal* x, size_t j) {
        return m_sprayPhase->cpl(x,j);
    }

    doublereal cpgf(const doublereal* x,size_t j) {
        // setGas(x,j);
        // doublereal Yr = 2.0/3.0*Yrs(x,j) + 1.0/3.0*Y(x,c_offset_fuel,j);
        // vector_fp cp_R = m_thermo->cp_R_ref();
        // return Yr*GasConstant*cp_R[c_offset_fuel] + (1.0-Yr)*m_cp[j];
        return m_cp[j];
    }
    
    std::vector<double> source(const doublereal* x, size_t j);

    doublereal Fr(const doublereal* x, size_t j) {
        return 3.0*Pi*dl(j)*m_visc[j]*(V(x,j)-Ul(j));
    }

    doublereal fz(const doublereal* x, size_t j) {
        return 3.0*Pi*dl(x,j)*m_visc[j]*(u(x,j)-vl(x,j));
    }

    //! @}
    
    //! @name convective spatial derivatives.
    //! These use upwind differencing, assuming vl(z) is negative
    //! @{
    doublereal dUldz(const doublereal* x, size_t j) const {
        size_t jloc = (vl(x,j) > 0.0 ? j : j + 1);
        return (Ul(x,jloc) - Ul(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dvldz(const doublereal* x, size_t j) const {
        size_t jloc = (vl(x,j) > 0.0 ? j : j + 1);
        return (vl(x,jloc) - vl(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dmlkdz(const doublereal* x, size_t k, size_t j) const {
        size_t jloc = (vl(x,j) > 0.0 ? j : j + 1);
 	if (k >= numFuelSpecies) return 0.0;	
        return (mlk(x,k,jloc) - mlk(x,k,jloc-1))/m_dz[jloc-1];
    }

    doublereal dnldz(const doublereal* x, size_t j) const {
        size_t jloc = (vl(x,j) > 0.0 ? j : j + 1);
        return (nl(x,jloc) - nl(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dTldz(const doublereal* x, size_t j) const {
        size_t jloc = (vl(x,j) > 0.0 ? j : j + 1);
        return (Tl(x,jloc) - Tl(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    //! @name artifitial viscosities
    //! @{
    doublereal av_mlk(const doublereal* x, size_t k, size_t j) const {
        doublereal m_visc_ml = 5.0e-0; // ml(x,j)/dl(x,j);
        doublereal c1 = m_visc_ml*(mlk(x,k,j) - mlk(x,k,j-1));
        doublereal c2 = m_visc_ml*(mlk(x,k,j+1) - mlk(x,k,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal av_nl(const doublereal* x, size_t j) const {
        doublereal m_visc_nl = 1.0e-5;
        doublereal c1 = m_visc_nl*(nl(x,j) - nl(x,j-1));
        doublereal c2 = m_visc_nl*(nl(x,j+1) - nl(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal av_Tl(const doublereal* x, size_t j) const {
        doublereal m_visc_Tl = 1.0e-2;
        doublereal c1 = m_visc_Tl*(Tl(x,j) - Tl(x,j-1));
        doublereal c2 = m_visc_Tl*(Tl(x,j+1) - Tl(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal av_Ul(const doublereal* x, size_t j) const {
        doublereal m_visc_Ul = 1.0e-5;
        doublereal c1 = m_visc_Ul*(Ul(x,j) - Ul(x,j-1));
        doublereal c2 = m_visc_Ul*(Ul(x,j+1) - Ul(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal av_vl(const doublereal* x, size_t j) const {
        doublereal m_visc_vl = 1.0e-5;
        doublereal c1 = m_visc_vl*(vl(x,j) - vl(x,j-1));
        doublereal c2 = m_visc_vl*(vl(x,j+1) - vl(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }
    //! @}

    // fuel species index
    std::vector<size_t> c_offset_fuel;
    size_t numFuelSpecies;
    std::vector<std::string> m_palette;

    // Evaporation model
    std::string m_evapModel;

    // Liquid heat capacity
    doublereal m_cpl;

    // Ideal gas constant
    doublereal R_u = 8.314; // J/mol/K

    // Underlying spray phase
    SprayPhase* m_SprayPhase;

};

class SprayPhase : public Sim1D {

public:
    SprayPhase();

    virtual void eval;
   
    virtual void advance;

    virtual void advance_to_convergence;

    virtual void match_grids;

    virtual void check_point_vapor;

    virtual bool check_grid_vapor;

    virtual void mdotVec(x,j);

    virtual void mdot(x,j);

    virtual void q(x,j);

    // Throw error if too small 
    virtual void reduce_timestep;

    virtual void reset_timestep; 
   
    virtual void compute_grid_source;

    virtual void nl(x,j);
 
private:
    SprayFlame* m_SprayFlame;
    

};

}

#endif
