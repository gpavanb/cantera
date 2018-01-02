//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/gc_wrap.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"

#include <algorithm> // accumulate, copy_n

using namespace std;

namespace Cantera
{

StFlow::StFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_epsilon_left(0.0),
    m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_kExcessLeft(0),
    m_kExcessRight(0)
{
    m_type = cFlowType;
    m_points = points;
    m_thermo = ph;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = c_offset_Y + m_nsp;

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(0, -1e20, 1e20); // no bounds on u
    setBounds(1, -1e20, 1e20); // V
    setBounds(2, 200.0, 1e9); // temperature bounds
    setBounds(3, -1e20, 1e20); // lambda should be negative

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    m_refiner->setActive(3, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());
    setID("stagnation flow");

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");

}

void StFlow::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);

}

void StFlow::setupGrid(size_t n, const doublereal* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("StFlow::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }

}

void StFlow::resetBadValues(double* xg) {
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void StFlow::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void StFlow::enableSoret(bool withSoret)
{
    if (m_do_multicomponent) {
        m_do_soret = withSoret;
    } else {
        throw CanteraError("setTransport",
                           "Thermal diffusion (the Soret effect) "
                           "requires using a multicomponent transport model.");
    }
}

void StFlow::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
    }
}

void StFlow::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void StFlow::_finalize(const doublereal* x)
{
    size_t nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (size_t j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            m_fixedtemp[j] = T(x, j);
        } else {
            double zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            double tt = linearInterp(zz, m_zfix, m_tfix);
            m_fixedtemp[j] = tt;
        }
    }
    if (e) {
        solveEnergyEqn();
    }
}

void StFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // if evaluating a Jacobian, compute the steady-state residual
    if (jg != npos) {
        rdt = 0.0;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    // ------------ update properties ------------

    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
    }
    if (jg == npos) {
        double* Yleft = x + index(c_offset_Y, jmin);
        m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
        double* Yright = x + index(c_offset_Y, jmax);
        m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);

    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const doublereal c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const doublereal c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            rsd[index(c_offset_T,0)] = T(x,0);
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                double wdot_;
                wdot_ = this->wdot(k,j);
                
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot_)
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {
                setGas(x,j);

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                double dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::updateTransport(const doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void StFlow::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

void StFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                doublereal sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            double wtm = m_wtm[j];
            double rho = density(j);
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum*Y(x,k,j);
            }
        }
    }

    if (m_do_soret) {
        for (size_t m = j0; m < j1; m++) {
            double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                              ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}

string StFlow::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "u";
    case 1:
        return "V";
    case 2:
        return "T";
    case 3:
        return "lambda";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const std::string& name) const
{
    if (name=="u") {
        return 0;
    } else if (name=="V") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
    }
    return npos;
}

void StFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    vector<string> ignored;
    size_t nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);
    vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_fp x;
    size_t np = 0;
    bool readgrid = false, wrote_header = false;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            readgrid = true;
            setupGrid(np, x.data());
        }
    }
    if (!readgrid) {
        throw CanteraError("StFlow::restore",
                           "domain contains no grid points.");
    }

    debuglog("Importing datasets:\n", loglevel >= 2);
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(0,j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V") {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(1,j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(2,j)] = x[j];
            }

            // For fixed-temperature simulations, use the imported temperature
            // profile by default.  If this is not desired, call
            // setFixedTempProfile *after* restoring the solution.
            vector_fp zz(np);
            for (size_t jj = 0; jj < np; jj++) {
                zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
            }
            setFixedTempProfile(zz, x);
        } else if (nm == "L") {
            debuglog("lambda   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "lambda arary size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(3,j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (loglevel >=2 && !ignored.empty()) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        size_t nn = ignored.size();
        for (size_t n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    if (loglevel >= 1) {
        for (size_t ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }
    }

    if (dom.hasChild("energy_enabled")) {
        getFloatArray(dom, x, false, "", "energy_enabled");
        if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_energy[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            throw CanteraError("StFlow::restore", "energy_enabled is length {}"
                               "but should be length {}", x.size(), nPoints());
        }
    }

    if (dom.hasChild("species_enabled")) {
        getFloatArray(dom, x, false, "", "species_enabled");
        if (x.size() == m_nsp) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_species[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            // This may occur when restoring from a mechanism with a different
            // number of species.
            if (loglevel > 0) {
                writelog("\nWarning: StFlow::restore: species_enabled is "
                    "length {} but should be length {}. Enabling all species "
                    "equations by default.", x.size(), m_nsp);
            }
            m_do_species.assign(m_nsp, true);
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }
}

XML_Node& StFlow::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    flow.addAttribute("type",flowType());

    if (m_desc != "") {
        addString(flow,"description",m_desc);
    }
    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");

    addFloatArray(gv,"z",m_z.size(), m_z.data(),
                  "m","length");
    vector_fp x(soln.nColumns());

    soln.getRow(0, x.data());
    addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

    soln.getRow(1, x.data());
    addFloatArray(gv,"V",x.size(),x.data(),"1/s","rate");

    soln.getRow(2, x.data());
    addFloatArray(gv,"T",x.size(),x.data(),"K","temperature");

    soln.getRow(3, x.data());
    addFloatArray(gv,"L",x.size(),x.data(),"N/m^4");

    for (size_t k = 0; k < m_nsp; k++) {
        soln.getRow(c_offset_Y+k, x.data());
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),x.data(),"","massFraction");
    }
    if (m_do_radiation) {
        addFloatArray(gv, "radiative_heat_loss", m_z.size(),
            m_qdotRadiation.data(), "W/m^3", "specificPower");
    }
    vector_fp values(nPoints());
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = m_do_energy[i];
    }
    addNamedFloatArray(flow, "energy_enabled", nPoints(), &values[0]);

    values.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        values[i] = m_do_species[i];
    }
    addNamedFloatArray(flow, "species_enabled", m_nsp, &values[0]);

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    return flow;
}

void StFlow::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
    m_refiner->setActive(0, true);
    m_refiner->setActive(1, true);
    m_refiner->setActive(2, true);
    if (changed) {
        needJacUpdate();
    }
}

void StFlow::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void StFlow::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    if (changed) {
        needJacUpdate();
    }
}

void AxiStagnFlow::evalRightBoundary(doublereal* x, doublereal* rsd,
                                     integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;
    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.
    rsd[index(0,j)] = rho_u(x,j);
    rsd[index(1,j)] = V(x,j);
    rsd[index(2,j)] = T(x,j);
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    doublereal sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
}

void AxiStagnFlow::evalContinuity(size_t j, doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
    //----------------------------------------------
    //    Continuity equation
    //
    //    Note that this propagates the mass flow rate information to the left
    //    (j+1 -> j) from the value specified at the right boundary. The
    //    lambda information propagates in the opposite direction.
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //------------------------------------------------
    rsd[index(c_offset_U,j)] =
        -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
        -(density(j+1)*V(x,j+1) + density(j)*V(x,j));

    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
}

FreeFlame::FreeFlame(IdealGasPhase* ph, size_t nsp, size_t points) :
    StFlow(ph, nsp, points),
    m_zfixed(Undef),
    m_tfixed(Undef)
{
    m_dovisc = false;
    setID("flame");
}

void FreeFlame::evalRightBoundary(doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    // zero gradient
    rsd[index(0,j)] = rho_u(x,j) - rho_u(x,j-1);
    rsd[index(1,j)] = V(x,j);
    rsd[index(2,j)] = T(x,j) - T(x,j-1);
    doublereal sum = 0.0;
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
}

void FreeFlame::evalContinuity(size_t j, doublereal* x, doublereal* rsd,
                               integer* diag, doublereal rdt)
{
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (grid(j) > m_zfixed) {
        rsd[index(c_offset_U,j)] =
            - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1]
            - (density(j-1)*V(x,j-1) + density(j)*V(x,j));
    } else if (grid(j) == m_zfixed) {
        if (m_do_energy[j]) {
            rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
        } else {
            rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                        - m_rho[0]*0.3);
        }
    } else if (grid(j) < m_zfixed) {
        rsd[index(c_offset_U,j)] =
            - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            - (density(j+1)*V(x,j+1) + density(j)*V(x,j));
    }
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
}

void FreeFlame::_finalize(const doublereal* x)
{
    StFlow::_finalize(x);
    // If the domain contains the temperature fixed point, make sure that it
    // is correctly set. This may be necessary when the grid has been modified
    // externally.
    if (m_tfixed != Undef) {
        for (size_t j = 0; j < m_points; j++) {
            if (z(j) == m_zfixed) {
                return; // fixed point is already set correctly
            }
        }

        for (size_t j = 0; j < m_points - 1; j++) {
            // Find where the temperature profile crosses the current
            // fixed temperature.
            if ((T(x, j) - m_tfixed) * (T(x, j+1) - m_tfixed) <= 0.0) {
                m_tfixed = T(x, j+1);
                m_zfixed = z(j+1);
                return;
            }
        }
    }
}

void FreeFlame::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    StFlow::restore(dom, soln, loglevel);
    getOptionalFloat(dom, "t_fixed", m_tfixed);
    getOptionalFloat(dom, "z_fixed", m_zfixed);
}

XML_Node& FreeFlame::save(XML_Node& o, const doublereal* const sol)
{
    XML_Node& flow = StFlow::save(o, sol);
    if (m_zfixed != Undef) {
        addFloat(flow, "z_fixed", m_zfixed, "m");
        addFloat(flow, "t_fixed", m_tfixed, "K");
    }
    return flow;
}

SprayFlame::SprayFlame(IdealGasPhase* ph, std::string fuel, std::vector<std::string> palette, std::string evapModel, size_t nsp, size_t points) :
    AxiStagnFlow(ph, nsp, points)
{
    initialize_gc(fuel);
    size_t ns = numSpecies();
    updateFuelSpecies(palette);
    m_evapModel = evapModel; 

    m_nv = c_offset_Y+m_nsp+c_offset_ml+ns;
    Domain1D::resize(m_nv,points);

    setBounds(c_offset_Y+m_nsp+c_offset_Ul, -1e20, 1e20); // no bounds on Ul
    setBounds(c_offset_Y+m_nsp+c_offset_vl, -1e20, 1e20); // no bounds on vl
    // TODO : Set more precise upper bound, especially for volatile liquids
    setBounds(c_offset_Y+m_nsp+c_offset_Tl, 200.0, 5000.0); // bounds on Tl
    setBounds(c_offset_Y+m_nsp+c_offset_nl, -1e-7, 1e20); // positivity for nl

    setID("spray flame");
    
    // Tighter lower bound for mass
    for (size_t i = 0; i < ns; i++)
        setBounds(c_offset_Y+m_nsp+c_offset_ml+i, 0.0, 1e20); // bounds on ml
}

void SprayFlame::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // Get the residual from the gaseous phase equations
    AxiStagnFlow::eval(jg,xg,rg,diagg,rdt);

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // Gaseous phase
    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
	    std::vector<doublereal> sourceVec_j = source(x,j);
	    std::vector<doublereal> sourceVec_jp1 = source(x,j+1);
	    std::vector<doublereal> mdot_j(numFuelSpecies);
	    std::vector<doublereal> mdot_jp1(numFuelSpecies);
	    std::copy_n(sourceVec_j.begin(),numFuelSpecies,mdot_j.begin());
	    std::copy_n(sourceVec_jp1.begin(),numFuelSpecies,mdot_jp1.begin());

	    doublereal sum_mdot_j = std::accumulate(mdot_j.begin(), mdot_j.end(),0.0);
	    doublereal sum_mdot_jp1 = std::accumulate(mdot_jp1.begin(), mdot_jp1.end(),0.0);

            rsd[index(c_offset_U,0)] += (nl(x,0)*sum_mdot_j + nl(x,1)*sum_mdot_jp1)/2.0;

        } else if (j == m_points - 1) {
            continue;

        } else { // interior points
	    // TODO : Modify comments to reflect new equations

	    std::vector<doublereal> sourceVec_j = source(x,j);
	    std::vector<doublereal> sourceVec_jp1 = source(x,j+1);
	    std::vector<doublereal> mdot_j(numFuelSpecies);
	    std::vector<doublereal> mdot_jp1(numFuelSpecies);
	    std::copy_n(sourceVec_j.begin(),numFuelSpecies,mdot_j.begin());
	    std::copy_n(sourceVec_jp1.begin(),numFuelSpecies,mdot_jp1.begin());
	    doublereal q_j = sourceVec_j[numFuelSpecies];
	    //doublereal q_jp1 = sourceVec_jp1[numFuelSpecies];

	    doublereal sum_mdot_j = std::accumulate(mdot_j.begin(), mdot_j.end(),0.0);
	    doublereal sum_mdot_jp1 = std::accumulate(mdot_jp1.begin(), mdot_jp1.end(),0.0);
            //------------------------------------------------
            //    Continuity equation
            //------------------------------------------------
		    rsd[index(c_offset_U,j)] += (nl(x,j)*sum_mdot_j + nl(x,j+1)*sum_mdot_jp1)/2.0;

		    //------------------------------------------------
		    //    Radial momentum equation
		    //
		    //    \rho dV/dt + \rho u dV/dz + \rho V^2
		    //       = d(\mu dV/dz)/dz - lambda
		    //         + nl mdot (Ul - Ug) - nl Fr
		    //-------------------------------------------------
		    //rsd[index(c_offset_V,j)] -= ( nl(x,j) * Fr(x,j) / m_rho[j] );
		    rsd[index(c_offset_V,j)] += 
			(nl(x,j) * sum_mdot_j * (Ul(x,j)-V(x,j)) - nl(x,j) * Fr(x,j)) / m_rho[j];

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //     + (\delta_kf*mdot[i] - Y_k*mdot) nl 
	    // Refer to Olguin, Gutheil (1) and (10)
            //-------------------------------------------------
            for (size_t k = 0; k < m_nsp; k++) {

                // Identify liquid species number of gas species
		size_t idx = 0;
		for (idx = 0; idx < numFuelSpecies; idx++) {
                    if (k == c_offset_fuel[idx]) 
			break;
                }
                if (idx < numFuelSpecies) {
                     rsd[index(c_offset_Y + k, j)] += 
                         (mdot_j[idx] - Y(x,k,j) * sum_mdot_j) * nl(x,j) / m_rho[j];
                }
                else {
                     rsd[index(c_offset_Y + k, j)] += 
                         (- Y(x,k,j) * sum_mdot_j) * nl(x,j) / m_rho[j];
                }
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //      + nl mdot cp (Tl - Tg) - nl mdot q
            //-----------------------------------------------
            rsd[index(c_offset_T, j)] += 
                (nl(x,j) * sum_mdot_j * cpgf(x,j) * (Tl(x,j) - T(x,j)) - 
                 nl(x,j) * q_j)/ (m_rho[j]*m_cp[j]);

        }
    }

    
    // Liquid phase
    for (size_t j = jmin; j <= jmax; j++) { 
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object
            
            // Number density. This propagates information right-to-left, since
            // ml_vl at point 0 is dependent on ml_vl at point 1, but not on
            // mdot from the inlet.
            // rsd[index(c_offset_Y+m_nsp+c_offset_nl,0)] = 
            //     -(nl_vl(x,1) - nl_vl(x,0))/m_dz[0]
            //     -(nl_Ul(x,1) + nl_Ul(x,0));
            rsd[index(c_offset_Y+m_nsp+c_offset_nl,0)] = nl(x,0); 
            diag[index(c_offset_Y+m_nsp+c_offset_nl, 0)] = 0;

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_Y+m_nsp+c_offset_vl,0)] = vl(x,0);
            rsd[index(c_offset_Y+m_nsp+c_offset_Ul,0)] = Ul(x,0);
            rsd[index(c_offset_Y+m_nsp+c_offset_Tl,0)] = Tl(x,0);
            for (size_t i = 0; i < numFuelSpecies; i++) 
                 rsd[index(c_offset_Y+m_nsp+c_offset_ml+i,0)] = mlk(x,i,0);

            diag[index(c_offset_Y+m_nsp+c_offset_Ul, 0)] = 0;
            diag[index(c_offset_Y+m_nsp+c_offset_vl, 0)] = 0;
            diag[index(c_offset_Y+m_nsp+c_offset_Tl, 0)] = 0;
            for (size_t i = 0; i < numFuelSpecies; i++) 
                 diag[index(c_offset_Y+m_nsp+c_offset_ml+i, 0)] = 0;

        } else if (j == m_points - 1) {
            evalRightBoundaryLiquid(x, rsd, diag, rdt);
        } else { // interior points
	    // TODO : Call source term function only once
	    std::vector<doublereal> sourceVec_j = source(x,j);
	    std::vector<doublereal> sourceVec_jp1 = source(x,j+1);
	    std::vector<doublereal> mdot_j(numFuelSpecies);
	    std::vector<doublereal> mdot_jp1(numFuelSpecies);
	    std::copy_n(sourceVec_j.begin(),numFuelSpecies,mdot_j.begin());
	    std::copy_n(sourceVec_jp1.begin(),numFuelSpecies,mdot_jp1.begin());
	    //doublereal q_j = sourceVec_j[numFuelSpecies];
	    //doublereal q_jp1 = sourceVec_jp1[numFuelSpecies];

	    //doublereal sum_mdot_j = std::accumulate(mdot_j.begin(), mdot_j.end(),0.0);
	    //doublereal sum_mdot_jp1 = std::accumulate(mdot_jp1.begin(), mdot_jp1.end(),0.0);

	    doublereal evap_term_j = sourceVec_j[numFuelSpecies + 1];
            doublereal cl_d_j = sourceVec_j[numFuelSpecies + 2];

            //------------------------------------------------
            //    Number density equation
            //------------------------------------------------
            evalNumberDensity(j, x, rsd, diag, rdt);

            //------------------------------------------------
            //    Mass equation
            //
            //    dm_l,i/dt + v_l dm_l,i/dz = -mdot_i
            //-------------------------------------------------
            for (size_t i = 0; i < numFuelSpecies; i++) {
                 rsd[index(c_offset_Y+m_nsp+c_offset_ml+i,j)] = -vl(x,j) * dmlkdz(x,i,j) - 
                      rdt * (mlk(x,i,j) - mlk_prev(i,j)) - mdot_j[i] + av_mlk(x,i,j);
                 diag[index(c_offset_Y+m_nsp+c_offset_ml+i, j)] = 1;
	    }

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    m_l dU_l/dt + m_l v_l dU_l/dz + m_l U_l^2
            //       = f_r/r
            //-------------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] = -vl(x,j) * dUldz(x,j) -
                Ul(x,j)*Ul(x,j) - rdt*(Ul(x,j)-Ul_prev(j)) + av_Ul(x,j);
            if (ml(x,j)>sqrt(numeric_limits<double>::min())) {
                rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] += Fr(x,j)/ml(x,j);
            }
            diag[index(c_offset_Y+m_nsp+c_offset_Ul, j)] = 1;

            //------------------------------------------------
            //    Axial momentum equation
            //
            //    m_l dv_l/dt + m_l v_l dv_l/dz = f_z
            //-------------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] = -vl(x,j)*dvldz(x,j) - 
                    rdt*(vl(x,j)-vl_prev(j)) + av_vl(x,j);
            if (ml(x,j)>sqrt(numeric_limits<double>::min())) {
                rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] += fz(x,j)/ml(x,j);
            }
            diag[index(c_offset_Y+m_nsp+c_offset_vl, j)] = 1;

            //-----------------------------------------------
            //    energy equation
            //
            //    m_l c_p_l dT_l/dt + m_l c_p_l v_l dT_l/dz
            //    = mdot_l (q - L) 
            //-----------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] = -vl(x,j)*dTldz(x,j) - 
                    rdt * (Tl(x,j) - Tl_prev(j)) + av_Tl(x,j);
            if (ml(x,j) > 0.0 && cl_d_j > 0.0) 
                rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] += evap_term_j/ml(x,j) / cl_d_j;
            else 
                rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] += 0.0;
 
            diag[index(c_offset_Y+m_nsp+c_offset_Tl, j)] = 1;

        }
    } 
    

}

void SprayFlame::evalNumberDensity(size_t j, doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
     //----------------------------------------------
     //    Number density equation
     //
     //    Note that this propagates the liquid mass flow rate information to the right
     //    (j+1 -> j) from the value specified at the left boundary.
     //
     //    d(n_l v_l)/dz + 2n_l U_l = 0
     //------------------------------------------------
     rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] = -vl(x,j) * dnldz(x,j) -
         nl(x,j) * (vl(x,j) - vl(x,j-1))/m_dz[j-1] -
         2.0 * nl_Ul(x,j) - rdt * (nl(x,j) - nl_prev(j)) + av_nl(x,j);
     // rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] =
     //     -(nl_vl(x,j) - nl_vl(x,j-1))/m_dz[j-1]
     //     -(nl_Ul(x,j) + nl_Ul(x,j-1))
     //     -rdt * (nl(x,j) - nl_prev(j));

     diag[index(c_offset_Y+m_nsp+c_offset_nl, j)] = 1;
}


// class is not const because transport and thermo is updated 
std::vector<double> SprayFlame::source(const doublereal* x, size_t j) {

    // Source term is mass(# species) + temperature + (non latent droplet term) + (droplet specific heat)
    // Last two quantities required for liquid phase equation
    std::vector<double> sourceTerm(numFuelSpecies + 3, 0.0);
    
    // Return zero if last element
    if (j == m_points - 1)
         return sourceTerm;    

    std::vector<double>::iterator mdotEnd = (sourceTerm.end() - 3);

    // Start the drama only if droplet exists
    if (ml(x,j)> 0.0 && dl(x,j) > 0.0) {

    // Reset boiling flag
    bool isBoiling = false;

    // Get ambient fractions of each liquid component
    std::vector<doublereal> Yspec(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)
	 Yspec[i] = Y(x,c_offset_fuel[i],j);


    // Get mass fractions of each liquid component
    std::vector<doublereal> Ycomp_l(numFuelSpecies);
    for (auto& f : Ycomp_l) f = mlk(x,&f-&Ycomp_l[0],j);
    doublereal sum_Ycomp_l = std::accumulate(Ycomp_l.begin(),Ycomp_l.end(), 0.0);
    if (sum_Ycomp_l != 0) {
      for (auto& f : Ycomp_l) f /= sum_Ycomp_l;
    }

    // Get mixture molecular weight
    std::vector<doublereal> MWVec(numFuelSpecies);
    MW(MWVec.data());
    // Convert to Cantera units
    for (auto &f : MWVec)
         f *= 1e3;

    // Get mole fraction
    std::vector<doublereal> Xcomp_l(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)
	 Xcomp_l[i] = Ycomp_l[i]/MWVec[i];
    doublereal sum_Xcomp_l = std::accumulate(Xcomp_l.begin(), Xcomp_l.end(), 0.0);
    if (sum_Ycomp_l != 0) {
      for (auto& f : Xcomp_l) f /= sum_Xcomp_l;
    }

    // Saturation pressure
    std::vector<doublereal> Psat_comp(numFuelSpecies);
    double tl = Tl(x,j);
    PSat(Psat_comp.data(),&tl);

    // Raoult's law
    std::vector<doublereal> Xsat_comp(numFuelSpecies);
    if (m_evapModel == "distillation") {
	 for (size_t i = 0; i < numFuelSpecies; i++)
	      Xsat_comp[i] = Xcomp_l[i]*Psat_comp[i]/m_press;
    }
    else if (m_evapModel == "diffusion") {
	 for (size_t i = 0; i < numFuelSpecies; i++)
	      Xsat_comp[i] = std::min(Xcomp_l[i], Xcomp_l[i]*Psat_comp[i]/m_press);
    }

    // Check if in boiling regime
    doublereal sum_Xsat_comp = std::accumulate(Xsat_comp.begin(), Xsat_comp.end(), 0.0);
    if (sum_Xsat_comp >= 1.0) {
	 isBoiling = true;
         sum_Xsat_comp = 1.0;
	 for (size_t i = 0; i < numFuelSpecies; i++)
	      Xsat_comp[i] = Xsat_comp[i]/sum_Xsat_comp;
    }     
    
    // Air composition
    doublereal X_air = 1.0 - sum_Xsat_comp;

    // Mass fraction at droplet surface
    std::vector<doublereal> Ysat_comp(numFuelSpecies);
    doublereal MWsat_comp = std::inner_product(Xsat_comp.begin(),Xsat_comp.end(),MWVec.begin(),0.0) + X_air*m_wtm[j];
    for (size_t i = 0; i < numFuelSpecies; i++)
	 Ysat_comp[i] = Xsat_comp[i]*MWVec[i]/MWsat_comp;
    //doublereal Y_air = X_air * m_wtm[j] / MWsat_comp;

    // Average specific heat capacity
    std::vector<doublereal> c_lVec(numFuelSpecies);
    c_l(c_lVec.data(),&tl);
    double cl_d = std::inner_product(Ycomp_l.begin(),Ycomp_l.end(),c_lVec.begin(),0.0);
    
    // 1/3rd rule
    double Tref = (2.0/3.0)*Tl(x,j) + (1.0/3.0)*T(x,j);

    // Reference fractions
    std::vector<doublereal> Yref_comp(numFuelSpecies);
    // Air fraction in ambient
    doublereal Yref_remain = 0.0;

    // Evaluate reference mass fraction
    for (size_t i = 0; i < numFuelSpecies; i++) {
	Yref_comp[i] = (2.0/3.0)*Ysat_comp[i] + (1.0/3.0)*Yspec[i];
	Yref_remain += Yspec[i];      
    }
    //doublereal Yref_air = (2.0/3.0)*Y_air + (1.0/3.0)*(1.0 - Yref_remain);

    // Evaluate reference mole fraction
    doublereal MW_ref = 0.0;
    std::vector<doublereal> Xref_comp(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)
	 MW_ref += Yref_comp[i]/MWVec[i];
    MW_ref = 1.0e-3/MW_ref;
    for (size_t i = 0; i < numFuelSpecies; i++)
	 Xref_comp[i] = Yref_comp[i]*MW_ref/MWVec[i];

    // Evaluate reference density	
    doublereal rho_ref = m_press*MW_ref/(R_u * Tref);
	    
    // Evaluate mixture viscosity, specific heat, thermal conductivity
    doublereal mu_ref = m_visc[j];
    doublereal cp_ref = m_cp[j];
    doublereal lambda_ref = m_tcon[j];

    // Get diffusion coefficients
    std::vector<double> D_ref(numFuelSpecies);
    double curr_pressure = m_press;
    D(D_ref.data(),&curr_pressure,&Tref);

    // Non-dimensional numbers
    doublereal relative_velocity = sqrt(pow(Ul(x,j) - u(x,j),2.0) + pow(vl(x,j) - V(x,j),2.0));
    doublereal Red_sl = dl(x,j)*m_rho[j]*relative_velocity/m_visc[j];
    doublereal tau_d = rhol(x,j)*pow(dl(x,j),2.0)/(18.0*m_visc[j]);

    doublereal Pr_ref_g = cp_ref*mu_ref/lambda_ref;
    std::vector<double> Sc_ref_g(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)
	 Sc_ref_g[i]  = mu_ref/(rho_ref*D_ref[i]);

    // Mixture-averaged Schmidt
    doublereal D_ref_bar = 0.0;
    for (size_t i = 0; i < numFuelSpecies; i++)
	D_ref_bar += Xref_comp[i]/D_ref[i];
    D_ref_bar = D_ref_bar/std::accumulate(Xref_comp.begin(),Xref_comp.end(),0.0);
    D_ref_bar = 1.0/D_ref_bar;
    doublereal Sc_ref_g_avg = mu_ref/(rho_ref*D_ref_bar);

    // Sherwood and Nusselt numbers
    doublereal Nu_ref_g = 2.0 + 0.552 * pow(Red_sl,0.5)*pow(std::max(0.0,Pr_ref_g),1.0/3.0);
    std::vector<double> Sh_ref_g(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)	
	 Sh_ref_g[i] = 2.0 + 0.552 * pow(Red_sl,0.5)*pow(std::max(0.0,Sc_ref_g[i]),1.0/3.0);

    // Mixture-averaged Sherwood
    doublereal Sh_ref_g_avg = 2.0 + 0.552 * pow(Red_sl,0.5)*pow(std::max(0.0,Sc_ref_g_avg),1.0/3.0);
    
    // Mass Spalding number
    doublereal sum_Yfuels_gas = std::accumulate(Yspec.begin(),Yspec.end(),0.0);
    std::vector<double> eta_k(numFuelSpecies);
    for (size_t i = 0; i < numFuelSpecies; i++)
	eta_k[i] = ( D_ref_bar * Sh_ref_g_avg ) / ( D_ref[i] * Sh_ref_g[i] );

    // Calculate mass split
    std::vector<double> eps_k(numFuelSpecies);
    std::vector<double> Bm_k(numFuelSpecies);
    doublereal Bm = 0.0;
    if (isBoiling == true) {
	if (m_evapModel == "distillation") {
	// BMi goes to infinity when droplet boils. Selective suppression effects absent
	// This is enforced unfortunately, but for a very short period - not much effect on total ignition
	     eps_k = Ysat_comp;
	}
	else if (m_evapModel == "diffusion") 
	     eps_k = Ycomp_l;
    }

    else {
	 doublereal sum_Ysat_comp = std::accumulate(Ysat_comp.begin(), Ysat_comp.end(), 0.0);
	 Bm = ( sum_Ysat_comp - sum_Yfuels_gas ) / ( 1.0 - sum_Ysat_comp);

	 // Handle the saturated ambient case
	 // <= -1.0 used to account for roundoff
	 if (Bm <= -1.0)
	      std::fill(Bm_k.begin(),Bm_k.end(),-1.0);
	 else {
	     for (size_t i = 0; i < numFuelSpecies; i++)
		  Bm_k[i] = pow(1.0 + Bm,eta_k[i]) - 1.0;
	 }

    // Calculate BM_k based on individual Spalding numbers for the distillation model
    if (m_evapModel == "distillation") {
	 for (size_t i = 0; i < numFuelSpecies; i++)
              if (Bm_k[i] == 0)
                   eps_k[i] = Ysat_comp[i];
              else 
	           eps_k[i] = Ysat_comp[i] + ( Ysat_comp[i] - Yspec[i] ) / Bm_k[i];
    }
    // For the diffusion model, you only get the entire source term from the Spalding numbers
    // The split is however decided only based on the liquid composition
    else if (m_evapModel == "diffusion") 
	 eps_k = Ycomp_l;

    // Normalize eps_k to 1 because sum must equal mdot
    // Correction not applied when eps_k is negative because
    // that corresponds to condensation on droplets.
    // Whether the correlation holds then is an open question
    if (std::all_of(eps_k.begin(),eps_k.end(),[](double i) {return i > 0.0;})) {
      doublereal sum_eps_k = std::accumulate(eps_k.begin(),eps_k.end(),0.0);
      for (auto& f : eps_k) 
	   f /= sum_eps_k;
    }

    }

    // Compute mass source term
    std::vector<double> LvVec(numFuelSpecies);
    Hv(LvVec.data(),&tl);
    
    // Mass source terms must be positive 
    if (isBoiling == true) {
        doublereal latentTerm = std::inner_product(LvVec.begin(),LvVec.end(),eps_k.begin(),0.0);
	// This term is obtained by forcing dTdt is zero
	for (size_t i = 0; i < numFuelSpecies; i++)
            sourceTerm[i] = - ml(x,j) * cp_ref * Nu_ref_g * (1.0/tau_d) * (Tl(x,j) - T(x,j)) * eps_k[i]/(3.0 * Pr_ref_g * latentTerm);
    }
    else {
	 if (Bm > -1.0) {
	     for (size_t i = 0; i < numFuelSpecies; i++)
		  sourceTerm[i] =  Sh_ref_g_avg * ml(x,j) /( 3.0 * Sc_ref_g_avg * tau_d) * log(1.0 + Bm) * eps_k[i];
	 }
    }
    
    // Compute correction coefficient for temperature source term
    doublereal beta = 0.0;
    if ( ml(x,j) > 0.0  && dl(x,j) > 0.0 )
	 beta = -1.5 * Pr_ref_g * std::accumulate(sourceTerm.begin(),mdotEnd,0.0) / (ml(x,j) / tau_d);
    doublereal f2 = 0.0;
    if ( beta != 0.0 ) f2 = beta / ( exp(beta) - 1.0 );

    // Latent heat of vaporization
    doublereal Evap_latent_heat = std::inner_product(sourceTerm.begin(),mdotEnd,LvVec.begin(),0.0);
    // Reverse latent heat sign due to mdot reversal
    Evap_latent_heat = -Evap_latent_heat;

    // Temperature source term
    if (isBoiling == true) { 
	 sourceTerm[numFuelSpecies] = 0.0;
	 sourceTerm[numFuelSpecies + 1] = 0.0;
    }
    else {
         if (ml(x,j) > 0.0 && dl(x,j) > 0.0) {
     	   sourceTerm[numFuelSpecies] = ml(x,j) * Nu_ref_g * cp_ref * f2 * ( T(x,j) - Tl(x,j) ) / ( 3.0 * Pr_ref_g * tau_d )
				        + Evap_latent_heat;
    }
	 sourceTerm[numFuelSpecies + 1] = sourceTerm[numFuelSpecies] - Evap_latent_heat; 	
    }

    // Set specific heat capacity. Used in flamelet equation in eval
    sourceTerm[numFuelSpecies + 2] = cl_d;

    }

    return sourceTerm;
}


std::vector<doublereal> SprayFlame::getMW() {
    std::vector<doublereal> MWVec(numFuelSpecies);
    MW(MWVec.data());
    // Convert to Cantera units
    for (auto &f : MWVec)
         f *= 1e3;
    return MWVec;
}

std::vector<doublereal> SprayFlame::getRhoL(double T) {
    std::vector<doublereal> rhol(numFuelSpecies);
    rho_l(rhol.data(),&T);
    return rhol;
}

doublereal SprayFlame::rhol(const doublereal* x, size_t j) const {
    // Calculate average density
    std::vector<doublereal> rhoL(numFuelSpecies);
    doublereal temp = Tl(x,j);
    rho_l(rhoL.data(),&temp);
    double ret = 0.0;
    for (size_t i = 0; i < numFuelSpecies; i++)
      ret += mlk(x,i,j)/rhoL[i];
    if (ret <= 0.) 
      return 0.0;
    else
      return ml(x,j)/ret; 
}

void SprayFlame::evalRightBoundaryLiquid(doublereal* x, doublereal* rsd,
                                         integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;

    rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] = nl(x,j) - nl(x,j-1);
    // rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] =
    //      -(nl_vl(x,j) - nl_vl(x,j-1))/m_dz[j-1]
    //      -(nl_Ul(x,j) + nl_Ul(x,j-1));

    diag[index(c_offset_Y+m_nsp+c_offset_nl, j)] = 0;

    // Neumann boundary condition for Ul, vl, Tl and ml
    rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] = Ul(x,j) - Ul(x,j-1);
    rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] = vl(x,j) - vl(x,j-1);
    rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] = Tl(x,j) - Tl(x,j-1);
    for (size_t i = 0; i < numFuelSpecies; i++)
         rsd[index(c_offset_Y+m_nsp+c_offset_ml+i,j)] = mlk(x,i,j) - mlk(x,i,j-1);
    diag[index(c_offset_Y+m_nsp+c_offset_Ul, j)] = 0;
    diag[index(c_offset_Y+m_nsp+c_offset_vl, j)] = 0;
    diag[index(c_offset_Y+m_nsp+c_offset_Tl, j)] = 0;
    for (size_t i = 0; i < numFuelSpecies; i++)
         diag[index(c_offset_Y+m_nsp+c_offset_ml+i, j)] = 0;
}

string SprayFlame::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "u";
    case 1:
        return "V";
    case 2:
        return "T";
    case 3:
        return "lambda";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else if (n >= (c_offset_Y + m_nsp) && n < (c_offset_Y + m_nsp + c_offset_ml + numFuelSpecies)) {
            size_t idx = n - c_offset_Y - m_nsp;
	    if (idx == 0)
	        return "Ul";
	    else if (idx == 1)
		return "vl";
	    else if (idx == 2) 
                return "Tl";
	    else if (idx == 3) 
                return "nl";
	    else
                return "ml_" + m_palette[idx - c_offset_ml];
        } else 
            return "<unknown>";
    }
}

size_t SprayFlame::componentIndex(const std::string& name) const
{
    if (name=="u") {
        return 0;
    } else if (name=="V") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else if (name=="Ul") {
        return c_offset_Y + m_nsp + c_offset_Ul;
    } else if (name=="vl") {
        return c_offset_Y + m_nsp + c_offset_vl;
    } else if (name=="Tl") {
        return c_offset_Y + m_nsp + c_offset_Tl;
    } else if (name=="ml") {
        return c_offset_Y + m_nsp + c_offset_ml;
    } else if (name=="nl") {
        return c_offset_Y + m_nsp + c_offset_nl;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
    }
    return npos;
}

XML_Node& SprayFlame::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = StFlow::save(o, sol);

    XML_Node& gv = flow.addChild("liq_grid_data");
    vector_fp x(soln.nColumns());

    soln.getRow(componentIndex("Ul"), x.data());
    addFloatArray(gv,"Ul",x.size(),x.data(),"m/s","liq. r velocity");

    soln.getRow(componentIndex("vl"), x.data());
    addFloatArray(gv,"vl",x.size(),x.data(),"m/s","liq. x velocity");

    soln.getRow(componentIndex("Tl"), x.data());
    addFloatArray(gv,"Tl",x.size(),x.data(),"K","liq. temperature");

    soln.getRow(componentIndex("nl"), x.data());
    addFloatArray(gv,"nl",x.size(),x.data(),"/m^3","number density");

    for (size_t i = 0; i < numFuelSpecies; i++) {
        soln.getRow(componentIndex("ml") + i, x.data());
        addFloatArray(gv,"ml",x.size(),x.data(),"kg","droplet mass");
    }

    return flow;
}

} // namespace
