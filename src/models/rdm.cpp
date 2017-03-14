#include "cantera/models/RDM.h"
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{

RDM::RDM(StFlow* const sf_, size_t delta) : 
     ModelGeneric(sf_),m_delta(delta)
{
    opt_interp_s = "spline";
    spline_bc_s = "parabolic-run-out";
    filter_type = "Tophat";
    alpha_amp = 1.0;
    use_constrained_qp = false;
    use_Landweber = false;

    update_grid();

}

void RDM::getWdot(const doublereal* x)
{
    size_t nsp = sf->phase().nSpecies();
    size_t nvar = c_offset_Y + nsp;

    doublereal* xnew = new doublereal[nvar*m_nz];

    // Deconvolving the variables
    for (size_t k=0; k<nvar; k++) {
        // Assemble the scalar k into one vector
        vector_fp sc;
        for (size_t i=0; i<m_npts; i++) {
            sc.push_back(x[k+i*nvar]);
        }
        // Interpolation
        vector_fp sc_bar = interpolation(sc);
        // Deconvolution
        vector_fp sc_dcv = deconvolution(sc_bar,k);
        // Subgrid reconstruction
        subgrid_reconstruction(sc_bar,sc_dcv);
        // Assemble the scalar k back into the mixed vector
        for (size_t i=0; i<m_nz; i++) {
            // xnew[k+i*nvar] = sc_dcv[i];
            xnew[k+i*nvar] = sc_bar[i];
        }
    }
    // Compute the source term
    for (size_t j=0; j<m_nz; j++) {
        sf->setGas(xnew,j);
        sf->kinetics().getNetProductionRates(&m_wdot_fine(0,j));
    }
    // Filter the source term and down-scale it to the original size
    for (size_t k=0; k<nsp; k++) {
        vector_fp src;
        for (size_t j=0; j<m_nz; j++) {
            src.push_back(m_wdot_fine(k,j));
        }
        vector_fp src_buf = filtering(src);
        // /* for development purpose only
        for (size_t j=0; j<m_nz; j++) {
            m_wdot_fine(k,j) = src_buf[j];
        }
        // */
        vector_fp src_buf2 = down_sampling(src_buf);
        for (size_t j=0; j<m_npts; j++) {
            // m_wdot(k,j) = src_buf2[j];
            doublereal sc_ = x[c_offset_Y+k+j*nvar];
            m_wdot(k,j) = ((src_buf2[j]>0.0 && std::abs(sc_-1.0)<1.0e-10) || 
                           (src_buf2[j]<0.0 && std::abs(sc_-0.0)<1.0e-10))
                          ? 0.0
                          : src_buf2[j];
        }
    }

    delete[] xnew;

}

void RDM::update_grid()
{
    vector_fp z = sf->grid();
    size_t refine_factor = (m_delta+1)/2;

    // fine mesh for interpolation
    m_npts = z.size();
    m_nz = (m_npts-1)*refine_factor+1;
    m_z.resize(m_nz);
    for (size_t i=0; i<m_npts-1; i++) {
        doublereal dz = (z[i+1]-z[i])/(doublereal)refine_factor;
        for (size_t j=0; j<refine_factor; j++) {
            m_z[j+i*refine_factor] = z[i] + (doublereal)j*dz; 
        }
    }
    m_z[m_nz-1] = z[m_npts-1];

    // interpolation matrix
    interp_factory(z);

    // filter and deconvolution operators
    operator_init();

    // solution arrays
    m_wdot.resize(sf->phase().nSpecies(),m_npts,0.0); m_wdot*=0.0;
    m_wdot_fine.resize(sf->phase().nSpecies(),m_nz,0.0); m_wdot_fine*=0.0;

}

void RDM::interp_factory(const vector_fp& z)
{
    if (opt_interp_s.compare("spline")==0) {
        opt_interp = 1;
        if (spline_bc_s.compare("free-run-out")==0) {
            spline_bc = 1;
        } else if(spline_bc_s.compare("parabolic-run-out")==0) {
            spline_bc = 2;
        } else if(spline_bc_s.compare("periodic")==0) {
            spline_bc = 3;
        } else {
            spline_bc = 1;
        }
    } else if (opt_interp_s.compare("Lagrange-polynomials")==0) {
        opt_interp = 2;
    } else if (opt_interp_s.compare("ENO")==0) {
        opt_interp = 3;
    }
    m_interp.resize(m_nz,m_npts,0.0); m_interp*=0.0;

    switch (opt_interp) {
    case 1: {
                // Spline interpolation
                vector_fp dz;
                dz.resize(m_npts-1);
                for (size_t i=0; i<m_npts-1; i++) {
                    dz[i] = z[i+1]-z[i];
                }
                
                DenseMatrix G(m_npts,m_npts,0.0);
                DenseMatrix rhs(m_npts,m_npts,0.0);
                G(0,0) = 1.0;
                G(m_npts-1,m_npts-1) = 1.0;
                for (size_t i = 1; i<m_npts-1; i++) {
                    G(i,i-1) = dz[i-1]/6.0;
                    G(i,i  ) = (dz[i-1]+dz[i])/3.0;
                    G(i,i+1) = dz[i]/6.0;
                    rhs(i,i-1) = 1.0/dz[i-1];
                    rhs(i,i  ) = -(1.0/dz[i-1]+1.0/dz[i]);
                    rhs(i,i+1) = 1.0/dz[i];
                }

                switch (spline_bc) {
                case 1: {
                            // free-run-out
                            // bc remain at 0
                            break;
                        }
                case 2: {
                            // parabolic-run-out
                            G(0,0) = G(1,0)+G(1,1);
                            G(0,2) = G(1,2);
                            doublereal* buf = new doublereal[m_npts];
                            rhs.getRow(1,buf);
                            rhs.setRow(0,buf);
                            G(m_npts-1,m_npts-1) = G(m_npts-2,m_npts-1)+G(m_npts-2,m_npts-2);
                            G(m_npts-1,m_npts-3) = G(m_npts-2,m_npts-3);
                            rhs.getRow(m_npts-2,buf);
                            rhs.setRow(m_npts-1,buf);
                            delete[] buf;
                            break;
                        }
                case 3: {
                            // periodic
                            break;
                        }
                default: break;
                }

                DenseMatrix C1(m_nz,m_npts,0.0);
                DenseMatrix C2(m_nz,m_npts,0.0);
                size_t i_src = 0;
                for (size_t i=0; i<m_nz; i++) {
                    doublereal z_ = m_z[i];
                    if (z_>z[i_src+1]) {
                        i_src++;
                    }
                    doublereal z_i  = z[i_src];
                    doublereal z_ip = z[i_src+1];
                    doublereal dz_i = dz[i_src];
                    C1(i,i_src  ) = (std::pow(z_ip-z_,3.0)/dz_i-dz_i*(z_ip-z_))/6.0;
                    C1(i,i_src+1) = (std::pow(z_-z_i ,3.0)/dz_i-dz_i*(z_-z_i ))/6.0;
                    C2(i,i_src  ) = (z_ip-z_)/dz_i;
                    C2(i,i_src+1) = (z_-z_i )/dz_i;
                }

                solve(G,rhs);
                C1.mult(rhs,m_interp);
                m_interp+=C2;

                break;
            }
    case 2: {
                order = 1;
                size_t idx = 0;
                for (size_t i = 0; i < m_nz; i++) {
                    if (m_z[i] > z[idx+1]) {
                        ++idx;
                    }
                    if (idx > m_npts-order-1) {
                        vector_fp coeff = Lagrange_polynomial(z,m_z[i],m_npts-order-1,m_npts-1);
                        auto it = coeff.begin();
                        for (size_t j = m_npts-order-1; j < m_npts; j++) {
                            m_interp(i,j) = (*it);
                            ++it;
                        }

                    } else {
                        vector_fp coeff = Lagrange_polynomial(z,m_z[i],idx,idx+order);
                        auto it = coeff.begin();
                        for (size_t j = idx; j <= idx+order; j++) {
                            m_interp(i,j) = (*it);
                            ++it;
                        }
                    }
                }
                break;
            }
    case 3: {
                order = 5;
                interp_pool.resize(order);
                for (auto it = interp_pool.begin(); it != interp_pool.end(); ++it) {
                    (*it).resize(m_nz,m_npts);
                    (*it) *= 0.0;
                }
                stencil_idx.resize(m_nz,order);
                stencil_idx *= 0.0;
                
                for (size_t s = 0; s < order; s++) {
                    size_t idx = 0;
                    for (size_t i = 0; i < m_nz; i++) {
                        if (m_z[i] > z[idx+1]) {
                            ++idx;
                        }
                        // left end index of stencil w.r.t current idx
                        int l = idx - s;
                        // right end index of stencil w.r.t current idx
                        int r = l + order;
                        // special case for stencil near the boundary
                        if (l < 0) {
                            vector_fp coeff = Lagrange_polynomial(z,m_z[i],0,order);
                            auto it = coeff.begin();
                            for (size_t j = 0; j <= order; j++) {
                                interp_pool[s](i,j) = (*it);
                                ++it;
                            }
                            
                        } else if (r >= m_npts) {
                            vector_fp coeff = Lagrange_polynomial(z,m_z[i],m_npts-order-1,m_npts-1);
                            auto it = coeff.begin();
                            for (size_t j = m_npts-order-1; j < m_npts; j++) {
                                interp_pool[s](i,j) = (*it);
                                ++it;
                            }

                        } else {
                            vector_fp coeff = Lagrange_polynomial(z,m_z[i],l,r);
                            auto it = coeff.begin();
                            for (size_t j = l; j <= r; j++) {
                                interp_pool[s](i,j) = (*it);
                                ++it;
                            }

                        }
                        stencil_idx(i,s) = std::min(std::max(l,0),(int)(m_npts-order-1));
                    }
                }
                    
                break;
            }
    default: {
                 break;
             }
    }
    /* for debug
    std::ofstream debug_file;
    debug_file.open("mat_debug_interp.dat");
    for (size_t i=0; i<m_interp.nRows(); i++) {
        for (size_t j=0; j<m_interp.nColumns(); j++) {
            debug_file << m_interp(i,j) << ",";
        }
        debug_file << std::endl;
    }
    debug_file.close();
    */

}

void RDM::operator_init()
{
    m_G.resize(m_nz,m_nz,0.0); m_G*=0.0;
    m_Q.resize(m_nz,m_nz,0.0); m_Q*=0.0;
    m_SG.resize(m_nz,m_nz,0.0); m_SG*=0.0;

    // filter operator
    vector_fp m_G_;
    int c = m_nz/2;
    doublereal zc = m_z[c];
    doublereal G_cum = 0.0;
    for (size_t j=0; j<m_nz; j++) {
        doublereal g_;
        if (filter_type.compare("Tophat")==0) {
            g_ = std::abs((int)j-c)<=m_delta/2? 1.0 : 0.0;
        } else if (filter_type.compare("Differential")==0) {
            doublereal delta_ = (doublereal)(m_delta-1)*(m_z[1]-m_z[0]);
            g_ = std::exp(-std::abs(m_z[j]-zc)/delta_)/(2*delta_);
        } else {
            g_ = 1.0;
        }
        G_cum += g_;
        m_G_.push_back(g_);
    }
    std::transform(m_G_.begin(),m_G_.end(),m_G_.begin(),std::bind1st(std::multiplies<doublereal>(),1.0/G_cum));;

    for (size_t i=0; i<m_nz; i++) {
        int l_ = i-m_nz/2;
        int r_ = i+m_nz/2-(m_nz+1)%2;
        G_cum = 0.0;
        for (size_t j=std::max(0,l_); j<=std::min((int)(m_nz)-1,r_); j++) {
            m_G(i,j) = m_G_[std::max(0,c-(int)i+(int)j)];
            G_cum += m_G(i,j);
        }
        if (l_<0) {
            m_G(i,0) += (1.0-G_cum);
        } else if (r_>m_nz-1) {
            m_G(i,m_nz-1) += (1.0-G_cum);
        }
    }

    // deconvolution operator
    DenseMatrix I(m_nz,m_nz,0.0);
    // filter and zeroth order regularization
    DenseMatrix GT(m_nz,m_nz,0.0);
    for (size_t i=0; i<m_nz; i++) {
        I(i,i) = 1.0;
        for (size_t j=0; j<m_nz; j++) {
            GT(j,i) = m_G(i,j);
        }
    }
    DenseMatrix G2(m_nz,m_nz,0.0);
    GT.mult(m_G,G2);
    alpha = 0.0;
    for (size_t i=0; i<m_nz; i++) {
        alpha += G2(i,i);
    }
    alpha /= (doublereal)m_nz;
    alpha *= alpha_amp;
    // second order derivative regularization
    DenseMatrix L(m_nz-2,m_nz,0.0);
    DenseMatrix LT(m_nz,m_nz-2,0.0);
    double dzmin = m_z[1]-m_z[0];
    for (size_t i = 0; i < m_nz-2; i++) {
        double dz_i  = m_z[i+1] - m_z[i];
        double dz_ip = m_z[i+2] - m_z[i+1];
        L(i,i  ) =  2.0/dz_i /(dz_i+dz_ip);
        L(i,i+1) = -2.0/dz_i / dz_ip;
        L(i,i+2) =  2.0/dz_ip/(dz_i+dz_ip);
        LT(i  ,i) = L(i,i  );
        LT(i+1,i) = L(i,i+1);
        LT(i+2,i) = L(i,i+2);
        dzmin = std::min(dzmin,dz_ip);
    }
    L *= std::pow(dzmin,2.0);
    LT *= std::pow(dzmin,2.0);
    DenseMatrix L2(m_nz,m_nz,0.0);
    LT.mult(L,L2);
    L2 *= alpha;
    // compute the deconvolution operator
    m_Q = I;
    m_Q *= alpha;
    m_Q += GT;
    DenseMatrix buf(m_nz,m_nz,0.0);
    buf = I;
    buf *= alpha;
    buf += G2;
    buf += L2;
    solve(buf,m_Q);

    // subgrid reconstruction operator
    m_Q.mult(m_G,m_SG);
    m_SG *= -1.0;
    m_SG += I;

    // Initialization for matrices used in constained optimization
    qp_D = new doublereal*[m_nz];
    qp_C = new doublereal*[m_nz];
    for (size_t i = 0; i < m_nz; i++) {
        qp_D[i] = new doublereal[m_nz];
        qp_C[i] = new doublereal[m_nz];
        for (size_t j = 0; j < m_nz; j++) {
            qp_D[i][j] = ( G2(i,j) + alpha*I(i,j))*2.0;
            qp_C[i][j] =-(m_G(i,j) + alpha*I(i,j))*2.0;
        }
    }
    qp_A = new doublereal*[m_nz];
    for (size_t j = 0; j < m_nz; j++) {
        qp_A[j] = new doublereal[1];
        qp_A[j][0] = 1.0;
    }
    qp_l  = new doublereal[m_nz];
    qp_u  = new doublereal[m_nz];
    qp_fl = new bool[m_nz];
    qp_fu = new bool[m_nz];
    for (size_t j = 0; j < m_nz; j++) {
        qp_l[j]  = 0.0;
        qp_u[j]  = 1.0;
        qp_fl[j] = true;
        qp_fu[j] = true;
    }

    // Initialization for Landweber deconvolution
    lw_niter = 5;
    lw_Greg.resize(m_nz,m_nz,0.0); lw_Greg *= 0.0;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G_tmp;
    G_tmp.resize(m_nz,m_nz);
    for (size_t i = 0; i < m_nz; i++) {
        doublereal lw_d = 0.0;
        for (size_t j = 0; j < m_nz; j++) {
            G_tmp(i,j) = m_G(i,j);
            lw_d += std::pow(m_G(i,j),2.0);
        }
        lw_d = 1.0/lw_d/(doublereal)m_nz;
        for (size_t j = 0; j < m_nz; j++) {
            lw_Greg(j,i) = m_G(j,i)*lw_d;
        }
    }
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G2_tmp;
    G2_tmp.resize(m_nz,m_nz);
    G2_tmp = G_tmp.transpose()*G_tmp;
    doublereal lw_w = 1.0/G2_tmp.lpNorm<2>();
    lw_Greg *= lw_w;

    /* for debug
    std::ofstream debug_file;
    debug_file.open("mat_debug_G.dat");
    for (size_t i=0; i<m_Q.nRows(); i++) {
        for (size_t j=0; j<m_Q.nColumns(); j++) {
            debug_file << m_G(i,j) << ",";
        }
        debug_file << std::endl;
    }
    debug_file.close();
    */
}

vector_fp RDM::mat_vec_multiplication(const DenseMatrix& A,const vector_fp& x)
{
    assert(A.nColumns()==x.size());
    double* x_ = new double[x.size()];
    for (size_t i=0; i<x.size(); i++) {
        x_[i] = x[i];
    }
    double* sol = new double[A.nRows()];
    multiply(A,x_,sol);
    vector_fp sol_;
    for (size_t i=0; i<A.nRows(); i++) {
        sol_.push_back(sol[i]);
    }
    delete[] x_;
    delete[] sol;
    return sol_;
}

void RDM::subgrid_reconstruction(const vector_fp& xbar, vector_fp& xdcv)
{
    assert(xbar.size()==xdcv.size());

    vector_fp xvar1,xvar2;

    // variance 1: compute from algebric model
    doublereal Cz = 0.8;
    for (vector_fp::const_iterator it = xbar.begin(); it != xbar.end(); ++it) {
        xvar1.push_back(std::pow(*it,2.0));
    }
    vector_fp buf1 = filtering(xvar1);
    vector_fp buf2 = filtering(xbar);
    vector_fp::iterator it1 = buf1.begin();
    vector_fp::iterator it2 = buf2.begin();
    for (vector_fp::iterator it = xvar1.begin(); it != xvar1.end(); ++it) {
        *it = Cz*std::max(*it1 - std::pow(*it2,2.0), 0.0);
        ++it1;
        ++it2;
    }
    
    // variance 2: compute from deconvolution results
    vector_fp::const_iterator it0 = xbar.begin();
    it2 = xdcv.begin();
    for (vector_fp::iterator it = xvar2.begin(); it != xvar2.end(); ++it) {
        *it = std::pow(*it0-*it2,2.0);
        ++it0;
        ++it2;
    }

    // compute the subgrid contribution
    buf1 = mat_vec_multiplication(m_SG,xbar); 
    it1 = xvar1.begin();
    it2 = xvar2.begin();
    vector_fp::iterator it3 = buf1.begin();
    for (vector_fp::iterator it = xdcv.begin(); it != xdcv.end(); ++it) {
        *it += std::sqrt((*it1)/(*it2))*(*it3);
        ++it1;
        ++it2;
        ++it3;
    }
       
}

vector_fp RDM::constrained_deconvolution(const vector_fp& x, size_t k)
{
    assert(x.size()==m_nz);

    // Update the constraints matrices
    doublereal* qp_c = new doublereal[m_nz];
    doublereal* qp_b = new doublereal[1];
    doublereal qp_c0 = 0.0;
    qp_b[0] = 0.0;
    for (size_t j = 0; j < m_nz; j++) {
        qp_c[j] = 0.0;
    }
    size_t i = 0;
    for (auto it = x.begin(); it != x.end(); ++it) {
        qp_c0 += std::pow(*it,2.0);
        qp_b[0] += (*it);
        for (size_t j = 0; j < m_nz; j++) {
            qp_c[j] += (*it)*qp_C[i][j];
        }
        i++;
    }
    qp_c0 *= (1.0+alpha);
    if (k<c_offset_Y) {
        for (size_t j = 0; j < m_nz; j++) {
            qp_fl[j] = false;
            qp_fu[j] = false;
        }
    } else {
        for (size_t j = 0; j < m_nz; j++) {
            qp_fl[j] = true;
            qp_fu[j] = true;
        }
    }

    // linear constraint is "="
    CGAL::Const_oneset_iterator<CGAL::Comparison_result> 
          r(    CGAL::EQUAL);                 	 
    // now construct the quadratic program; the first two parameters are
    // the number of variables and the number of constraints (rows of A)
    Program qp (m_nz, 1, qp_A, qp_b, r, qp_fl, qp_l, qp_fu, qp_u, qp_D, qp_c, qp_c0);
    // solve the program, using ET as the exact type
    Solution s = CGAL::solve_quadratic_program(qp, ET());
    // output solution
    // std::cout << s;
    vector_fp sol;
    for (auto it = s.variable_values_begin(); it != s.variable_values_end(); ++it) {
        sol.push_back(to_double(*it));
    }

    delete[] qp_c;
    delete[] qp_b;

    return sol;
}

vector_fp RDM::Landweber_deconvolution(const vector_fp& x, size_t k)
{
    vector_fp sol(x.begin(),x.end());

    for (size_t i = 0; i < lw_niter; i++) {
        vector_fp buf = mat_vec_multiplication(m_G,sol);
        auto it_x = x.begin();
        for (auto it = buf.begin(); it != buf.end(); ++it) {
            *it = (*it_x) - (*it);
            ++it_x;
        }
        buf = mat_vec_multiplication(lw_Greg,buf);
        it_x = buf.begin();
        for (auto it = sol.begin(); it != sol.end(); ++it) {
            if (k < c_offset_Y && k != c_offset_T) {
                *it = (*it) + (*it_x);
            } else if (k == c_offset_T) {
                *it = std::min(std::max( (*it) + (*it_x), 300.0), 2500.0);
            } else {
                *it = std::min(std::max( (*it) + (*it_x), 0.0), 1.0);
            }
            ++it_x;
        }
    }
    return sol;

}

void RDM::update_ENO_interp(const vector_fp& x)
{
    vector_fp v;
    for (auto it = x.begin(); it != x.end(); ++it) {
        v.push_back(*it);
    }
    for (size_t s = 0; s < order; s++) {
        for (size_t i = 0; i < v.size()-1; i++) {
            v[i] = (v[i+1]-v[i])/(sf->grid(i+s+1)-sf->grid(i));
        }
        v.resize(v.size()-1);
    }
    for (size_t i = 0; i < m_nz; i++) {
        vector_fp v_;
        for (size_t s = 0; s < order; s++) {
            v_.push_back(std::abs(v[stencil_idx(i,s)]));
        }
        size_t idx = std::distance(v_.begin(),std::min_element(v_.begin(),v_.end()));
        for (size_t j = 0; j < x.size(); j++) {
            m_interp(i,j) = interp_pool[idx](i,j);
        }
    }
       
}

} // Cantera

