#ifndef LATTICE_COUPLING_H
#define LATTICE_COUPLING_H

#include "custom_dynamics.h"
#include <cmath>
#include <palabos2D.h>
#include <palabos2D.hh>

#define PHI_GRAD_FIELD 21
#define PHI_LAPLACE_FIELD 22
using namespace plb;

template <typename T, template <typename U> class Descriptor>
class lattice_coupling
    : public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
  // constructor initialization (order matches member declaration)
  lattice_coupling(BlockLattice2D<T, Descriptor> &phi, T omega, T chi, T mu,
                   T a, T b, T epsilon, T c_bulk_k, T tau1, T tau2)
      : omega_(omega), chi_(chi), mu_(mu), a_(a), b_(b), epsilon_(epsilon),
        c_bulk_k_(c_bulk_k), tau1_(tau1), tau2_(tau2),
        dyna_(custom_dynamics<T, Descriptor>(omega, chi, mu)), phi_(phi) {}

  void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice1,
               BlockLattice2D<T, Descriptor> &lattice2) override {
    const T density_floor = (T)1e-12;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint jY = domain.y0; jY <= domain.y1; ++jY) {
        Cell<T, Descriptor> &cell1 = lattice1.get(iX, jY);
        Cell<T, Descriptor> &cell2 = lattice2.get(iX, jY);
        Cell<T, Descriptor> &cellPhi = phi_.get(iX, jY);

        // densities with floor
        T c1 = std::max(cell1.computeDensity(), density_floor);
        T c2 = std::max(cell2.computeDensity(), density_floor);
        T phi_val = std::max(cellPhi.computeDensity(), density_floor);

        // compute velocity and momentum for lattice1
        Array<T, Descriptor<T>::d> u1;
        cell1.computeVelocity(u1);               // fills u1
        Array<T, Descriptor<T>::d> j1 = c1 * u1; // momentum
        T jSqr1 = dot(j1, j1);

        // compute velocity and momentum for lattice2
        Array<T, Descriptor<T>::d> u2;
        cell2.computeVelocity(u2);
        Array<T, Descriptor<T>::d> j2 = c2 * u2;
        T jSqr2 = dot(j2, j2);

        // reaction / source terms
        T Sj1 = (T)0;
        T Sj2 = (T)0;

        if (phi_val >= (T)0.5) {
          T J1 = (1.0 / epsilon_) *
                 (c1 * (c1 - (T)1.0) - ((b_ * c2 * (c1 - a_)) / (c1 + a_)));
          T J2 = c1 - c2;
          Sj1 = J1;
          Sj2 = J2;
        } else {
          Sj1 = -(c1 - c_bulk_k_);
          Sj2 = -(c2 - c_bulk_k_);
        }

        // collision step: use computed j and j^2 for equilibrium
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
          T feq1 = dyna_.computeEquilibrium(iPop, c1, j1, jSqr1);
          T feq2 = dyna_.computeEquilibrium(iPop, c2, j2, jSqr2);

          T f1 = cell1[iPop];
          T f2 = cell2[iPop];

          cell1[iPop] = f1 - (f1 - feq1) / tau1_ + Sj1;
          cell2[iPop] = f2 - (f2 - feq2) / tau2_ + Sj2;
        }
      }
    }
  }

  lattice_coupling<T, Descriptor> *clone() const override {
    return new lattice_coupling<T, Descriptor>(
        phi_, omega_, chi_, mu_, a_, b_, epsilon_, c_bulk_k_, tau1_, tau2_);
  }

private:
  // members (declaration order matters for initialization order)
  T omega_, chi_, mu_;
  T a_, b_;
  T epsilon_;
  T c_bulk_k_;
  T tau1_, tau2_;
  custom_dynamics<T, Descriptor> dyna_;
  BlockLattice2D<T, Descriptor> &phi_;
};

// class CouplePhiMomentum
#define PHI_GRAD_FIELD 21
#define PHI_LAPLACE_FIELD 22

template <typename T, template <typename U> class Descriptor>
class PhiPcoupling2D
    : public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
  PhiPcoupling2D(T beta, T kappa) : beta_(beta), kappa_(kappa) {}

  // this method need to be overrided , 
  void process(Box2D domain, BlockLattice2D<T, Descriptor> &phiLattice,
               BlockLattice2D<T, Descriptor> &pLattice) override {

    const T cs2 = Descriptor<T>::cs2;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint iY = domain.y0; iY <= domain.y1; ++iY) {

        // ---- Step 1: Get phi, grad(phi), laplacian(phi) ----
        T phi = phiLattice.get(iX, iY).computeDensity();

        Array<T, 2> gradPhi =
            phiLattice.get(iX, iY).template getExternal<Array<T, 2>>(
                PHI_GRAD_FIELD);

        T lapPhi =
            phiLattice.get(iX, iY).template getExternal<T>(PHI_LAPLACE_FIELD);

        // ---- Step 2: Compute mu_phi and surface-tension force ----
        T muPhi =
            4.0 * beta_ * phi * (phi - 0.5) * (phi - 1.0) - kappa_ * lapPhi;

        Array<T, 2> Fs;
        Fs[0] = muPhi * gradPhi[0];
        Fs[1] = muPhi * gradPhi[1];

        // ---- Step 3: Apply Fs as a source to the P lattice ----
        T rho = pLattice.get(iX, iY).computeDensity();
        Array<T, 2> j;
        pLattice.get(iX, iY).computeMomentum(j);

        Array<T, 2> u;
        for (int d = 0; d < 2; ++d)
          u[d] = j[d] / (rho + 1e-18);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
          const Array<T, 2> &e_i = Descriptor<T>::c[iPop];
          T w_i = Descriptor<T>::t[iPop];

          T eu = dot(e_i, u);
          T uu = dot(u, u);
          T Gamma_i = w_i * (1.0 + eu / cs2 + 0.5 * (eu * eu) / (cs2 * cs2) -
                             0.5 * uu / cs2);

          // S_j = dt * Gamma_j * (e_j - u) · Fs
          Array<T, 2> diff = e_i - u;
          T Sj = Gamma_i * dot(diff, Fs);

          T f_old = pLattice.get(iX, iY)[iPop];
          T feq = pLattice.get(iX, iY).getDynamics().computeEquilibrium(
              iPop, rho, j, dot(j, j));
          T omega = 1.0 / pLattice.get(iX, iY).getDynamics().getTau();

          // collision
          pLattice.get(iX, iY)[iPop] = f_old - (f_old - feq) * omega + Sj;
        }
      }
    }
  }

  PhiPcoupling2D<T, Descriptor> *clone() const override {
    return new PhiPcoupling2D<T, Descriptor>(*this);
  }

private:
  T beta_, kappa_;
};

#endif
