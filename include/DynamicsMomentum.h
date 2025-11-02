
#ifndef DYNAMICS_MOMENTUM_H
#define DYNAMICS_MOMENTUM_H
#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

#define FORCE_FIELD 25 // external field index for Fs

template <typename T, template <typename U> class Descriptor>
class DynamicsMomentum : public plb::BGKdynamics<T, Descriptor> {
public:
  DynamicsMomentum(T omega) : plb::BGKdynamics<T, Descriptor>(omega) {}

  // ---- Collision step with force-corrected velocity ----
  void collide(Cell<T, Descriptor> &cell,
               BlockStatistics &statistics) override {

    // --- Step 1: Compute density and momentum ---
    T rho = cell.computeDensity();
    Array<T, Descriptor<T>::d> j;
    cell.computeMomentum(j); // j = sum_i f_i * e_i

    // --- Step 2: Retrieve force field from external field ---
    Array<T, Descriptor<T>::d> F_s =
        cell.template getExternal<Array<T, Descriptor<T>::d>>(FORCE_FIELD);

    // --- Step 3: Compute corrected velocity (Guo's formula) ---
    Array<T, Descriptor<T>::d> u;
    for (int d = 0; d < Descriptor<T>::d; ++d) {
      u[d] = j[d] / rho + 0.5 * F_s[d] / rho; // Δt = 1 in lattice units
    }

    // --- Step 4: Collision step with equilibrium and forcing ---
    T omega = this->getOmega();
    const T cs2 = Descriptor<T>::cs2;
    T jSqr = dot(j, j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
      const Array<T, Descriptor<T>::d> &e_i = Descriptor<T>::c[iPop];
      T w_i = Descriptor<T>::t[iPop];

      // Compute equilibrium
      T feq = this->computeEquilibrium(iPop, rho, j, jSqr);

      // Guo forcing term
      Array<T, Descriptor<T>::d> e_minus_u = e_i - u;
      T S_i = w_i * dot(e_minus_u, F_s) / cs2;

      // Standard BGK + forcing
      T f_old = cell[iPop];
      cell[iPop] = f_old - omega * (f_old - feq) + (1.0 - 0.5 * omega) * S_i;
    }
  }

  // ---- Equilibrium computation ----
  T computeEquilibrium(plint iPop, T rhoBar,
                       Array<T, Descriptor<T>::d> const &j,
                       T jSqr) const override {
    const T cs2 = Descriptor<T>::cs2;
    const T invCs2 = Descriptor<T>::invCs2;
    const T w_i = Descriptor<T>::t[iPop];
    const Array<T, Descriptor<T>::d> &e_i = Descriptor<T>::c[iPop];

    Array<T, 2> u(j[0] / rhoBar, j[1] / rhoBar);
    T eu = dot(e_i, u);
    T uu = dot(u, u);

    // Compute equilibrium distribution
    T Gamma_i =
        (1.0 + eu / cs2 + 0.5 * (eu * eu) / (cs2 * cs2) - 0.5 * uu / cs2);

    T feq = w_i * rhoBar * Gamma_i;
    return feq;
  }

  DynamicsMomentum<T, Descriptor> *clone() const override {
    return new DynamicsMomentum<T, Descriptor>(*this);
  }
};

#endif
