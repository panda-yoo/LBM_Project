#ifndef PHI_H
#define PHI_H

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

#define PHI_NORMGRAD_FIELD 20

template <typename T, template <typename U> class Descriptor>
class phi : public BGKdynamics<T, Descriptor> {

public:
  phi(T M, T zeta) : M_(M), zeta_(zeta) {}

  // must override this method ,if inhereted from dynamics class
  phi<T, Descriptor> *clone() const override {
    return new phi<T, Descriptor>(*this);
  }

  T computeEquilibrium(plint iPop, T rhoBar,
                       Array<T, Descriptor<T>::d> const &j,
                       T jSqr) const override {

    const T rho_eps = (std::abs(rhoBar) < (T)1e-18) ? (T)1e-18 : rhoBar;

    // calculate the usual equilibrium feq
    const T cs2 = Descriptor<T>::cs2;
    const T invCs2 = Descriptor<T>::invCs2;

    const T w_i = Descriptor<T>::t[iPop];
    const Array<T, Descriptor<T>::d> &e_i = Descriptor<T>::c[iPop];
    Array<T, 2> u(j[0] / rho_eps, j[1] / rho_eps);

    T eu = dot(e_i, u);
    T uu = dot(u, u);
    T Gamma_i =
        w_i * (1.0 + eu / cs2 + 0.5 * (eu * eu) / (cs2 * cs2) - 0.5 * uu / cs2);

    // phi i.e  rhobar -> 0th momentum
    T term1 = rhoBar * Gamma_i;
    Array<T, 2> nHat =
        this->getCell().template getExternal<Array<T, 2>>(PHI_NORMGRAD_FIELD);

    T enHat = dot(e_i, nHat);

    T term2 =
        w_i * M_ * invCs2 * (4.0 / zeta_) * rhoBar * (1.0 - rhoBar) * enHat;

    return term1 + term2;
  }

private:
  T M_, zeta_;
};

#endif