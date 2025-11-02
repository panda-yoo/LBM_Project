

#ifndef CUSTOM_DYNAMICS_H
#define CUSTOM_DYNAMICS_H
#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

// #define DESCRIPTOR plb::descriptors::D2Q9Descriptor

template <typename T, template <typename U> class Descriptor>
class custom_dynamics : public plb::BGKdynamics<T, Descriptor>
{

public:
    custom_dynamics(T omega, T chi, T mu) : BGKdynamics<T, Descriptor>(omega), chi_(chi), mu_(mu) {}

    custom_dynamics<T, Descriptor> *clone() const
    {
        return new custom_dynamics<T, Descriptor>(*this);
    }

    T computeEquilibrium(plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr) const override
    {
        const T rho_eps = (std::abs(rhoBar) < (T)1e-18) ? (T)1e-18 : rhoBar;

        // calculate the usual equilibrium feq
        const T cs2 = Descriptor<T>::cs2;
        const T invCs2 = Descriptor<T>::invCs2;

        const T w_i = Descriptor<T>::t[iPop];
        const Array<T, Descriptor<T>::d> &e_i = Descriptor<T>::c[iPop];

        Array<T, 2> u(j[0] / rho_eps, j[1] / rho_eps);

        T eu = dot(e_i, u);
        T uu = dot(u, u);
        T Gamma_i = w_i * (1.0 + eu / cs2 +
                           0.5 * (eu * eu) / (cs2 * cs2) -
                           0.5 * uu / cs2);

        // T feq_i = w_i * ((chi_ * mu_ * invCs2) + rhoBar * (Gamma_i - w_i));
        if (iPop != 0)
        {
            T feq_i = w_i * ((chi_ * mu_ * invCs2) + rhoBar * (Gamma_i - w_i));
            return feq_i;
        }
        else
        {
            T sumFeq = (T)0;
            for (plint i = 1; i < Descriptor<T>::q; i++)
            {
                const T w = Descriptor<T>::t[i];
                const Array<T, Descriptor<T>::d> &e = Descriptor<T>::c[i];

                T eu_i = dot(e, u);
                T Gamma_i_local = w * ((T)1.0 +
                                       eu_i / cs2 +
                                       (T)0.5 * (eu_i * eu_i) / (cs2 * cs2) -
                                       (T)0.5 * uu / cs2);

                sumFeq += w * ((chi_ * mu_) * invCs2 + rhoBar * (Gamma_i_local - w));
            }
            T feq_i = rhoBar - sumFeq;

            return feq_i;
        }
    }

private:
    T mu_, chi_;
};

#endif