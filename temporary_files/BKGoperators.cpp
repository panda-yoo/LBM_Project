#include <palabos2D.h>
#include <palabos2D.hh>
#include "BKGoperators.h"



template <typename T, template <typename U> class Descriptor>
CA_Dynamics<T,Descriptor>::CA_Dynamics(T omega_):omega_val_(omega_);

template <typename T, template <typename U> class Descriptor>
T CA_Dynamics<T, Descriptor>::getOmega() const {
        return omega_val_};


//standard polymorphic “virtual copy constructor” so Palabos can duplicate your dynamics object across the lattice
template <typename T, template <typename U> class Descriptor>
CA_Dynamics<T, Descriptor>* CA_Dynamics<T, Descriptor>::clone() const {
    return new CA_Dynamics<T, Descriptor>(*this);

};

template <typename T, template <typename U> class Descriptor>
void CA_Dynamics<T, Descriptor>::computeVelocity(Cell<T, Descriptor>& cell, Array<T, Descriptor<T>::d>& u) const {
    u.resetToZero();

    T cs2 = Descriptor<T>::cs2
    
        for (plint ipop = 0; ipop < Descriptor<T>::q; ++ipop)
    {
        T fi = cell[ipop];
        u += fi * Descriptor<T>::c[ipop];


    }
    T rho = computeDensity(cell);

    u = u / rho;


}



template <typename T, template <typename U> class Descriptor>
T CA_Dynamics<T, Descriptor>::computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const& j,
        T jSqr,T thetaBar = T()) const{
        // your formula for one direction
        const Array<T, Descriptor<T>::d>& c = Descriptor<T>::c[iPop];
        T eu = dotProduct(c, j) / rhoBar;
        return Descriptor<T>::t[iPop] * rhoBar * (1 + eu + 0.5 * eu * eu - 0.5 * jSqr / (rhoBar * rhoBar));
    }

template <typename T, template <typename U> class Descriptor>
void CA_Dynamics<T, Descriptor>::collide(Cell<T, Descriptor>& cell, plb::BlockStatistics& statistics) {
    T omega = this->getOmega();
    T rhoBar = this->computeRhoBar(cell);
    Array<T, Descriptor<T>::d> j;

    this->computeJ(cell, j);



}
   


