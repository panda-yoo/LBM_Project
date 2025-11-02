
#include "palabos2D.hh"
#include "palabos2D.hh"

using namespace plb;
typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

T mag(Array<T, 2> vec)
{
    T mag_ = vec[0] * vec[0] + vec[1] * vec[1];
    return sqrt(mag_);
}
template <typename T, template <typename U> class Descriptor>
class InitializePhiFunctional : public BoxProcessingFunctional2D_L<T, Descriptor>
{
private:
    plint r, x0, y0;
    T zeta;

public:
    InitializePhiFunctional(plint r_, T zeta_, plint x0_, plint y0_)
        : r(r_), zeta(std::max(T(1e-3), zeta_)), x0(x0_), y0(y0_) {}

    void process(Box2D domain, BlockLattice2D<T, Descriptor> &phi) override
    {
        Array<T, 2> u(0.0, 0.0);

        for (plint i = domain.x0; i <= domain.x1; ++i)
            for (plint j = domain.y0; j <= domain.y1; ++j)
            {
                // compute phi_val
                T dx = i - x0;
                T dy = j - y0;
                T mag = std::sqrt(dx * dx + dy * dy);

                T arg = (2.0 / zeta) * (mag - r);
                const T maxArg = 50.0;
                if (arg > maxArg)
                    arg = maxArg;

                T phi_val = 0.5 * (1.0 - std::tanh(arg));
                phi_val = std::max(T(0.0), std::min(phi_val, T(1.0)));

                // Optionally guard tiny densities to avoid divide-by-zero later
                if (phi_val < T(1e-15))
                    phi_val = T(0.0);

                // get the cell and initialize to equilibrium
                Cell<T, Descriptor> &cell = phi.get(i, j);

                // Preferred style inside functors:
                cell.defineDensity(phi_val);
                cell.defineVelocity(u);
                iniCellAtEquilibrium(cell, phi_val, u);

                // If you need to set external fields or other scalars:
                // cell.setExternalField(pos, size, extArray);
            }
    }

    InitializePhiFunctional<T, Descriptor> *clone() const override
    {
        return new InitializePhiFunctional<T, Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT> &modified) const override
    {
        modified.push_back(modif::staticVariables);
    }
};
template <typename T, template <typename U> class Descriptor>
class InitializeDensityFunctional : public BoxProcessingFunctional2D_L<T, Descriptor>
{
private:
    T rho0, rhoGradient;

public:
    InitializeDensityFunctional(T rho0_, T rhoGradient_)
        : rho0(rho0_), rhoGradient(rhoGradient_) {}

    void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) override
    {
        Array<T, 2> u(0.0, 0.0); // start with zero velocity
        for (plint i = domain.x0; i <= domain.x1; ++i)
            for (plint j = domain.y0; j <= domain.y1; ++j)
            {
                T rho = rho0 + rhoGradient * j; // e.g., linear profile in y
                Cell<T, Descriptor> &cell = lattice.get(i, j);
                cell.defineDensity(rho);
                cell.defineVelocity(u);
                iniCellAtEquilibrium(cell, rho, u);
            }
    }

    void getTypeOfModification(std::vector<modif::ModifT> &modified) const override
    {
        modified.push_back(modif::staticVariables);
    }

    InitializeDensityFunctional<T, Descriptor> *clone() const override
    {
        return new InitializeDensityFunctional<T, Descriptor>(*this);
    }
};

// NEED A FUNCTOR IMPLEMENTATION
/*

void initialize_phi(MultiBlockLattice2D<T, DESCRIPTOR> &phi, plint r, T zeta,
plint x0, plint y0)
{

const plint nx = phi.getNx();
const plint ny = phi.getNy();
Array<T, 2> u(0.0, 0.0);

// Safety: zeta should not be too small
if (zeta < 1e-3)
zeta = 1e-3;

// Optional: open file to save phi
plb_ofstream file("../data/data_phi.dat");
file << "i j phi" << std::endl;

for (plint i = 0; i < nx; ++i)
{
    for (plint j = 0; j < ny; ++j)
    {
        T dx = i - x0;
        T dy = j - y0;
        T mag = std::sqrt(dx * dx + dy * dy);

        T arg = (2.0 / zeta) * (mag - r);
        const T maxArg = 50.0;
        if (arg > maxArg)
        arg = maxArg;

        T phi_val = 0.5 * (1.0 - std::tanh(arg));
        phi_val = std::max(T(0.0), std::min(phi_val, T(1.0)));

        Box2D cellBox(i, i, j, j);
        initializeAtEquilibrium(phi, cellBox, phi_val, Array<T, 2>(0.0, 0.0));

        file << i << ' ' << j << ' ' << phi_val << std::endl;
    }
}
}

*/