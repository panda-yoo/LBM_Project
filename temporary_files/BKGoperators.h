#ifndef BKG_OPERATORS
#define BKG_OPERATORS
#include <palabos2D.h>
#include <palabos2D.hh>
using namespace plb;



    /// A dynamics which reads the relaxation parameter from external scalar before collision.
template <typename T, template <typename U> class Descriptor>
class CA_Dynamics : public BasicBulkDynamics<T, Descriptor> {
public:
    CA_Dynamics(T omega_);
    virtual T getOmega() const;


    virtual CA_Dynamics<T, Descriptor>* clone() const ;

    virtual void collide(Cell<T,Descriptor> &cell, BlockStatistics &statistics);

    computeEquil

    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const& j,
        T jSqr,T thetaBar = T()) const override;
    virtual void computeVelocity(Cell<T,Descriptor> &cell,Array<T,Descriptor<T>::d> &u) const override;
       /* void computeEquilibria(
            Array<T, Descriptor<T>::q>& fEq, T rhoBar,
            Array<T, Descriptor<T>::d> const& j, T jSqr,
            T thetaBar = T()) const override
        {
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                fEq[iPop] = computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
        }*/
    private:
        T omega_val_;
    };



// In a new file: gradientProcessor.h

template<typename T>
class GradientProcessor : public BoxProcessingFunctional2D_S<T> {
public:
    GradientProcessor(MultiTensorField2D<T, 2>& gradientField_) :
        gradientField(gradientField_)
    {
    }

    // This is the core function that Palabos will execute on each cell
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField)
    {
        // Use a 2nd-order central difference scheme for the gradient
        for (plint iX = domain.x0 + 1; iX <= domain.x1 - 1; ++iX) {
            for (plint iY = domain.y0 + 1; iY <= domain.y1 - 1; ++iY) {
                

                    // Get scalar values at neighboring nodes
                    T phi_plus_x = scalarField.get(iX + 1, iY);
                    T phi_minus_x = scalarField.get(iX - 1, iY);
                    T phi_plus_y = scalarField.get(iX, iY + 1);
                    T phi_minus_y = scalarField.get(iX, iY - 1);
                  

                    // Compute partial derivatives. Since delta_x = 1 in lattice units,
                    // the formula is simplified.
                    T gradX = (phi_plus_x - phi_minus_x) / 2.0;
                    T gradY = (phi_plus_y - phi_minus_y) / 2.0;
                    
                    Array<T, 2>& grad = gradientField.get(iX, iY);
                    // Store the resulting gradient vector in the output field
                    T* grad = gradientField.get(iX, iY).getData();
                    grad = gradX;
                    grad[1] = gradY;
                    
                
            }
        }
    }

    // This method is required by the Palabos functional interface
    virtual GradientProcessor<T>* clone() const {
        return new GradientProcessor<T>(*this);
    }

private:
    MultiTensorField2D<T, 2>& gradientField;
};

template<typename T>
class LaplacianProcessor : public BoxProcessingFunctional2D_S<T> {
public:
    LaplacianProcessor(MultiScalarField2D<T>& laplacianField_) :
        laplacianField(laplacianField_)
    {
    }

    // The core function that Palabos executes on each cell
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField)
    {
        // Iterate over the domain, avoiding the boundaries where the stencil is not defined
        for (plint iX = domain.x0 + 1; iX <= domain.x1 - 1; ++iX) {
            for (plint iY = domain.y0 + 1; iY <= domain.y1 - 1; ++iY) {

                // Get scalar values at the central node and its four neighbors
                T phi_center = scalarField.get(iX, iY);
                T phi_plus_x = scalarField.get(iX + 1, iY);
                T phi_minus_x = scalarField.get(iX - 1, iY);
                T phi_plus_y = scalarField.get(iX, iY + 1);
                T phi_minus_y = scalarField.get(iX, iY - 1);

                // Compute the Laplacian using the 5-point stencil
                T laplacianValue = phi_plus_x + phi_minus_x + phi_plus_y + phi_minus_y - 4.0 * phi_center;

                // Store the result in the output scalar field
                laplacianField.get(iX, iY) = laplacianValue;
            }
        }
    }

    // Required clone method for the Palabos functional interface
    virtual LaplacianProcessor<T>* clone() const {
        return new LaplacianProcessor<T>(*this);
    }

private:
    MultiScalarField2D<T>& laplacianField;
};

//#endif // LAPLACIAN_PROCESSOR_H
//#endif // GRADIENT_PROCESSOR_H





#endif // !BKG_OPERATORS
