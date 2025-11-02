#### 1. If you want to compute some scalar quantity that depends on lattice variables (e.g., local density ρ) and write it into a scalar field φ, you would do:

``` 
template <typename T, template <typename U> class Descriptor>
class ComputeDensityFieldFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T, Descriptor>& lattice,
                         ScalarField2D<T>& field) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                field.get(iX, iY) = lattice.get(iX, iY).computeDensity();
    }

    virtual ComputeDensityFieldFunctional2D<T, Descriptor>* clone() const {
        return new ComputeDensityFieldFunctional2D<T, Descriptor>(*this);
    }
};

```
Then call:
```applyProcessingFunctional<T, Descriptor, T>(
    new ComputeDensityFieldFunctional2D<T, Descriptor>(),
    lattice.getBoundingBox(),
    lattice, phiField);
```