#ifndef COMPUTE_LAPLACIAN_H
#define COMPUTE_LAPLACIAN_H

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;

template <typename T>
class BoxLaplacianFunctional2D
    : public BoundedBoxProcessingFunctional2D_ST<T, T, 1> {
public:
  void processBulk(Box2D domain, ScalarField2D<T> &phi,
                   ScalarField2D<T> &laplacian) override {
    for (plint iX = domain.x0 + 1; iX <= domain.x1 - 1; ++iX) {
      for (plint iY = domain.y0 + 1; iY <= domain.y1 - 1; ++iY) {
        T center = phi.get(iX, iY);
        T lap = (phi.get(iX + 1, iY) + phi.get(iX - 1, iY) +
                 phi.get(iX, iY + 1) + phi.get(iX, iY - 1) - 4.0 * center);
        laplacian.get(iX, iY) = lap;
      }
    }
  }

  void processEdge(int, int, Box2D domain, ScalarField2D<T> &phi,
                   ScalarField2D<T> &laplacian) override {
    processBulk(domain, phi, laplacian);
  }

  void processCorner(int, int, Box2D domain, ScalarField2D<T> &phi,
                     ScalarField2D<T> &laplacian) override {
    processBulk(domain, phi, laplacian);
  }

  BoxLaplacianFunctional2D<T> *clone() const override {
    return new BoxLaplacianFunctional2D<T>(*this);
  }

  void
  getTypeOfModification(std::vector<modif::ModifT> &modified) const override {
    modified.push_back(modif::staticVariables);
  }

  BlockDomain::DomainT appliesTo() const override {
    return BlockDomain::bulkAndEnvelope;
  }
};

#endif
