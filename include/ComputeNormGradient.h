#ifndef COMPUTE_NORM_GRADIENT_H
#define COMPUTE_NORM_GRADIENT_H

#include "palabos2D.h"
#include "palabos2D.hh"
#include <algorithm>
#include <cmath>

using namespace plb;

template <typename T>
class BoxNormGradientFunctional2D
    : public BoundedBoxProcessingFunctional2D_ST<T, T, 2> {
public:
  // ------------------- BULK REGION -------------------
  virtual void processBulk(Box2D domain, ScalarField2D<T> &phi,
                           TensorField2D<T, 2> &normGrad) override {
    for (plint iX = domain.x0 + 1; iX <= domain.x1 - 1; ++iX) {
      for (plint jY = domain.y0 + 1; jY <= domain.y1 - 1; ++jY) {

        // Central difference approximation
        T dphidx = (phi.get(iX + 1, jY) - phi.get(iX - 1, jY)) / (T)2;
        T dphidy = (phi.get(iX, jY + 1) - phi.get(iX, jY - 1)) / (T)2;

        // Compute gradient magnitude
        T mag = std::sqrt(dphidx * dphidx + dphidy * dphidy + (T)1e-16);

        // Normalize gradient components
        normGrad.get(iX, jY)[0] = dphidx / mag;
        normGrad.get(iX, jY)[1] = dphidy / mag;
      }
    }
  }

  // ------------------- EDGE REGIONS -------------------
  virtual void processEdge(int direction, int orientation, Box2D domain,
                           ScalarField2D<T> &phi,
                           TensorField2D<T, 2> &normGrad) override {
    // Simple one-sided finite differences at boundaries
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint jY = domain.y0; jY <= domain.y1; ++jY) {
        T dphidx =
            (phi.get(std::min((plint)(iX + 1), (plint)(phi.getNx() - 1)), jY) -
             phi.get(std::max((plint)(iX - 1), (plint)0), jY)) /
            (T)2;

        T dphidy =
            (phi.get(iX, std::min((plint)(jY + 1), (plint)(phi.getNy() - 1))) -
             phi.get(iX, std::max((plint)(jY - 1), (plint)0))) /
            (T)2;

        // Compute normalized components
        T mag = std::sqrt(dphidx * dphidx + dphidy * dphidy + (T)1e-16);
        normGrad.get(iX, jY)[0] = dphidx / mag;
        normGrad.get(iX, jY)[1] = dphidy / mag;
      }
    }
  }

  // ------------------- CORNERS -------------------
  virtual void processCorner(int normalX, int normalY, Box2D domain,
                             ScalarField2D<T> &phi,
                             TensorField2D<T, 2> &normGrad) override {
    processEdge(0, 0, domain, phi, normGrad); // reuse one-sided gradient
  }

  // Clone function required by Palabos
  virtual BoxNormGradientFunctional2D<T> *clone() const override {
    return new BoxNormGradientFunctional2D<T>(*this);
  }

  // Palabos metadata: this functional modifies static (non-distribution)
  // variables
  virtual void
  getTypeOfModification(std::vector<modif::ModifT> &modified) const override {
    modified.push_back(modif::staticVariables);
  }

  virtual BlockDomain::DomainT appliesTo() const override {
    return BlockDomain::bulkAndEnvelope;
  }
};

#endif
