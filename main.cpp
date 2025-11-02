
// #include "atomicBlock/reductiveDataProcessorWrapper2D.h"

#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "io/parallelIO.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiDataField2D.h"

#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "palabos2D.h"
#include "palabos2D.hh"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

#include "ComputeNormGradient.h"

using namespace plb;
typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

template <typename U>
void writeField(plb_ofstream &file, TensorField2D<U, 2> &nhat) {

  plint nx = nhat.getNx();
  plint ny = nhat.getNy();
  file << 'x' << ' ' << 'y' << ' ' << "nx" << ' ' << "ny" << std::endl;
  for (plint iX = 0; iX < nx; ++iX) {
    for (plint iY = 0; iY < ny; ++iY) {

      file << iX << ' ' << iY << ' ' << nhat.get(iX, iY)[0] << ' '
           << nhat.get(iX, iY)[1] << std::endl;
    }
  }
}

template <typename U>
void writeDatFile(plb_ofstream &file, ScalarField2D<U> &field) {
  plint nx = field.getNx();
  plint ny = field.getNy();

  file << 'x' << ' ' << 'y' << ' ' << "field" << std::endl;

  for (plint i = 0; i < nx; i++) {
    for (plint j = 0; j < ny; j++) {
      file << i << ' ' << j << ' ' << field.get(i, j) << std::endl;
    }
  }
}

// ---------------- Initialization functional ----------------
template <typename T, template <typename U> class Descriptor>
class InitializePhiFunctional
    : public BoxProcessingFunctional2D_L<T, Descriptor> {
private:
  T r0, zeta;
  plint cx, cy;

public:
  InitializePhiFunctional(T r0_, T zeta_, plint cx_, plint cy_)
      : r0(r0_), zeta(zeta_), cx(cx_), cy(cy_) {}

  void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) override {
    Array<T, 2> u(0.0, 0.0);

    for (plint i = domain.x0; i <= domain.x1; ++i) {
      for (plint j = domain.y0; j <= domain.y1; ++j) {
        T dx = i - cx;
        T dy = j - cy;
        T mag = std::sqrt(dx * dx + dy * dy);
        T arg = (2.0 / zeta) * (mag - r0);
        if (arg > 50.0)
          arg = 50.0;
        T phi_val = 0.5 * (1.0 - std::tanh(arg));

        phi_val = std::max(T(1e-8), std::min(phi_val, T(1.0)));

        Cell<T, Descriptor> &cell = lattice.get(i, j);
        cell.defineDensity(phi_val);
        cell.defineVelocity(u);
        iniCellAtEquilibrium(cell, phi_val, u);
      }
    }
  }

  void
  getTypeOfModification(std::vector<modif::ModifT> &modified) const override {
    modified.push_back(modif::staticVariables); // <== Required by Palabos
  }

  InitializePhiFunctional<T, Descriptor> *clone() const override {
    return new InitializePhiFunctional<T, Descriptor>(*this);
  }
};

// ---------------- Main program ----------------
int main() {
  const T OMEGA = 1.0;
  const plint nx = 200, ny = 200;
  const T r0 = 40.0, zeta = 2.0;

  std::filesystem::create_directories("./data");
  global::directories().setOutputDir("./data");

  BlockLattice2D<T, DESCRIPTOR> phi(nx, ny,
                                    new BGKdynamics<T, DESCRIPTOR>(OMEGA));
  ScalarField2D<T> phi_s(nx, ny);
  TensorField2D<T, 2> nhat(nx, ny);

  std::cout << "Before initialization..." << std::endl;

  applyProcessingFunctional(
      new InitializePhiFunctional<T, DESCRIPTOR>(r0, zeta, nx / 2, ny / 2),
      phi.getBoundingBox(), phi);

  // Ensure lattice-level initialization/communication BEFORE computing
  // densities
  phi.initialize();

  // compute density into multi-scalar field
  computeDensity(phi, phi_s);
  std::vector<AtomicBlock2D *> args;
  args.push_back(&phi_s);
  args.push_back(&nhat);
  applyProcessingFunctional(new BoxNormGradientFunctional2D<T>(),
                            phi_s.getBoundingBox(), args, 0);

  std::cout << "After initialization." << std::endl;
  std::cout << "Center density = " << phi.get(nx / 2, ny / 2).computeDensity()
            << std::endl;

  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledPpm("phi_field", *computeDensity(phi));

  plb_ofstream file("../data/data_phi_field.dat");
  writeDatFile(file, *computeDensity(phi));

  plb_ofstream file2("../data/nhat_field.dat");

  writeField(file2, nhat);

  //   for (plint i = 0; i < nx; i++) {
  //     for (plint j = 0; j < ny; j++) {
  //       file << i << ' ' << j << ' ' << phi.get(i, j).computeDensity()
  //            << std::endl;
  //     }
  //   }
  file2.close();

  file.close();
  return 0;
}
