// Copyright [2025] [Pranav Shinde]
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

#include "palabos2D.h"
#include "palabos2D.hh"
#include "diffusion_d2q5.h"



void sim_d2q5() {

    plint nx = 100;
    plint ny = 100;
    plint itr = 1000;

    //MultiBlockLattice2D(plint nx, plint ny, Dynamics<T, Descriptor> *backgroundDynamics_);
    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>(1.0));
    lattice.periodicity().toggle(0, true);
    lattice.periodicity().toggle(1, true);


    initializeAtEquilibrium(lattice, Box2D(0, nx, 0, ny), 1.0, Array<T, 2>(0.0, 0.0));
    initializeAtEquilibrium(lattice, Box2D(nx / 2, nx / 2, ny / 2, ny / 2), 2.0, Array<T, 2>(0.0, 0.0));


    for (plint t = 0; t < itr; ++t) {

        lattice.collideAndStream();

    }

    plb_ofstream out("data.dat");

    for (plint ix = 0; ix < nx; ++ix) {
        for (plint jy = 0; jy < ny; ++jy) {

            Cell<T, DESCRIPTOR> cell = lattice.get(ix, jy);

            out << ix << " " << jy << " " << cell.computeDensity() << std::endl;


        }
    }
    out.close();
}
