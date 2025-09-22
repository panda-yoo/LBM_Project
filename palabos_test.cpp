
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


//#include "diffusion_d2q5.h"
//#include "vortex_shredding.h"
//#include "cylinder2d.h"
#include "palabos2D.h"
#include "palabos2D.hh"
#include <omp.h>
#include "obstacle.h"
// #include "droplet1.h" 
using namespace plb;


int main() {
    pcout << "MPI ranks: " << global::mpi().getSize() << std::endl;
    //pcout << "Threads  : " << plb::global::openMPlb().getSize() << std::endl;
    // sim();
    sim_obstacle();
    //sim_d2q5();
    //sim_obstacle();
    //sim_2d_cylinder();
    return 0;
}

