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

#include <palabos2D.h>
#include <palabos2D.hh>
#include "obstacle.h"
#include <filesystem>

template <typename T>
Obstacle<T>::Obstacle(plb::plint x_start, plb::plint x_end, plb::plint y_start, plb::plint y_end) 
	:x_0(x_start), x_1(x_end), y_0(y_start), y_1(y_end) {}


template<typename T>
bool Obstacle<T>::operator()(plb::plint iX, plb::plint iY) const {
	return iX >= x_0 && iX <= x_1 && iY >= y_0 && iY <= y_1;
}

template<typename T>
Obstacle<T>* Obstacle<T>::clone() const{
	return new Obstacle<T>(*this);

}

void setup(plb::MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
	plb::OnLatticeBoundaryCondition2D<T, DESCRIPTOR>& boundarycondtions,
	plb::Array<T,2> velocity,
	plb::Array<plb::plint,4> obstacle_idx1,
	plb::Array<plb::plint, 4> obstacle_idx2) {

	plb::plint nx = lattice.getNx();
	plb::plint ny = lattice.getNx();
	
	plb::Box2D inlet(0, 0, 1, ny - 2);
	plb::Box2D outlet(nx - 1, nx - 1, 1, ny - 2);

	plb::Box2D up_wall(0, nx - 1,0, 0);
	plb::Box2D down_wall(0, nx - 1, ny - 1, ny - 1);

	plb::Array<T, 2> zero_velocity(0., 0.);

	plb::setBoundaryVelocity(lattice, up_wall, zero_velocity);
	boundarycondtions.setVelocityConditionOnBlockBoundaries(lattice, up_wall);

	plb::setBoundaryVelocity(lattice, down_wall, zero_velocity);
	boundarycondtions.setVelocityConditionOnBlockBoundaries(lattice, down_wall);

	boundarycondtions.setVelocityConditionOnBlockBoundaries(lattice, inlet);
	plb::setBoundaryVelocity(lattice, inlet, velocity);

	
	boundarycondtions.setPressureConditionOnBlockBoundaries(lattice, outlet);
	plb::setBoundaryDensity(lattice, outlet, 1.0);


	plb::defineDynamics(lattice, lattice.getBoundingBox(),
		new Obstacle<T>(obstacle_idx1[0],obstacle_idx1[1], obstacle_idx1[2], obstacle_idx1[3]),
		new plb::BounceBack<T,DESCRIPTOR>);


	plb::defineDynamics(lattice, lattice.getBoundingBox(),
		new Obstacle<T>(obstacle_idx2[0], obstacle_idx2[1], obstacle_idx2[2], obstacle_idx2[3]),
		new plb::BounceBack<T, DESCRIPTOR>);

	plb::initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1.0, velocity);
	lattice.initialize();
}


void sim_obstacle() {
	
	std::filesystem::create_directory("sim_obstacle");

	std::filesystem::create_directory("sim_obstacle/2_obstacles");

	std::filesystem::create_directory("sim_obstacle/2_obstacles/velocity");
	std::filesystem::create_directory("sim_obstacle/2_obstacles/vorticity");


	plb::global::directories().setOutputDir("./sim_obstacle/");
	plb::plb_ofstream file("./sim_obstacle/2_obstacles/data.dat");
	plb::plb_ofstream file_ene("./sim_obstacle/2_obstacles/data_mass.dat");


	plb::plint lx = 200;
	plb::plint ly = 60;
	plb::plint ntimes= 1000;

	T omega = 1.0;
	plb::Array<T, 2> u(0.002, 0.0);

	//obstacle 1
	plb::Array<T, 4> obstacle_idx1(40,80,0,30);
	//obstacle 2
	plb::Array<T, 4> obstacle_idx2(100, 140, 30, 59);


	//lattice initialization 

	plb::MultiBlockLattice2D<T, DESCRIPTOR> lattice(lx,ly,
		new plb::BGKdynamics<T,DESCRIPTOR>(omega));

	plb::OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundary_condition = 
		plb::createLocalBoundaryCondition2D<T,DESCRIPTOR>();

	setup(lattice,*boundary_condition,u,obstacle_idx1, obstacle_idx2);

	plb::ImageWriter<T> img_("leeloo");
	//plb::ImageWriter<T> img_2("./obstacle/2_obstacles/vorticity/");


	//plb::TensorField2D<T, 2> *velocity_field;
	file_ene << "iteration" << ' ' << "Mass" << ' ' << std::endl;

	for (plb::plint tm = 0; tm < ntimes; ++tm)
	{
		if (tm%10==0)
		{
			
			plb::pcout << "step " << tm
					<< "; av energy =" << plb::getStoredAverageEnergy<T>(lattice)
					<< "; av rho =" << plb::getStoredAverageDensity<T>(lattice) << std::endl;
			
			file_ene << tm << ' ' << plb::getStoredAverageDensity<T>(lattice) << std::endl;
			img_.writeScaledGif(plb::createFileName("./2_obstacles/velocity/Velocity_", tm, 6), *plb::computeVelocityNorm(lattice));
			img_.writeScaledGif(plb::createFileName("./2_obstacles/vorticity/Vorticity_", tm, 6), *plb::computeVorticity(*plb::computeVelocity(lattice)));

			
		}
		lattice.collideAndStream();
	}


	plb::Array<T, 2> velocity(0.0,0.0);
	file <<  "x" << " " << "y" << " " <<"ux" << " " << "uy" << std::endl;

	for (plb::plint ix = 0; ix < lx; ++ix) {
	
		for (plb::plint jy = 0; jy < ly; ++jy) {
			plb::Cell cell = lattice.get(ix, jy);
			cell.computeVelocity(velocity);
			
			file << ix << " " << jy << " " << velocity[0] << " " << velocity[1] << std::endl;

		}
	}
	file.close();

}