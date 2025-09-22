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


#include "droplet1.h"
#include <palabos2D.h>
#include <palabos2D.hh>
#include <filesystem>

template<typename T>
droplet1<T>::droplet1(plb::plint cx_, plb::plint cy_, plb::plint radius) :cx(cx_), cy(cy_), radiusSqr(plb::util::sqr(radius)){};

template<typename T>
bool droplet1<T>::operator()(plb::plint ix, plb::plint jy)const {

	return plb::util::sqr(ix - cx) + plb::util::sqr(jy - cy) >= radiusSqr;
}
template<typename T>
droplet1<T>* droplet1<T>::clone() const {
	return new droplet1<T>(*this);
}


void sim() {

	
	plb::plint nx = 400;
	plb::plint ny = 400;
	T omega = 1.0;
	plb::plint ntimes = 800;
	plb::plint image_per = 4;

	std::filesystem::create_directory("../droplet");
	std::filesystem::create_directory("../droplet/density");
	std::filesystem::create_directory("../droplet/velocity_norm");


	plb::global::directories().setOutputDir("../droplet/");

	plb::MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx,ny,new plb::BGKdynamics<T,DESCRIPTOR>(omega));
	
	plb::Box2D domain(nx/2, nx / 2, nx / 2, nx / 2);

	plb::defineDynamics(lattice,lattice.getBoundingBox(),new droplet1<T>(nx / 2, nx / 2, 180), new plb::BounceBack<T, DESCRIPTOR>);

	plb::initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1.0, plb::Array<T,2>(0.0,0.0));
	plb::initializeAtEquilibrium(lattice, domain, (T)40.0, plb::Array<T, 2>(0.0, 0.0));
	lattice.initialize();
	plb::ImageWriter<T> imgs("images");
	plb::plb_ofstream data("../droplet/data.dat");

	data << "iterations" << " " << "mass" << std::endl;

	for (plb::plint t = 0; t < ntimes; t++)
	{
		lattice.collideAndStream();

		if (t % image_per == 0) {
			imgs.writeScaledGif(plb::createFileName("./density/density", t, 4), *plb::computeDensity(lattice));
			imgs.writeScaledGif(plb::createFileName("./velocity_norm/velocity_norm", t, 4), *plb::computeVelocityNorm(lattice));

			plb::pcout << "iter -> " << t << " mass -> " << plb::computeAverageDensity(lattice)<<std::endl;
			data << t << " " << plb::computeAverageDensity(lattice) << std::endl;
		}
	}


}