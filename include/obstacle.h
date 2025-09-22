#ifndef OBSTACLE_h
#define OBSTACLE_h

#include <palabos2D.hh>
#include <palabos2D.h>

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor


template <typename T>
class ConstantDensity {
public:
	ConstantDensity(T density_) : density(density_) {}
	T operator()(plb::plint, plb::plint) const { return density; }

private:
	T density;
};

template <typename T>
class Obstacle :public plb::DomainFunctional2D {
public:
	Obstacle(plb::plint x_start, plb::plint x_end, plb::plint y_start, plb::plint y_end);
	virtual bool operator()(plb::plint iX, plb::plint iY) const override;
	virtual Obstacle<T>* clone() const override;

private:
	plb::plint x_0, x_1, y_0, y_1;
};


void setup (plb::MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
	plb::OnLatticeBoundaryCondition2D<T, DESCRIPTOR>& boundarycondtions,
	plb::Array<T, 2> velocity,
	plb::Array<plb::plint, 4> obstacle_idx1,
	plb::Array<plb::plint, 4> obstacle_idx2);


void sim_obstacle();

#endif // !OBSTACLE_h
