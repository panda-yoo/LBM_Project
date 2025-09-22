#ifndef DROPLET_1_H
#define DROPLET_1_H

#include <palabos2D.h>
#include <palabos2D.hh>

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

template<typename T>
class droplet1 : public plb::DomainFunctional2D
{

public :
	plb::plint cx, cy, radiusSqr;

public:
	droplet1(plb::plint cx_, plb::plint cy_, plb::plint radius);
	virtual bool operator()(plb::plint ix, plb::plint jy)const;
	virtual droplet1* clone() const;


};
void sim();
/*

class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
	CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
		: cx(cx_), cy(cy_), radiusSqr(plb::util::sqr(radius)) {
	}

	virtual bool operator()(plb::plint iX, plb::plint iY) const {
		return plb::util::sqr(iX - cx) + plb::util::sqr(iY - cy) <= radiusSqr;
	}

	virtual CylinderShapeDomain2D* clone() const {
		return new CylinderShapeDomain2D(*this);
	}

*/
#endif

