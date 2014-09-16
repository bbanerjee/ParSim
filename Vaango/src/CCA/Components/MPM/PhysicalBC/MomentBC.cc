/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Geometry/BBox.h>
#include <Core/Math/Matrix3.h>
#include <iostream>

using namespace Uintah;
using namespace std;

// Store the geometry object and the load curve
MomentBC::MomentBC(ProblemSpecP& ps, const GridP& grid, const MPMFlags* flags)
{
	// First read the geometry information
	// d_surface is the geometry object containing the surface to be loaded.
	// **WARNING** Currently allows only for box, cylinder or sphere.
	if (flags->d_useCBDI){
		ps->require("outward_normal",d_outwardNormal);
	}
	d_dxpp = Vector(1.,1.,1.);  // Only needed for axisymmetric end, see below
	ProblemSpecP geom_adult = ps->findBlock("geom_object");
	ProblemSpecP geom_child = geom_adult->findBlock();
	std::string go_type = geom_child->getNodeName();
	// Read normal plane information.
	ProblemSpecP norm_adult = ps->findBlock("normal_plane");
	norm_adult->require("normal", d_norm_norm);
	Vector d_norm_pt(0.0, 0.0, 0.0);
	norm_adult->require("point", d_norm_pt);
	d_norm_d = Dot(-d_norm_norm, d_norm_pt);

	// TODO sphere and cylinder
	if (go_type == "box") {
		d_surface = scinew BoxGeometryPiece(geom_child);
		//Box box = d_surface->getBoundingBox();
		d_surfaceType = "box";

		// Compute the length of the edge the moment is being applied to.
		BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
		Vector normal(0.0, 0.0, 0.0);
		normal[gp->thicknessDirection()] = 1.0; // Normal to surface where force is applied.

		// Find the direction of the normal and surface planes, and the length along where the force is applied.
		Vector upper(0.0, 0.0, 0.0), lower(0.0, 0.0, 0.0);
		// Normal in the x direction.
		if (d_norm_norm.x()!=0 && d_norm_norm.y()==0 && d_norm_norm.z()==0) {
			upper[0] = gp->getBoundingBox().upper().x();
			lower[0] = gp->getBoundingBox().lower().x();
		}
		// Normal in the y direction.
		else if (d_norm_norm.x()==0 && d_norm_norm.y()!=0 && d_norm_norm.z()==0) {
			upper[1] = gp->getBoundingBox().upper().y();
			lower[1] = gp->getBoundingBox().lower().y();
		}
		// Normal in the z direction.
		else if (d_norm_norm.x()==0 && d_norm_norm.y()==0 && d_norm_norm.z()!=0) {
			upper[2] = gp->getBoundingBox().upper().z();
			lower[2] = gp->getBoundingBox().lower().z();
		}
		else {
			cout << "error" << endl; // TODO error.
		}

		// Compute L1 and L2, distance from normal plane to either side of the boundary.
		double d_norm_L1 = (Dot(d_norm_norm, lower)+d_norm_d) / sqrt(Dot(d_norm_norm, d_norm_norm));
		double d_norm_L2 = (Dot(d_norm_norm, upper)+d_norm_d) / sqrt(Dot(d_norm_norm, d_norm_norm));
		d_norm_L1L2 = pow(abs(d_norm_L1),3) + pow(abs(d_norm_L2),3);

		// TODO DC comment
		/* cout << "plane normal: " << d_norm_norm.x() << " " << d_norm_norm.y() << " " << d_norm_norm.z() << endl;
		cout << "face normal: " << normal.x() << " " << normal.y() << " " << normal.z() << endl;
		cout << "box: [" << gp->getBoundingBox().lower().x() << ", " << gp->getBoundingBox().lower().y() << ", " << gp->getBoundingBox().lower().z() <<
				"], [" << gp->getBoundingBox().upper().x() << ", " << gp->getBoundingBox().upper().y() << ", " << gp->getBoundingBox().upper().z() << "]" << endl;
		cout << "lower: " << lower.x() << " " << lower.y() << " " << lower.z() << " " << endl;
		cout << "upper: " << upper.x() << " " << upper.y() << " " << upper.z() << " " << endl;
		cout << "d_norm_L1: " << d_norm_L1 << " d_norm_L2: " << d_norm_L2 << endl; */

	} else {
		throw ParameterNotFound("** ERROR ** No surface specified for moment BC.",
				__FILE__, __LINE__);
	}

	d_numMaterialPoints = 0;  // this value is read in on a restart
	ps->get("numberOfParticlesOnLoadSurface",d_numMaterialPoints);

	// Read and save the load curve information
	d_loadCurve = scinew LoadCurve<double>(ps);

	//__________________________________
	//   Bulletproofing
	// user shouldn't specify a geometry object that is bigger than the domain
	Box boundingBox = d_surface->getBoundingBox();
	BBox compDomain;
	grid->getSpatialRange(compDomain);

	Point BB_min = boundingBox.lower();
	Point CD_min = compDomain.min();
	Point BB_max = boundingBox.upper();
	Point CD_max = compDomain.max();

	if( ( BB_min.x() < CD_min.x() ) ||
			( BB_min.y() < CD_min.y() ) ||
			( BB_min.z() < CD_min.z() ) ||
			( BB_max.x() > CD_max.x() ) ||
			( BB_max.y() > CD_max.y() ) ||
			( BB_max.z() > CD_max.z() ) ){
		if(!d_axisymmetric_end && !d_axisymmetric_side){
			proc0cout <<"_________________________________________________________\n";
			proc0cout << "\n Input File WARNING: <PhysicalBC : MPM : Moment> \n"
					<< " The geometry Object ["<<d_surface->getType() << "] exceeds the dimensions of the computational domain.\n"
					<< " \n Please change the parameters so it doesn't. \n\n"
					<< " There is a flaw in the surface area calculation for the geometry object,\n"
					<< " it does not take into account that the object exceeds the domain\n";
			proc0cout <<"_________________________________________________________\n";
		}
	}
}

// Destroy the moment BCs
MomentBC::~MomentBC()
{
	delete d_surface;
	delete d_loadCurve;
}

void MomentBC::outputProblemSpec(ProblemSpecP& ps)
{
	ProblemSpecP momen_ps = ps->appendChild("moment");
	ProblemSpecP geom_ps = momen_ps->appendChild("geom_object");
	d_surface->outputProblemSpec(geom_ps);
	momen_ps->appendElement("numberOfParticlesOnLoadSurface",d_numMaterialPoints);
	d_loadCurve->outputProblemSpec(momen_ps);
	momen_ps->appendElement("res",d_res);
}

// Get the type of this object for BC application
std::string 
MomentBC::getType() const
{
	return "Moment";
}

// Locate and flag the material points to which this moment BC is
// to be applied. Assumes that the "checkForSurface" function in ParticleCreator.cc
// has been used to identify this material point as being on the surface of the body.
// WARNING : For this logic to work, the surface object should be a 
// box (zero volume), cylinder, sphere geometry piece that touches
// contains the surface on which the moment is to be applied.
bool
MomentBC::flagMaterialPoint(const Point& p,
		const Vector& dxpp)
{
	bool flag = false;
	// TODO sphere and cylinder
	if (d_surfaceType == "box") {
		// Create box that is min-dxpp, max+dxpp;
		Box box = d_surface->getBoundingBox();
		GeometryPiece* volume = scinew BoxGeometryPiece(box.lower()-dxpp,
				box.upper()+dxpp);

		if (volume->inside(p)) flag = true;
		delete volume;

	} else {
		throw ParameterNotFound("ERROR: Unknown surface specified for moment BC",
				__FILE__, __LINE__);
	}

	return flag;
}

// Calculate the force per particle at a certain time
double 
MomentBC::forcePerParticle(double time) const
{
	if (d_numMaterialPoints < 1) return 0.0;

	// Get the initial moment that is applied ( t = 0.0 )
	return moment(time);
}

// Calculate the force vector to be applied to a particular
// material point location
Vector
MomentBC::getForceVector(const Point& px, double forcePerParticle,
		const double time) const
{
	Vector force(0.0,0.0,0.0);
	// TODO sphere and cylinder
	if (d_surfaceType == "box") {
		BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
		Vector normal(0.0, 0.0, 0.0);
		normal[gp->thicknessDirection()] = 1.0; // Normal to surface where force is applied.

		/////////////////////////////////////////////////////////////////
		// Compute the force for the particle. The denominator should always equal 1 if the normal matrix is normalised (it should be...).
		double px_dist = (Dot(px, d_norm_norm)+d_norm_d) / sqrt(Dot(d_norm_norm, d_norm_norm));
		double force_m = (3*px_dist*forcePerParticle) / d_norm_L1L2;

		force = normal*force_m;

		// TODO DC comment
		/* cout << "px_dist: " << px_dist << endl;
		cout << "force_m: " << force_m << endl;
		cout << "force: " << force.x() << " " << force.y() << " " << force.z() << endl;
		cout << "L1^3+L2^3: " << d_norm_L1L2 << endl;
		cout << "forcePerParticle: " << forcePerParticle << endl; */

	} else {
		throw ParameterNotFound("ERROR: Unknown surface specified for moment BC",
				__FILE__, __LINE__);
	}
	return force;
}

// TODO DC
// Calculate the force vector to be applied to a particular
// material point location
Vector
MomentBC::getForceVectorCBDI(const Point& px, const Matrix3& psize,
		const Matrix3& pDeformationMeasure,
		double forcePerParticle,const double time,
		Point& pExternalForceCorner1,
		Point& pExternalForceCorner2,
		Point& pExternalForceCorner3,
		Point& pExternalForceCorner4,
		const Vector& dxCell) const
{
	Vector force(0.0,0.0,0.0);
	Vector normal(0.0, 0.0, 0.0);
	// TODO sphere and cylinder
	if (d_surfaceType == "box") {
		BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
		normal[gp->thicknessDirection()] = 1.0;
		force = normal*forcePerParticle;
	} else {
		throw ParameterNotFound("ERROR: Unknown surface specified for moment BC",
				__FILE__, __LINE__);
	}
	// 25% of total particle force goes to each corner
	force = force*0.25;
	// modify the sign of force if outward normal is not correctly defined
	if (!d_outwardNormal) {
		force = force*(-1.0);
	}
	// determine four boundary-corners of the particle
	int i1=0,i2=0;
	Matrix3 dsize=pDeformationMeasure*psize;
	Point px1;
	for (int i = 0; i < 3; ++i) {
		Vector dummy=Vector(dsize(0,i)*dxCell[0],dsize(1,i)*dxCell[1],
				dsize(2,i)*dxCell[2])/2.0;
		if (abs(Dot(normal,dummy)/(normal.length()*dummy.length())-1.0)<0.1) {
			px1=Point(px.x()+dummy[0],px.y()+dummy[1],px.z()+dummy[2]);
			i1=(i+1)%3;
			i2=(i+2)%3;
		} else if (abs(Dot(normal,dummy)/(normal.length()*dummy.length())+1.0)<0.1) {
			Point px1(px.x()-dummy[0],px.y()-dummy[1],px.z()-dummy[2]);
			i1=(i+1)%3;
			i2=(i+2)%3;
		}
	}
	// px1 is the position of the center of the boundary particle face that is on the physical boundary.
	pExternalForceCorner1=Point(px1.x()-dsize(0,i1)*dxCell[0]/2.0-dsize(0,i2)*dxCell[0]/2.0,
			px1.y()-dsize(1,i1)*dxCell[1]/2.0-dsize(1,i2)*dxCell[1]/2.0,
			px1.z()-dsize(2,i1)*dxCell[2]/2.0-dsize(2,i2)*dxCell[2]/2.0);
	pExternalForceCorner2=Point(px1.x()+dsize(0,i1)*dxCell[0]/2.0-dsize(0,i2)*dxCell[0]/2.0,
			px1.y()+dsize(1,i1)*dxCell[1]/2.0-dsize(1,i2)*dxCell[1]/2.0,
			px1.z()+dsize(2,i1)*dxCell[2]/2.0-dsize(2,i2)*dxCell[2]/2.0);
	pExternalForceCorner3=Point(px1.x()-dsize(0,i1)*dxCell[0]/2.0+dsize(0,i2)*dxCell[0]/2.0,
			px1.y()-dsize(1,i1)*dxCell[1]/2.0+dsize(1,i2)*dxCell[1]/2.0,
			px1.z()-dsize(2,i1)*dxCell[2]/2.0+dsize(2,i2)*dxCell[2]/2.0);
	pExternalForceCorner4=Point(px1.x()+dsize(0,i1)*dxCell[0]/2.0+dsize(0,i2)*dxCell[0]/2.0,
			px1.y()+dsize(1,i1)*dxCell[1]/2.0+dsize(1,i2)*dxCell[1]/2.0,
			px1.z()+dsize(2,i1)*dxCell[2]/2.0+dsize(2,i2)*dxCell[2]/2.0);
	// Recalculate the force based on area changes (current vs. initial)
	Vector iniVec1(psize(0,i1),psize(1,i1),psize(2,i1));
	Vector iniVec2(psize(0,i2),psize(1,i2),psize(2,i2));
	Vector curVec1(dsize(0,i1),dsize(1,i1),dsize(2,i1));
	Vector curVec2(dsize(0,i2),dsize(1,i2),dsize(2,i2));
	Vector iniA = Cross(iniVec1,iniVec2);
	Vector curA = Cross(curVec1,curVec2);
	double iniArea=iniA.length();
	double curArea=curA.length();
	force=force*(curArea/iniArea);
	return force;
}

namespace Uintah {
// A method to print out the moment bcs
ostream& operator<<(ostream& out, const MomentBC& bc)
{
	out << "Begin MPM Moment BC # = " << bc.loadCurveID() << endl;
	std::string surfType = bc.getSurfaceType();
	out << "    Surface of application = " << surfType << endl;
	// TODO sphere and cylinder
	if (surfType == "box") {
		Box box = (bc.getSurface())->getBoundingBox();
		out << "        " << box << endl;
	}
	out << "    Time vs. Load = " << endl;
	LoadCurve<double>* lc = bc.getLoadCurve();
	int numPts = lc->numberOfPointsOnLoadCurve();
	for (int ii = 0; ii < numPts; ++ii) {
		out << "        time = " << lc->getTime(ii)
        		 << " moment = " << lc->getLoad(ii) << endl;
	}
	out << "End MPM Moment BC # = " << bc.loadCurveID() << endl;
	return out;
}

} // end namespace Uintah
