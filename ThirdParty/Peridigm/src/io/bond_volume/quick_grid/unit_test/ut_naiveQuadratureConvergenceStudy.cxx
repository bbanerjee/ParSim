//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "mesh_output/Field.h"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "bond_volume/quick_grid/calculators.h"
#include "material_utilities.h"
#include "pdneigh/NeighborhoodList.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/BondFilter.h"
#include "utilities/Array.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace boost::unit_test;
using std::size_t;
using std::tr1::shared_ptr;
using UTILITIES::Array;
using UTILITIES::Vector3D;
using std::pair;
using std::cout;
using std::endl;

static int nx;
static int ny;
static int nz;
static double horizon;
const size_t numProcs=1;
const size_t myRank=0;

void probe_shear
(
		MATERIAL_EVALUATION::PURE_SHEAR mode,
		Array<int> neighborhoodPtr,
		Array<double> X,
		Array<double> xPtr,
		Array<double> Y,
		Array<double> yPtr,
		Array<double> bondVolume,
		double horizon,
		double gamma,
		double m_code
);

void scf_probe(const std::string& json_filename);

void set_static_data(const std::string& json_filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	try {
		read_json(json_filename, pt);
	} catch(std::exception& e){
		std::cerr << e.what();
		std::exit(1);
	}

	/*
	 * Get Discretization
	 */
	ptree discretization_tree=pt.find("Discretization")->second;
	std::string path=discretization_tree.get<std::string>("Type");
	horizon=pt.get<double>("Discretization.Horizon");

	if("QuickGrid.TensorProduct3DMeshGenerator"==path){

		nx = pt.get<int>(path+".Number Points X");
		ny = pt.get<int>(path+".Number Points Y");
		nz = pt.get<int>(path+".Number Points Z");

	} else {
		std::string s;
		s = "Error-->ut_naiveQuadratureConvergenceStudy\n";
		s += "\tTest only works for Discretization.Type==QuickGrid.TensorProduct3DMeshGenerator\n";
		throw std::runtime_error(s);
	}

}


void write_table_1_header(const std::string& output_tex_table){
	std::stringstream table_out;

	table_out << "\\begin{table}[ht]" << "\n";
	table_out << "\\centering" << "\n";
	table_out << "\\bigskip" << "\n";
	table_out << "\\begin{tabular}{|c|c|c|c|}" << "\n";
	table_out << "\\hline" << "\n";
	table_out << "$n$ "
			    << "& $\\frac{|m-m_n|}{m}$ "
			    << "& $\\frac{\\Vert e^d\\Vert^2-\\Vert e^d_n\\Vert^2}{\\Vert e^d\\Vert^2}$ "
			    << "& $\\frac{\\Vert e^d\\Vert^2}{\\Vert e^d_n\\Vert^2}$ \\\\" << "\n";
	table_out << "\\hline" << "\n";


	std::ofstream file_stream;
	file_stream.open(output_tex_table.c_str(),std::ios::app|std::ios::out);

	file_stream << table_out.str();
	file_stream.close();

}

void close_table_1(const std::string& output_tex_table) {
	std::stringstream table_out;

	table_out << "\\hline" << "\n";
	table_out << "\\end{tabular}" << "\n";
	table_out << "\\end{table}" << "\n";
	std::ofstream file_stream;
	file_stream.open(output_tex_table.c_str(),std::ios::app|std::ios::out);
	file_stream << table_out.str();
	file_stream.close();
}


QUICKGRID::QuickGridData getGrid(const string& _json_filename) {
	shared_ptr<QUICKGRID::QuickGridMeshGenerationIterator> g;
	g = QUICKGRID::getMeshGenerator(numProcs,_json_filename);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, *g);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

void ut_naiveQuadratureConvergenceStudy_n3() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=3.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n5() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=5.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n7() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=7.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n9() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=9.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n11() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=11.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n13() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=13.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n17() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=17.json";
	set_static_data(file);
	scf_probe(file);
}

void ut_naiveQuadratureConvergenceStudy_n33() {
	std::string file = "./input_files/ut_bondVolumeConvergenceStudy_n=33.json";
	set_static_data(file);
	scf_probe(file);
}



/*
 * Dave's Influence Function
 * "x < 0.5 ? 1.0 : -4.0*x*x + 4.0*x"
 */

void scf_probe(const std::string& json_filename) {

	QUICKGRID::QuickGridData gridData = getGrid(json_filename);

	// This load-balances
	gridData = PDNEIGH::getLoadBalancedDiscretization(gridData);

	/*
	 * Create neighborhood with an enlarged horizon
	 */

	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon);


	/*
	 * Unit test looks exclusively at the cell at the center of cube;
	 * This cell ID depends upon nx, ny, nz
	 *
	 * MESH INPUT MUST HAVE EVEN NUMBER OF CELLS
	 */
	BOOST_CHECK(0==(nx+1)%2);
	/*
	 * mesh must have equal number of cells along each axis
	 */
	BOOST_CHECK(nx==ny);
	BOOST_CHECK(nx==nz);

	// coordinate indices of center cell
	size_t ic = (nx -1)/2;
	size_t jc = ic;
	size_t kc = ic;
	size_t gId = nx * ny * kc + nx * jc + ic;
//	std::cout << "ut_scf::center cell gID = " << gId << std::endl;

	/**
	 * WARNING: note following ASSUMPTION -- gId == local id
	 * CAUTION: this test only makes sense in 'serial' -- local id
	 * and gId are not the same in parallel.
	 */
	size_t num_neigh = list.get_num_neigh(gId);
//	std::cout << "ut_scf::center cell num_neigh = " << num_neigh << std::endl;
	Array<int> neighborhoodPtr(1+num_neigh);
	{
		/*
		 * copy neighborhood list for center point over to array
		 */
		const int *neighborhood = list.get_neighborhood(gId);
		BOOST_CHECK((int)num_neigh == *neighborhood);
		for(size_t j=0;j<num_neigh+1;j++,neighborhood++)
			neighborhoodPtr[j]=*neighborhood;
	}
	Array<double> xPtr(list.get_num_owned_points()*3,list.get_owned_x());
//	shared_ptr<double> xPtr = list.get_owned_x();


	/*
	 * expectation is that cell center is at origin
	 */
//	std::cout << "ut_scf::cell center coordinate X = " << *(xPtr.get()+3*gId) << std::endl;
//	std::cout << "ut_scf::cell center coordinate Y = " << *(xPtr.get()+3*gId + 1) << std::endl;
//	std::cout << "ut_scf::cell center coordinate Z = " << *(xPtr.get()+3*gId + 2) << std::endl;
	BOOST_CHECK_SMALL(xPtr[3*gId+0],1.0e-15);
	BOOST_CHECK_SMALL(xPtr[3*gId+1],1.0e-15);
	BOOST_CHECK_SMALL(xPtr[3*gId+2],1.0e-15);
	/*
	 * X is the center of the sphere
	 */
	Array<double> X(3); X.set(0.0);
	/*
	 * Y = X since we are fixing the center of the sphere
	 */
	Array<double> Y(3); Y.set(0.0);

	Array<double> cellVolume(gridData.numPoints,gridData.cellVolume);
	double m_analytical = 4.0 * M_PI * pow(horizon,5) / 5.0;
	double m_code = MATERIAL_EVALUATION::computeWeightedVolume(X.get(),xPtr.get(),cellVolume.get(),neighborhoodPtr.get(),horizon);
	double rel_diff = std::abs(m_analytical-m_code)/m_analytical;
	std::cout << std::scientific;
	std::cout.precision(3);
	std::cout << "ut_scf::analytical value for weighted volume on sphere = " << m_analytical << std::endl;
	std::cout << "ut_scf::code computed weighted volume on sphere = " << m_code << std::endl;
	std::cout << "ut_scf::% relative error weighted volume = " << 100*rel_diff << std::endl;

	double gamma = 1.0e-6;
	Array<double> yPtr(3*list.get_num_owned_points());

	/*
	 * PROBE XY
	 */
	probe_shear(MATERIAL_EVALUATION::XY,neighborhoodPtr,X,xPtr,Y,yPtr,cellVolume,horizon,gamma,m_code);
	/*
	 * PROBE YZ
	 */
//	probe_shear(MATERIAL_EVALUATION::YZ,neighborhoodPtr,X,xPtr,Y,yPtr,volPtr,horizon,gamma,m_code);
	/*
	 * PROBE ZX
	 */
//	probe_shear(MATERIAL_EVALUATION::ZX,neighborhoodPtr,X,xPtr,Y,yPtr,volPtr,horizon,gamma,m_code);


}


void probe_shear
(
	MATERIAL_EVALUATION::PURE_SHEAR mode,
	Array<int> neighborhoodPtr,
	Array<double> X,
	Array<double> xPtr,
	Array<double> Y,
	Array<double> yPtr,
	Array<double> cellVolume,
	double horizon,
	double gamma,
	double m_code
)
{

	/*
	 * This is the reference value for the weighted volume
	 */
	double m_analytical = 4.0 * M_PI * pow(horizon,5) / 5.0;
	double m_err = std::fabs(m_analytical-m_code)/m_analytical;
	/*
	 * NOTE: X is center of sphere and there no displacement at this point
	 * therefore, Y=X
	 */
	MATERIAL_EVALUATION::set_pure_shear(neighborhoodPtr.get(),X.get(),xPtr.get(),yPtr.get(),mode,gamma);
	double theta = MATERIAL_EVALUATION::computeDilatation(neighborhoodPtr.get(),X.get(),xPtr.get(),X.get(),yPtr.get(),cellVolume.get(),m_code);
	std::cout << "ut_naiveQuadratureConvergenceStudy::probe_shear dilatation = " << theta << std::endl;
	double tolerance=1.0e-12;
	BOOST_CHECK_SMALL(theta,tolerance);

	/*
	 * compute shear correction factor
	 */
	/*
	 * This is the reference value for ed_squared
	 */
	double reference = 4.0 * M_PI * gamma * gamma * pow(horizon,5) / 75.0;
	double ed2 = MATERIAL_EVALUATION::compute_norm_2_deviatoric_extension(neighborhoodPtr.get(),X.get(),xPtr.get(),Y.get(),yPtr.get(),cellVolume.get(),m_code);
	double scf = reference/ed2;
	double ed_err = fabs(reference-ed2)/reference;
	std::cout << "ut_scf::probe_shear MODE = " << mode << std::endl;
	std::cout << "ut_scf::ed^2 = " << reference << std::endl;
	cout.precision(2);
	std::cout << std::scientific << "ut_scf::probe_shear computed % ed_err in pure shear = " << 100*ed_err << std::endl;
	std::cout << "ut_scf::probe_shear computed scf in pure shear = " << scf << std::endl;

	std::stringstream table_1_out;
	table_1_out << nx << " & ";
	table_1_out.precision(4);
	table_1_out << m_err*100 << "\\% & ";
	table_1_out << ed_err*100 << "\\% & ";
	table_1_out.precision(3);
	table_1_out << scf << " \\\\ \n";

	/*
	 * write latex table
	 */
	std::ofstream file_stream;
	file_stream.open("naive_table_1.tex",std::ios::app|std::ios::out);
	file_stream << table_1_out.str();
	file_stream.close();

	/*
	 * write raw data
	 */
	file_stream.open("ut_naiveQuadratureConvergenceStudy.dat",std::ios::app|std::ios::out);
	file_stream << nx << " ";
	file_stream << std::scientific;
	file_stream.precision(12);
	file_stream << 2.0*horizon/nx << " ";
	file_stream << m_code << " ";
	file_stream << ed2 << "\n";
	file_stream.close();

}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "ut_naiveQuadratureConvergenceStudy" );
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n3 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n5 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n7 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n9 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n11 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n13 ));
	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n17 ));
//	proc->add(BOOST_TEST_CASE( &ut_naiveQuadratureConvergenceStudy_n33 ));
	framework::master_test_suite().add( proc );
	return success;

}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(
		int argc,
		char* argv[]
)
{


	write_table_1_header("naive_table_1.tex");
//	write_table_2_header("table_2.tex");

	// Initialize UTF
	int flag = unit_test_main( init_unit_test, argc, argv );

	close_table_1("naive_table_1.tex");
	return flag;
}

