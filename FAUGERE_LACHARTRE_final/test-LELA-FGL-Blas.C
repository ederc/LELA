/*
 * test-LELA-FGL-Blas.C
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 *  Created on: 24 August 2012
 *      Author: martani (UPMC University Paris 06 / INRIA)
 */

#include "matrix-utils.h"
#include "lela/matrix/sparse.h"
#include "lela/matrix/dense.h"
#include "lela/util/commentator.h"
#include "../util/support.h"

#include "lela/algorithms/elimination.h"
#include "lela/algorithms/gauss-jordan.h"

#include "lela/algorithms/faugere-lachartre.h"

using namespace LELA;
using namespace std;

template <typename Ring>
void test_LELA_FGL_double_elts(const Ring& R, SparseMatrix<typename Ring::Element>& A, bool reduced)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	typedef Modular<double> Ring_Double;
	Ring_Double	R_double (R._modulus);

	//transform a 16bits matrix to a double represented matrix
	SparseMatrix<typename Ring_Double::Element> A_double (A.rowdim (), A.coldim ());
	MatrixUtils::invertMatrixRows(A);

	MatrixUtils::show_mem_usage("");

	commentator.start("uint16 -> double");
	for(uint32 i=0; i<A.rowdim (); ++i)
	{
		for(uint32 j=0; j<A[i].size (); ++j)
		{
			//A_double[i].push_back(typename Vector<Ring_Double>::Sparse::value_type (A[i][j].first, Ring_Double::Element ()));
			//R_double.init (A_double[i].back ().second, A[i][j].second);
			A_double[i].push_back(typename Vector<Ring_Double>::Sparse::value_type (A[i][j].first, A[i][j].second));
		}
	}
	commentator.stop(MSG_DONE);

	MatrixUtils::show_mem_usage("");

	for(uint32 i=0; i<A.rowdim (); ++i)
	{
		for(uint32 j=0; j<A[i].size (); ++j)
		{
			if(A[i][j].first != A_double[i][j].first || A[i][j].second != A_double[i][j].second)
			{
				report << "MATRIX NOT EQUAL at index" << A[i][j].first << endl;
				report << "elt A " << A[i][j].second << " | elt A_double " << A_double[i][j].second << endl;
				return;
			}
		}
	}

	MatrixUtils::show_mem_usage("");
	Context<Ring_Double> ctx (R_double);
	FaugereLachartre<Ring_Double> Solver (ctx);

	size_t rank;
	typename Ring_Double::Element det;

	commentator.start("Solver.echelonize");
		Solver.echelonize (A_double, A_double, rank, det, reduced);
	commentator.stop("Solver.echelonize");

	MatrixUtils::show_mem_usage("");

	report << "True Rank : " << rank << endl;
}

template <typename Ring>
void test_LELA_FGL_uint16_elts(const Ring& R, SparseMatrix<typename Ring::Element>& A, bool reduced)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	MatrixUtils::invertMatrixRows(A);
	MatrixUtils::show_mem_usage("");

	Context<Ring> ctx (R);
	FaugereLachartre<Ring> Solver (ctx);

	size_t rank;
	typename Ring::Element det;

	commentator.start("Solver.echelonize");
		Solver.echelonize (A, A, rank, det, reduced);
	commentator.stop("Solver.echelonize");

	MatrixUtils::show_mem_usage("");

	report << "True Rank : " << rank << endl;
}

int main (int argc, char **argv)
{
	const char *file_name = "";
	bool reduced = false;
	bool use_uint16 = false;
	bool use_double = !use_uint16;

	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &file_name},
		{ 'r', "-r", "compute reduced echelon form?", TYPE_NONE, &reduced},
		{ 'u', "-u", "compute over uint16 ring", TYPE_NONE, &use_uint16},
		{ 'd', "-d", "compute over bouble ring", TYPE_NONE, &use_double},
		{ '\0' }
	};

	use_double = !use_uint16;

	parseArguments (argc, argv, args, "", 0);

	commentator.setReportStream (std::cout);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtils::loadF4Modulus(file_name);
	Ring R (modulus);

	Context<Ring> ctx (R);

	SparseMatrix<Ring::Element> A;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);


#ifdef __LELA_BLAS_AVAILABLE
	report << "BLAS IS WORKING" << endl;
#else
	report << "BLAS IS NOT WORKING" << endl;
#endif

#ifdef	__LELA_HAVE_CBLAS
	report << "__LELA_HAVE_CBLAS" << endl;
#endif

	commentator.start("Loading matrix");
		MatrixUtils::loadF4Matrix__low_memory_syscall_no_checks(R, A, file_name);
	commentator.stop(MSG_DONE);
	report << endl;
	SHOW_MATRIX_INFO_SPARSE(A);

	MatrixUtils::show_mem_usage("Loading matrix");
	report << endl;

	if(use_uint16)
		test_LELA_FGL_uint16_elts(R, A, reduced);
	else if(use_double)
		test_LELA_FGL_double_elts(R, A, reduced);


	return 0;
}
