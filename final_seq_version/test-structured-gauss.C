/*
 * test-structured-gauss.C
 *
 *  Created on: 24 mai 2012
 *      Author: martani  (LIP6 / UPMC University Paris06)
 * 
 * Tests the structured Gaussian elmination method
 * run with; ./test-structured-gauss - -f path/to/matrix/file
 */

#include "consts-macros.h"
#include "types.h"
#include "matrix-utils.h"
//#include "level3-ops.h"
//#include "indexer.h"

#include "lela/matrix/sparse.h"
#include "lela/matrix/dense.h"
#include "lela/util/commentator.h"
#include "../util/support.h"

#include "lela/algorithms/elimination.h"
#include "lela/algorithms/gauss-jordan.h"

#include "structured-gauss-lib.h"


using namespace LELA;
using namespace std;

template <typename Ring, typename Matrix>
void test_structured_gauss_standard(const Ring& R, Matrix& A)
{
	size_t rank;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Structured Gauss Rref");
		rank = StructuredGauss::echelonize_reduced(R, A);
		report << "True rank " << rank << endl;

	commentator.stop(MSG_DONE);
}


template <typename Ring, typename Matrix>
void test_structured_gauss_acc64(const Ring& R, Matrix& A)
{
	size_t rank;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Structured Gauss Rref Accumulator 64");
		rank = StructuredGauss::echelonize_reduced_uint16(R, A);
		report << "True rank " << rank << endl;

	commentator.stop(MSG_DONE);
}

template <typename Ring, typename Matrix>
void test_RREF_LELA(const Ring& R, Matrix& A)
{
	Context<Ring> ctx (R);
	Elimination<Ring> elim (ctx);

	size_t rank;
	typename Ring::Element det;

	DenseMatrix<typename Ring::Element> L (A.rowdim (), A.rowdim ());
	typename GaussJordan<Ring>::Permutation P;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("LELA echelonize_reduced");
		elim.echelonize_reduced (A, L, P, rank, det);
		report << "True rank " << rank << endl;

	commentator.stop(MSG_DONE);
}

int main (int argc, char **argv)
{
	char *file_name = NULL;

	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &file_name},
		{ '\0' }
	};

	parseArguments (argc, argv, args, "", 0);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
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
	
	commentator.start("Loading matrix");
		MatrixUtils::loadF4Matrix__low_memory_syscall_no_checks(R, A, file_name);
	commentator.stop(MSG_DONE);
	report << endl;
	SHOW_MATRIX_INFO_SPARSE(A);

	MatrixUtils::show_mem_usage("Loading matrix");
	report << endl;
	/*commentator.start("Writing matrix to file");
		report << "\t" << getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-orig.pbm")).c_str() << endl;
		dumpMatrixAsPbmImage(A, getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-orig.pbm")).c_str());
	commentator.stop(MSG_DONE);
	report << endl;*/

	SparseMatrix<Ring::Element> C (A.rowdim (), A.coldim ()), D (A.rowdim (), A.coldim ());

	BLAS3::copy(ctx, A, C);
	BLAS3::copy(ctx, A, D);

	test_structured_gauss_acc64(R, A);

	//test_structured_gauss_standard(R, D);

	test_RREF_LELA(R, C);

	if(BLAS3::equal(ctx, A, C))
		report << "MATRIX RREF OK" << endl;
	else
		report << "MATRIX RREF WRONG" << endl;

//	if(BLAS3::equal(ctx, D, C))
//		report << "MATRIX RREF OK" << endl;
//	else
//		report << "MATRIX RREF WRONG" << endl;

//	commentator.start("Structured Gauss Rref", "STRUCTURED GAUSS");
//		rank = StructuredGauss::echelonize_reduced(R, A);
//		report << "Rank " << rank << endl;
//	commentator.stop(MSG_DONE, "STRUCTURED GAUSS");
	report << endl;
	SHOW_MATRIX_INFO_SPARSE(A);

	MatrixUtils::show_mem_usage("STRUCTURED GAUSS");



	/*commentator.start("Writing matrix to file");
		report << "\t" << getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-structured_rref.pbm")).c_str() << endl;
		dumpMatrixAsPbmImage(A, getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-structured_rref.pbm")).c_str());
	commentator.stop(MSG_DONE);*/


	return 0;
}

