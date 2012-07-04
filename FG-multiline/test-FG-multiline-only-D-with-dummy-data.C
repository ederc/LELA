/*
 * test-FG-multiline-only-D-with-dummy-data.C
 *
 *  Created on: 25 juin 2012
 *      Author: martani
 */




/*
 * test-FG-multiline-only-D.C
 *
 *  Created on: 25 juin 2012
 *      Author: martani
 */



#include <iostream>

#include "FG-types.h"
/*#include "../new-faugere-lachartre/matrix-util.h"
#include "../new-faugere-lachartre/matrix-op.h"
#include "../new-faugere-lachartre/indexer.h"*/
#include "../only-D/matrix-util.h"
#include "../only-D/matrix-op.h"
#include "../only-D/indexer.h"

#include "multiline-indexer.h"
#include "matrix-util-m.C"
#include "matrix-op-m.C"

#include "../util/support.h"


using namespace LELA;
using namespace std;

bool testFaugereLachartre(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	//modulus = 7;
	Modular<modulus_type> R(modulus);
	Context<Ring> ctx(R);

	MultiLineIndexer multilineIndexer;
	Indexer<uint32> simpleIndexer;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	/*SparseMatrix<Ring::Element> A (4, 5);
	A[0].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(1,1));
	A[0].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(3,1));
	A[0].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(4,2));

	A[1].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(0,1));
	A[1].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(1,2));
	A[1].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(3,2));
	A[1].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(4,3));

	A[2].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(1,3));
	A[2].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(3,3));
	A[2].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(4,4));

	A[3].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(1,2));
	A[3].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(2,5));
	A[3].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(3,4));
	A[3].push_back(typename SparseMatrix<Ring::Element>::Row::value_type(4,5));*/
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");


commentator.start("GAUSS ROUND 1", "GAUSS_ROUND_1");
	commentator.start("[MultiLineIndexer] constructing indexes");
		multilineIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);

	commentator.start("[SimpleIndexer] constructing indexes");
		simpleIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);
	//MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	/*SparseMultilineMatrix<Ring::Element> sub_A(2, 5), // (multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(2, 2), // (multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 //sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(2, 2),// (A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(2, 5);*/
	SparseMultilineMatrix<Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);




	//SparseMatrix<Ring::Element> sub_D_copy(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

	commentator.start("[multilineIndexer] constructing sub matrices");
		multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);
	report << endl;

	/*sub_A[0].push_back(0, 1, 0);
	sub_A[0].push_back(1, 2, 1);
	sub_A[0].push_back(2, 3, 3);
	sub_A[0].push_back(3, 4, 2);
	sub_A[0].push_back(4, 5, 4);

	sub_C[0].push_back(0, 1, 2);
	sub_C[0].push_back(1, 4, 1);
	sub_C[0].push_back(2, 3, 3);
	sub_C[0].push_back(3, 2, 1);
	sub_C[0].push_back(4, 1, 4);

	sub_B[0].push_back(0, 1, 4);
	sub_B[0].push_back(1, 4, 0);

	sub_D[0].push_back(0, 6, 0);
	sub_D[0].push_back(1, 2, 3);

	NS::write(sub_A);
	NS::write(sub_B);
	NS::write(sub_C);
	NS::write(sub_D);*/

////////////////////////////////////////////////////////////////////////////////////////////////////
	SparseMatrix<Ring::Element> sub_A_(simpleIndexer.Npiv, simpleIndexer.Npiv),
								sub_B_(simpleIndexer.Npiv, A.coldim() - simpleIndexer.Npiv),
								sub_C_(A.rowdim() - simpleIndexer.Npiv, simpleIndexer.Npiv),
								sub_D_(A.rowdim() - simpleIndexer.Npiv, A.coldim() - simpleIndexer.Npiv),
								sub_D_copy(A.rowdim() - simpleIndexer.Npiv, A.coldim() - simpleIndexer.Npiv);

	commentator.start("[SimpleIndexer] constructing sub matrices");
		simpleIndexer.constructSubMatrices(A, sub_A_, sub_B_, sub_C_, sub_D_, false);
	commentator.stop(MSG_DONE);

////////////////////////////////////////////////////////////////////////////////////////////////////

	/*if (NS::equal(R, sub_A, sub_A_)
			&& NS::equal(R, sub_B, sub_B_)
			&& NS::equal(R, sub_C, sub_C_)
			&& NS::equal(R, sub_D, sub_D_))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;

	if (NS::equal_reverse(R, sub_A, sub_A_)
			&& NS::equal_reverse(R, sub_B, sub_B_)
			&& NS::equal_reverse(R, sub_C, sub_C_)
			&& NS::equal_reverse(R, sub_D, sub_D_))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;*/

	//dumpMatrixAsPbmImage(sub_A, "sub_A_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_A_, "sub_A_simple.pbm");

	//dumpMatrixAsPbmImage(sub_C, "sub_C_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_C_, "sub_C_simple.pbm");

	/*commentator.start("B = A^-1 x B [multiline]");
		NS::reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop(MSG_DONE);*/

	/*commentator.start("B = A^-1 x B [simple]");
		MatrixOp::reducePivotsByPivots(R, sub_A_, sub_B_, false, 8);
	commentator.stop(MSG_DONE);

	if (NS::equal(R, sub_B, sub_B_))
		cout << "sub_B EQUAL" << endl;
	else
		cout << "sub_B NOT EQUAL" << endl;*/

	/*commentator.start("D = D - CB [multiline]");
		NS::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
	commentator.stop(MSG_DONE);*/

	/*commentator.start("D = D - CB [Simple]");
		MatrixOp::reduceNonPivotsByPivots(R, sub_C_, sub_B_, sub_D_);
	commentator.stop(MSG_DONE);

	if (NS::equal(R, sub_D, sub_D_))
		cout << "sub_D EQUAL" << endl;
	else
		cout << "sub_D NOT EQUAL" << endl;*/

	commentator.start("Reducing CD", "REDUCE_CD");
		commentator.start("D <- D - CB [MODIFIED]");
			NS::reduceCD(R, sub_A, sub_B, sub_C, sub_D);
		commentator.stop(MSG_DONE);

		//NS::write(sub_C);

		commentator.start("Copy D");
			NS::copy(R, sub_D, sub_D_copy);
		commentator.stop(MSG_DONE);

		//BLAS3::write(ctx, report, sub_D_copy, FORMAT_SAGE);
		/*if (NS::equal(R, sub_D, sub_D_copy))
			cout << "COPY EQUAL" << endl;
		else
			cout << "COPY NOT EQUAL" << endl;*/
		size_t rank;
		commentator.start("Echelonize (D)");
			rank = MatrixOp::echelonize(R, sub_D_copy);
			MatrixUtil::makeRowsUnitary(R, sub_D_copy);		//D pivots are not necessarily unitary
		commentator.stop(MSG_DONE);
		report << "Rank of D: " << rank << endl;

		//BLAS3::write(ctx, report, sub_D_copy, FORMAT_SAGE);

	commentator.stop("REDUCE_CD");


	commentator.start("[SIMPLE] reduce CD", "SIMPLE_REDUCE_CD");
		commentator.start("D <- D - CB [MODIFIED]");
			MatrixOp::reduceCD(R, sub_A_, sub_B_, sub_C_, sub_D_);
		commentator.stop(MSG_DONE);

		commentator.start("Echelonize (D)");
			rank = MatrixOp::echelonize(R, sub_D_);
			MatrixUtil::makeRowsUnitary(R, sub_D_);		//D pivots are not necessarily unitary
		commentator.stop(MSG_DONE);
		report << "Rank of D: " << rank << endl;

	commentator.stop("SIMPLE_REDUCE_CD");

commentator.stop("GAUSS_ROUND_1");

	if (BLAS3::equal(ctx, sub_D_, sub_D_copy))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;

	report << endl << "True Rank " << multilineIndexer.Npiv + rank << endl;
	//dumpMatrixAsPbmImage(sub_B, "sub_B_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_B_, "sub_B_simple.pbm");

	return true;
}

int main(int argc, char **argv)
{
	char *fileName = NULL;

	bool pass = true;

	static Argument args[] =
	{
	{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING,
			&fileName },
	{ '\0' } };

	parseArguments(argc, argv, args, "", 0);

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(
			Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDepth(3);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDetailLevel(
			Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(PROGRESS_REPORT).setMaxDepth(3);

	commentator.start("Faugère-Lachartre Bloc Version", "FaugereLachartre");

	pass = testFaugereLachartre(fileName);

	commentator.stop(MSG_STATUS (pass));

	return pass ? 0 : -1;
}

