/*
 * test-FG-multiline.C
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#include <iostream>

#include "FG-types.h"
#include "../new-faugere-lachartre/matrix-util.h"
#include "../new-faugere-lachartre/matrix-op.h"
#include "../new-faugere-lachartre/indexer.h"
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
	Modular<modulus_type> R(modulus);
	Context<Ring> ctx(R);

	MultiLineIndexer multilineIndexer;
	Indexer<uint32> simpleIndexer;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);

	commentator.start("Loading matrix");
	SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	commentator.start("[MultiLineIndexer] constructing indexes");
	multilineIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);

	commentator.start("[SimpleIndexer] constructing indexes");
	simpleIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);
	//MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	SparseMultilineMatrix<Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

	commentator.start("[multilineIndexer] constructing sub matrices");
	multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);
	report << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////
	SparseMatrix<Ring::Element> sub_A_(simpleIndexer.Npiv, simpleIndexer.Npiv),
			sub_B_(simpleIndexer.Npiv, A.coldim() - simpleIndexer.Npiv), sub_C_(
					A.rowdim() - simpleIndexer.Npiv, simpleIndexer.Npiv),
			sub_D_(A.rowdim() - simpleIndexer.Npiv,
					A.coldim() - simpleIndexer.Npiv);

	commentator.start("[SimpleIndexer] constructing sub matrices");
		simpleIndexer.constructSubMatrices(A, sub_A_, sub_B_, sub_C_, sub_D_, false);
	commentator.stop(MSG_DONE);

////////////////////////////////////////////////////////////////////////////////////////////////////

	if (NS::equal(R, sub_A, sub_A_)
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
		cout << "NOT EQUAL" << endl;

	//dumpMatrixAsPbmImage(sub_A, "sub_A_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_A_, "sub_A_simple.pbm");

	//dumpMatrixAsPbmImage(sub_C, "sub_C_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_C_, "sub_C_simple.pbm");

	commentator.start("B = A^-1 x B [multiline]");
		NS::reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop(MSG_DONE);

	commentator.start("B = A^-1 x B [simple]");
		MatrixOp::reducePivotsByPivots(R, sub_A_, sub_B_, false, 8);
	commentator.stop(MSG_DONE);

	if (NS::equal(R, sub_B, sub_B_))
		cout << "sub_B EQUAL" << endl;
	else
		cout << "sub_B NOT EQUAL" << endl;

	commentator.start("D = D - CB [multiline]");
		NS::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
	commentator.stop(MSG_DONE);

	commentator.start("D = D - CB [Simple]");
		MatrixOp::reduceNonPivotsByPivots(R, sub_C_, sub_B_, sub_D_);
	commentator.stop(MSG_DONE);

	if (NS::equal(R, sub_D, sub_D_))
		cout << "sub_D EQUAL" << endl;
	else
		cout << "sub_D NOT EQUAL" << endl;

	commentator.start("gauss(D) [Simple]");
		MatrixOp::echelonize(R, sub_D_, true);
	commentator.stop(MSG_DONE);

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

