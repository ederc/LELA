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
#include "../only-D/structured-gauss-lib.h"

#include "multiline-indexer.h"
#include "matrix-util-m.C"
#include "matrix-op-m.C"

#include "../util/support.h"


using namespace LELA;
using namespace std;

bool testFaugereLachartre(const char *file_name, bool parallel, bool verify_results)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R(modulus);
	Context<Ring> ctx(R);
	size_t rank;

	MultiLineIndexer multilineIndexer;
	Indexer<uint32> simpleIndexer;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	SparseMatrix<Ring::Element> C (A.rowdim(), A.coldim());
	SparseMatrix<Ring::Element> M_orig (A.rowdim(), A.coldim());

	BLAS3::copy(ctx, A, C);
	BLAS3::copy(ctx, A, M_orig);

	commentator.start("Reducing CD", "REDUCE_CD");

		commentator.start("[MultiLineIndexer] constructing indexes");
			multilineIndexer.processMatrix(A);
		commentator.stop(MSG_DONE);

		report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

		SparseMultilineMatrix<Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
											 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
											 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
											 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

		report << "sub_A (" << sub_A.rowdim() << ", " << sub_A.coldim() << ")" << endl;
		report << "sub_B (" << sub_B.rowdim() << ", " << sub_B.coldim() << ")" << endl;
		report << "sub_C (" << sub_C.rowdim() << ", " << sub_C.coldim() << ")" << endl;
		report << "sub_D (" << sub_D.rowdim() << ", " << sub_D.coldim() << ")" << endl;

		/*sub_D(5, 6);

		sub_D[0].push_back(0, 0, 1);
		sub_D[0].push_back(1, 1, 2);
		sub_D[0].push_back(3, 1, 2);
		sub_D[0].push_back(4, 2, 3);
		sub_D[0].push_back(5, 6, 1);

		sub_D[1].push_back(1, 3, 2);
		sub_D[1].push_back(2, 0, 5);
		sub_D[1].push_back(3, 3, 4);
		sub_D[1].push_back(4, 4, 5);
		sub_D[1].push_back(5, 2, 0);

		sub_D[2].push_back(0, 1, 0);
		sub_D[2].push_back(1, 2, 0);
		sub_D[2].push_back(2, 3, 0);
		sub_D[2].push_back(3, 4, 0);
		sub_D[2].push_back(4, 5, 0);
		sub_D[2].push_back(5, 6, 0);*/


		SparseMatrix<Ring::Element>  sub_D_copy(sub_D.rowdim(), sub_D.coldim());

		commentator.start("[multilineIndexer] constructing sub matrices");
			multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
		commentator.stop(MSG_DONE);
		report << endl;


		/*SparseMatrix<Ring::Element> sub_A__(sub_A.rowdim(), sub_A.coldim()),
									sub_B__(sub_B.rowdim(), sub_B.coldim ());

		NS::copy(R, sub_A, sub_A__);
		NS::copy(R, sub_B, sub_B__);

		commentator.start("[Multiline] B = A^-1 B");
			NS::reducePivotsByPivots(R, sub_A, sub_B);
		commentator.stop(MSG_DONE);


		commentator.start("[Simple] B = A^-1 B");
			MatrixOp::reducePivotsByPivots(R, sub_A__, sub_B__);
		commentator.stop(MSG_DONE);

		if(NS::equal(R, sub_B, sub_B__, false))
			cout << "OK" << endl;
		else
			cout << "WRONG" << endl;*/

		commentator.start("D <- D - CB [MODIFIED]");
			if(parallel)
				NS::reduceCDParallel(R, sub_A, sub_B, sub_C, sub_D, 8);
			else
				NS::reduceCD(R, sub_A, sub_B, sub_C, sub_D);
		commentator.stop(MSG_DONE);

		/*commentator.start("Copy D");
			NS::copy(R, sub_D, sub_D_copy);
		commentator.stop(MSG_DONE);*/

		//NS::write(sub_D);

		commentator.start("[Multiline] Echelonize (D)");
			rank = NS::echelonize(R, sub_D);
		commentator.stop(MSG_DONE);
		report << "Rank of D: " << rank << endl;
		/*commentator.start("Copy D");
			NS::copy(R, sub_D, sub_D_copy);
		commentator.stop(MSG_DONE);
		MatrixUtil::dumpMatrixAsPbmImage(sub_D_copy, "copy.pbm");*/

		/*commentator.start("Echelonize (D)");
			rank = MatrixOp::echelonize(R, sub_D_copy);
			MatrixUtil::makeRowsUnitary(R, sub_D_copy);		//D pivots are not necessarily unitary
		commentator.stop(MSG_DONE);*/
		/*NS::write(sub_D);
		BLAS3::write(ctx, report, sub_D_copy, FORMAT_SAGE);*/

		MultiLineIndexer multiIdxrCD;
		/*commentator.start("Processing new matrix D");
			multiIdxrCD.processMatrix(sub_D_copy);
		commentator.stop(MSG_DONE);*/

		commentator.start("[Multiline] Processing new matrix D");
			multiIdxrCD.processMatrix(sub_D);
		commentator.stop(MSG_DONE);
		report << "Pivots found: " << multiIdxrCD.Npiv << endl << endl;

		commentator.start("[Multiline] Combine inner indexer");
			multilineIndexer.combineInnerIndexer(multiIdxrCD, true);
		commentator.stop(MSG_DONE);

		commentator.start("[Multiline] Reconstructing matrix");
			multilineIndexer.reconstructMatrix(A, sub_A, sub_B, sub_D);
		commentator.stop(MSG_DONE);

		//MatrixUtil::dumpMatrixAsPbmImage(A, "A_gauss.pbm");

		report << "------------------------------------------------" << endl;
		MultiLineIndexer idx2;
		commentator.start("[MultiLineIndexer] constructing indexes");
			idx2.processMatrix(A);
		commentator.stop(MSG_DONE);

		report << "Pivots found: " << idx2.Npiv << endl << endl;

		SparseMultilineMatrix<Ring::Element> sub_A_prime(idx2.Npiv, idx2.Npiv),
											 sub_B_prime(idx2.Npiv, A.coldim() - idx2.Npiv),
											 sub_C_prime(A.rowdim() - idx2.Npiv, idx2.Npiv),
											 sub_D_prime(A.rowdim() - idx2.Npiv, A.coldim() - idx2.Npiv);

		report << "sub_A (" << sub_A_prime.rowdim() << ", " << sub_A_prime.coldim() << ")" << endl;
		report << "sub_B (" << sub_B_prime.rowdim() << ", " << sub_B_prime.coldim() << ")" << endl;
		report << "sub_C (" << sub_C_prime.rowdim() << ", " << sub_C_prime.coldim() << ")" << endl;
		report << "sub_D (" << sub_D_prime.rowdim() << ", " << sub_D_prime.coldim() << ")" << endl;


		commentator.start("[multilineIndexer] constructing sub matrices 2 ");
			idx2.constructSubMatrices(A, sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime, true);
		commentator.stop(MSG_DONE);

		commentator.start("[Multiline] B = A^-1 B");
			NS::reducePivotsByPivots(R, sub_A_prime, sub_B_prime);
		commentator.stop(MSG_DONE);

		MultiLineIndexer inner_dummy_idxr;
		inner_dummy_idxr.processMatrix(sub_D_prime);

		idx2.combineInnerIndexer(inner_dummy_idxr, true);

		commentator.start("[Multiline] Reconstructing final matrix");
			idx2.reconstructMatrix(A, sub_B_prime);
		commentator.stop(MSG_DONE);

		//MatrixUtil::dumpMatrixAsPbmImage(A, "A_rref.pbm");

	commentator.stop("REDUCE_CD");

	/*commentator.start("Structured rref Original matrix", "STRUCTURED_RREF");
		StructuredGauss::echelonize_reduced(R, M_orig);
	commentator.stop(MSG_DONE, "STRUCTURED_RREF");
	MatrixUtil::dumpMatrixAsPbmImage(M_orig, "M0_rref.pbm");

	if(BLAS3::equal(ctx, M_orig, A))
	{
		report << "A [Multiline] - Result CORRECT" << std::endl;
		report << endl;
	}
	else
	{
		report << "A [Multiline] - Result NOT OK" << std::endl;
		report << endl;
	}

	return true;*/

	report << endl << "*******************************************" << endl << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////
	commentator.start("[SIMPLE] reduce CD", "SIMPLE_REDUCE_CD");

		commentator.start("[SimpleIndexer] constructing indexes");
			simpleIndexer.processMatrix(C);
		commentator.stop(MSG_DONE);
		report << "Pivots found: " << simpleIndexer.Npiv << endl << endl;

		SparseMatrix<Ring::Element> sub_A_(simpleIndexer.Npiv, simpleIndexer.Npiv),
									sub_B_(simpleIndexer.Npiv, C.coldim() - simpleIndexer.Npiv),
									sub_C_(C.rowdim() - simpleIndexer.Npiv, simpleIndexer.Npiv),
									sub_D_(C.rowdim() - simpleIndexer.Npiv, C.coldim() - simpleIndexer.Npiv);

		commentator.start("[SimpleIndexer] constructing sub matrices");
			simpleIndexer.constructSubMatrices(C, sub_A_, sub_B_, sub_C_, sub_D_, true);
		commentator.stop(MSG_DONE);

		commentator.start("D <- D - CB [MODIFIED]");
			MatrixOp::reduceCD(R, sub_A_, sub_B_, sub_C_, sub_D_);
		commentator.stop(MSG_DONE);

		commentator.start("Echelonize (D)");
			rank = MatrixOp::echelonize(R, sub_D_);
			MatrixUtil::makeRowsUnitary(R, sub_D_);		//D pivots are not necessarily unitary
		commentator.stop(MSG_DONE);
		report << "Rank of D: " << rank << endl;

		Indexer<uint32> simpleIdxrCD;

		commentator.start("Processing new matrix D");
			simpleIdxrCD.processMatrix(sub_D_);
		commentator.stop(MSG_DONE);
		report << "Pivots found: " << simpleIdxrCD.Npiv << endl << endl;


		commentator.start("Combine inner indexer");
			simpleIndexer.combineInnerIndexer(simpleIdxrCD, true);
		commentator.stop(MSG_DONE);

		commentator.start("Reconstructing matrix");
			simpleIndexer.reconstructMatrix(C, sub_A_, sub_B_, sub_D_);
		commentator.stop(MSG_DONE);

	commentator.stop("SIMPLE_REDUCE_CD");


///////////////////////////////////////////////////////////////////////////////////////////////////
	bool pass = true;

	if(verify_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Checks on intermediary matrices" << endl;

		/*if(multilineIndexer.compare(simpleIndexer))
			report << "Indexer multiline - simple OK" << endl;
		else
		{
			report << "Indexer multiline - simple NOT OK" << endl;
			pass = false;
		}*/

		//NS::dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");
		//MatrixUtil::dumpMatrixAsPbmImage(sub_A_, "sub_A_.pbm");


		if(NS::equal(R, sub_A, sub_A_, false))
			report << "sub_A ---> CORRECT" << endl;
		else
		{
			report << "sub_A ---> NOT CORRECT" << endl;
			pass = false;
		}

		if(NS::equal(R, sub_B, sub_B_, false))
			report << "sub_B ---> CORRECT" << endl;
		else
		{
			report << "sub_B ---> NOT CORRECT" << endl;
			pass = false;
		}

		if(NS::equal(R, sub_C, sub_C_, false))
			report << "sub_C ---> CORRECT" << endl;
		else
		{
			report << "sub_C ---> NOT CORRECT" << endl;
			pass = false;
		}

		StructuredGauss::echelonize_reduced(R, sub_D_);
		NS::copy(R, sub_D, sub_D_copy);
		StructuredGauss::echelonize_reduced(R, sub_D_copy);

		//MatrixUtil::dumpMatrixAsPbmImage(sub_D_, "sub_D.pbm");
		//MatrixUtil::dumpMatrixAsPbmImage(sub_D_copy, "sub_D_multi.pbm");

		/*if(NS::equal(R, sub_D, sub_D_))*/
		if (BLAS3::equal(ctx, sub_D_, sub_D_copy))
			report << "RREF(sub_D) ---> CORRECT" << endl;
		else
		{
			report << "RREF(sub_D) ---> NOT CORRECT" << endl;
			pass = false;
		}

		report << endl;

		//if(!pass)
			//return false;

///////////////////////////////////////////////////////////////////////////////////////////////////
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Computing reduced echelon form of the matrix using structured Gaussian elimination" << endl;

		//MatrixUtil::dumpMatrixAsPbmImage(A, "A.pbm");
		//MatrixUtil::dumpMatrixAsPbmImage(C, "C.pbm");

		commentator.start("Structured rref Original matrix", "STRUCTURED_RREF");
			StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		commentator.start("Structured rref A [Multiline]", "STRUCTURED_RREF");
			StructuredGauss::echelonize_reduced(R, C);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		commentator.start("Structured rref C [Simple]", "STRUCTURED_RREF");
			StructuredGauss::echelonize_reduced(R, A);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		report << endl;

		if(BLAS3::equal(ctx, M_orig, C))
		{
			report << "C [Simple] - Result CORRECT" << std::endl;
			report << endl;
		}
		else
		{
			report << "C [Simple] - Result NOT OK" << std::endl;
			report << endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, M_orig, A))
		{
			report << "A [Multiline] - Result CORRECT" << std::endl;
			report << endl;
		}
		else
		{
			report << "A [Multiline] - Result NOT OK" << std::endl;
			report << endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, C, A))
		{
			report << "A [Multiline] == C [Simple]" << std::endl;
			report << endl;
		}
		else
		{
			report << "A [Multiline] != C [Simple]" << std::endl;
			report << endl;
			pass = false;
		}
	}

	report << endl;

	return pass;
}

int main(int argc, char **argv)
{
	char *fileName = NULL;

	bool pass = true;
	bool validate_results = false;
	bool parallel = false;

	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ 'p', "-p", "Use parallel version", TYPE_NONE, &parallel },
		{ '\0' }
	};

	parseArguments(argc, argv, args, "", 0);

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(
			Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDepth(3);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDetailLevel(
			Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(PROGRESS_REPORT).setMaxDepth(3);

	commentator.start("Faugère-Lachartre Bloc Version");

	pass = testFaugereLachartre(fileName, parallel, validate_results);

	commentator.stop(MSG_STATUS (pass));

	return pass ? 0 : -1;
}

