/*
 * test-FG-multiline.C
 *
 *  Created on: 14 juin 2012
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

#define SHOW_MATRIX_INFO(M)			\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << "(" << right << setw(8) << M.rowdim() << "x " << left << setw(8) << M.coldim() << ")\t"	\
			<< "- density: " << MatrixUtil::getMatrixSizeAndDensity(M).second << "% \t*****" << endl;

template <typename Ring>
bool testFaugereLachartre_invert_steps(const Ring& R, SparseMatrix<typename Ring::Element>& A, bool only_D, bool parallel, bool validate_results)
{
	Context<Ring> ctx(R);

	MultiLineIndexer multilineIndexer;
	size_t rank;
	bool reduced = false;
	//std::pair<uint64, double> size_density;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	SparseMatrix<typename Ring::Element> M_orig (A.rowdim(), A.coldim());

	if(validate_results)
		BLAS3::copy(ctx, A, M_orig);

commentator.start("FG_LACHARTRE_INVERT_STEPS", "FG_LACHARTRE_INVERT_STEPS");
commentator.start("ROUND 1", "ROUND 1");

	commentator.start("[MultiLineIndexer] constructing indexes");
		multilineIndexer.processMatrix(A);
	commentator.stop("[MultiLineIndexer] constructing indexes");

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	SparseMultilineMatrix<typename Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

	commentator.start("[multilineIndexer] constructing sub matrices");
		multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
	commentator.stop("[multilineIndexer] constructing sub matrices"); report << endl;

	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);
	report << endl;

	commentator.start("D <- D - CB [MODIFIED]");
		if(parallel)
			NS::reduceCDParallel(R, sub_A, sub_B, sub_C, sub_D, 8);
		else
		{
			NS::reduceCD(R, sub_A, sub_B, sub_C, sub_D);
			/*commentator.start("Reduce C by A");
				NS::reduceCD_onlyC(R, sub_A, sub_C);
			commentator.stop(MSG_DONE);

			commentator.start("D <- D - C*B");
				NS::reduceCD_onlyD(R, sub_C, sub_B, sub_D);
			commentator.stop(MSG_DONE);*/
		}
	commentator.stop("D <- D - CB [MODIFIED]");

	SHOW_MATRIX_INFO(sub_C); report << endl;
	//NS::dumpMatrixAsPbmImage(sub_D, "sub_D_CD.pbm");

	commentator.start("[Multiline] Echelonize (D)");
		rank = NS::echelonize(R, sub_D);
	commentator.stop("[Multiline] Echelonize (D)"); report << endl;
	report << "Rank of D: " << rank << endl; report << endl;


	MultiLineIndexer multiIdxrCD;

	commentator.start("[Multiline] Processing new matrix D");
		multiIdxrCD.processMatrix(sub_D);
	commentator.stop(MSG_DONE); report << endl;
	//report << "Pivots found: " << multiIdxrCD.Npiv << endl << endl;

	commentator.start("[Multiline] Combine inner indexer");
		multilineIndexer.combineInnerIndexer(multiIdxrCD, true);
	commentator.stop(MSG_DONE); report << endl;

	commentator.start("[Multiline] Reconstructing matrix");
		multilineIndexer.reconstructMatrix(A, sub_A, sub_B, sub_D);
	commentator.stop(MSG_DONE); report << endl;

commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "------------------------------------------------" << endl;
commentator.start("ROUND 2");
		MultiLineIndexer idx2;
		commentator.start("[MultiLineIndexer] constructing indexes");
			idx2.processMatrix(A);
		commentator.stop(MSG_DONE); report << endl;

		report << "Pivots found: " << idx2.Npiv << endl << endl;

		SparseMultilineMatrix<typename Ring::Element> sub_A_prime(idx2.Npiv, idx2.Npiv),
											 sub_B_prime(idx2.Npiv, A.coldim() - idx2.Npiv),
											 sub_C_prime(A.rowdim() - idx2.Npiv, idx2.Npiv),
											 sub_D_prime(A.rowdim() - idx2.Npiv, A.coldim() - idx2.Npiv);

		/*report << "sub_A (" << sub_A_prime.rowdim() << ", " << sub_A_prime.coldim() << ")" << endl;
		report << "sub_B (" << sub_B_prime.rowdim() << ", " << sub_B_prime.coldim() << ")" << endl;
		report << "sub_C (" << sub_C_prime.rowdim() << ", " << sub_C_prime.coldim() << ")" << endl;
		report << "sub_D (" << sub_D_prime.rowdim() << ", " << sub_D_prime.coldim() << ")" << endl;*/


		commentator.start("[multilineIndexer] constructing sub matrices 2 ");
			idx2.constructSubMatrices(A, sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime, true);
		commentator.stop(MSG_DONE); report << endl;

		commentator.start("[Multiline] B = A^-1 B");
			NS::reducePivotsByPivots(R, sub_A_prime, sub_B_prime);
		commentator.stop("[Multiline] B = A^-1 B"); report << endl;

		MultiLineIndexer inner_dummy_idxr;
		inner_dummy_idxr.processMatrix(sub_D_prime);

		idx2.combineInnerIndexer(inner_dummy_idxr, true);

		commentator.start("[Multiline] Reconstructing final matrix");
			idx2.reconstructMatrix(A, sub_B_prime);
		commentator.stop("[Multiline] Reconstructing final matrix"); report << endl;
		reduced = true;

commentator.stop("ROUND 2");
	}
commentator.stop("FG_LACHARTRE_INVERT_STEPS");

	report << endl << "True Rank " << multilineIndexer.Npiv + rank << endl;
	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Computing reduced echelon form of the matrix using structured Gaussian elimination" << endl;

		commentator.start("Structured rref of original matrix", "STRUCTURED_RREF");
			size_t rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			commentator.start("Structured rref of A", "STRUCTURED_RREF");
				StructuredGauss::echelonize_reduced(R, A);
			commentator.stop(MSG_DONE, "STRUCTURED_RREF");
		}

		report << endl << ">> Rank " << rank_strucutured_gauss << endl;

		report << endl;

		if(BLAS3::equal(ctx, M_orig, A))
		{
			report << "Result CORRECT" << std::endl;
			report << endl;
		}
		else
		{
			report << "Result NOT OK" << std::endl;
			report << endl;
			pass = false;
		}
	}

	report << endl;

	return pass;
}

template <typename Ring>
bool testFaugereLachartre_normal(const Ring& R, SparseMatrix<typename Ring::Element>& A, bool only_D, bool validate_results)
{
	Context<Ring> ctx(R);

	MultiLineIndexer multilineIndexer;
	size_t rank;
	bool reduced = false;
	//std::pair<uint64, double> size_density;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	SparseMatrix<typename Ring::Element> M_orig (A.rowdim(), A.coldim());
	if(validate_results)
		BLAS3::copy(ctx, A, M_orig);


commentator.start("FG_LACHARTRE", "FG_LACHARTRE");
commentator.start("ROUND 1", "ROUND 1");

	commentator.start("[MultiLineIndexer] constructing indexes");
		multilineIndexer.processMatrix(A);
	commentator.stop("[MultiLineIndexer] constructing indexes");

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	SparseMultilineMatrix<typename Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

	commentator.start("[multilineIndexer] constructing sub matrices");
		multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
	commentator.stop("[multilineIndexer] constructing sub matrices"); report << endl;

	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);
	report << endl;

	commentator.start("[Multiline] B = A^-1 B");
		NS::reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop("[Multiline] B = A^-1 B");

	SHOW_MATRIX_INFO(sub_B); report << endl;

	commentator.start("[Multiline] D = D - CB");
		NS::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
	commentator.stop("[Multiline] D = D - CB"); report << endl;

	//NS::dumpMatrixAsPbmImage(sub_D, "sub_D.pbm");
	
	commentator.start("[Multiline] Echelonize (D)");
		rank = NS::echelonize(R, sub_D);
	commentator.stop("[Multiline] Echelonize (D)"); report << endl;
	report << "Rank of D: " << rank << endl; report << endl;

	if(only_D)
	{
		MultiLineIndexer idx2;
		idx2.processMatrix(sub_D);
		multilineIndexer.combineInnerIndexer(idx2, true);

		commentator.start("[Multiline] Reconstructing matrix");
			multilineIndexer.reconstructMatrix(A, sub_B, sub_D);
		commentator.stop("[Multiline] Reconstructing matrix"); report << endl;
	}
commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "------------------------------------------------" << endl;
commentator.start("ROUND 2", "ROUND 2");

	SparseMatrix<typename Ring::Element> sub_D_sparse (sub_D.rowdim(), sub_D.coldim());
	commentator.start("Copy D");
		NS::copy(R, sub_D, sub_D_sparse);
	commentator.stop(MSG_DONE);

	MultiLineIndexer inner_indexer;
	commentator.start("[MultiLineIndexer] constructing indexes");
		inner_indexer.processMatrix(sub_D_sparse);
	commentator.stop(MSG_DONE); report << endl;

	SparseMultilineMatrix<typename Ring::Element> D1(inner_indexer.Npiv, inner_indexer.Npiv),
										 D2(inner_indexer.Npiv, sub_D.coldim() - inner_indexer.Npiv),
										 B1(sub_B.rowdim(), inner_indexer.Npiv),
										 B2(sub_B.rowdim(), sub_B.coldim() - inner_indexer.Npiv);

	commentator.start("[MultiLineIndexer] constructing submatrices B1, B1, D1, D2");
		inner_indexer.constructSubMatrices(sub_B, sub_D_sparse, B1, B2, D1, D2);
	commentator.stop("[MultiLineIndexer] constructing submatrices B1, B1, D1, D2"); report << endl;

	commentator.start("D2 = D1^-1 x D2");
		NS::reducePivotsByPivots(R, D1, D2);
	commentator.stop("D2 = D1^-1 x D2"); report << endl;

	commentator.start("B2 <- B2 - D2 D1");
		NS::reduceNonPivotsByPivots(R, B1, D2, B2);
	commentator.stop("B2 <- B2 - D2 D1"); report << endl;

	commentator.start("[MultiLineIndexer] Reconstructing indexes");
		multilineIndexer.combineInnerIndexer(inner_indexer);
	commentator.stop(MSG_DONE); report << endl;

	commentator.start("[Indexer] Reconstructing matrix");
		multilineIndexer.reconstructMatrix(A, B2, D2);
	commentator.stop("[Indexer] Reconstructing matrix"); report << endl;
	reduced = true;

commentator.stop("ROUND 2");
	}
commentator.stop("FG_LACHARTRE", "FG_LACHARTRE");

report << endl << ">> True Rank " << multilineIndexer.Npiv + rank << endl;
	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Computing reduced echelon form of the matrix using structured Gaussian elimination" << endl;

		size_t rank_strucutured_gauss;
		commentator.start("Structured rref Original matrix", "STRUCTURED_RREF");
			rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			commentator.start("Structured rref of A", "STRUCTURED_RREF");
				StructuredGauss::echelonize_reduced(R, A);
			commentator.stop(MSG_DONE, "STRUCTURED_RREF");
		}

		report << endl << ">> Rank " << rank_strucutured_gauss << endl;
		report << endl;

		//MatrixUtil::dumpMatrixAsPbmImage(A, "A-rref.pbm");
		//MatrixUtil::dumpMatrixAsPbmImage(M_orig, "M0-rref.pbm");

		if(BLAS3::equal(ctx, M_orig, A))
		{
			report << "Result CORRECT" << std::endl;
			report << endl;
		}
		else
		{
			report << "Result NOT OK" << std::endl;
			report << endl;
			pass = false;
		}
	}

	report << endl;

	return pass;
}

template <typename Ring>
bool testFaugereLachartre_invert_steps_separated_C_D(const Ring& R, SparseMatrix<typename Ring::Element>& A, bool only_D, bool validate_results)
{
	Context<Ring> ctx(R);

	MultiLineIndexer multilineIndexer;
	size_t rank;
	bool reduced = false;
	//std::pair<uint64, double> size_density;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	SparseMatrix<typename Ring::Element> M_orig (A.rowdim(), A.coldim());

	if(validate_results)
		BLAS3::copy(ctx, A, M_orig);

commentator.start("FG_LACHARTRE_INVERT_STEPS", "FG_LACHARTRE_INVERT_STEPS");
commentator.start("ROUND 1", "ROUND 1");

	commentator.start("[MultiLineIndexer] constructing indexes");
		multilineIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	SparseMultilineMatrix<typename Ring::Element> sub_A(multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_B(multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv),
										 sub_C(A.rowdim() - multilineIndexer.Npiv, multilineIndexer.Npiv),
										 sub_D(A.rowdim() - multilineIndexer.Npiv, A.coldim() - multilineIndexer.Npiv);

	commentator.start("[multilineIndexer] constructing sub matrices");
		multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
	commentator.stop("[multilineIndexer] constructing sub matrices"); report << endl;

	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);
	report << endl;

	commentator.start("CD STEP");
		commentator.start("NS::reduceCD_onlyC");
			NS::reduceCD_onlyC(R, sub_A, sub_C);
		commentator.stop("NS::reduceCD_onlyC");
		SHOW_MATRIX_INFO(sub_C); report << endl;

		commentator.start("[Multiline] D = D - CB");
			NS::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D, false);
		commentator.stop("[Multiline] D = D - CB"); report << endl;
	commentator.stop(MSG_DONE);

	//NS::dumpMatrixAsPbmImage(sub_D, "sub_D_C_then_D.pbm");

	commentator.start("[Multiline] Echelonize (D)");
		rank = NS::echelonize(R, sub_D);
	commentator.stop("[Multiline] Echelonize (D)"); report << endl;
	report << "Rank of D: " << rank << endl; report << endl;


	MultiLineIndexer multiIdxrCD;

	commentator.start("[Multiline] Processing new matrix D");
		multiIdxrCD.processMatrix(sub_D);
	commentator.stop(MSG_DONE); report << endl;
	//report << "Pivots found: " << multiIdxrCD.Npiv << endl << endl;

	commentator.start("[Multiline] Combine inner indexer");
		multilineIndexer.combineInnerIndexer(multiIdxrCD, true);
	commentator.stop(MSG_DONE); report << endl;

	commentator.start("[Multiline] Reconstructing matrix");
		multilineIndexer.reconstructMatrix(A, sub_A, sub_B, sub_D);
	commentator.stop("[Multiline] Reconstructing matrix"); report << endl;

commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "------------------------------------------------" << endl;
commentator.start("ROUND 2");
		MultiLineIndexer idx2;
		commentator.start("[MultiLineIndexer] constructing indexes");
			idx2.processMatrix(A);
		commentator.stop(MSG_DONE); report << endl;

		report << "Pivots found: " << idx2.Npiv << endl << endl;

		SparseMultilineMatrix<typename Ring::Element> sub_A_prime(idx2.Npiv, idx2.Npiv),
											 sub_B_prime(idx2.Npiv, A.coldim() - idx2.Npiv),
											 sub_C_prime(A.rowdim() - idx2.Npiv, idx2.Npiv),
											 sub_D_prime(A.rowdim() - idx2.Npiv, A.coldim() - idx2.Npiv);

		/*report << "sub_A (" << sub_A_prime.rowdim() << ", " << sub_A_prime.coldim() << ")" << endl;
		report << "sub_B (" << sub_B_prime.rowdim() << ", " << sub_B_prime.coldim() << ")" << endl;
		report << "sub_C (" << sub_C_prime.rowdim() << ", " << sub_C_prime.coldim() << ")" << endl;
		report << "sub_D (" << sub_D_prime.rowdim() << ", " << sub_D_prime.coldim() << ")" << endl;*/


		commentator.start("[multilineIndexer] constructing sub matrices 2 ");
			idx2.constructSubMatrices(A, sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime, true);
		commentator.stop("[multilineIndexer] constructing sub matrices 2 "); report << endl;

		commentator.start("[Multiline] B = A^-1 B");
			NS::reducePivotsByPivots(R, sub_A_prime, sub_B_prime);
		commentator.stop("[Multiline] B = A^-1 B"); report << endl;

		MultiLineIndexer inner_dummy_idxr;
		inner_dummy_idxr.processMatrix(sub_D_prime);

		idx2.combineInnerIndexer(inner_dummy_idxr, true);

		commentator.start("[Multiline] Reconstructing final matrix");
			idx2.reconstructMatrix(A, sub_B_prime);
		commentator.stop("[Multiline] Reconstructing final matrix"); report << endl;
		reduced = true;

commentator.stop("ROUND 2");
	}
commentator.stop("FG_LACHARTRE_INVERT_STEPS");

	report << endl << "True Rank " << multilineIndexer.Npiv + rank << endl;
	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Computing reduced echelon form of the matrix using structured Gaussian elimination" << endl;

		commentator.start("Structured rref of original matrix", "STRUCTURED_RREF");
			size_t rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			commentator.start("Structured rref of A", "STRUCTURED_RREF");
				StructuredGauss::echelonize_reduced(R, A);
			commentator.stop(MSG_DONE, "STRUCTURED_RREF");
		}

		report << endl << ">> Rank " << rank_strucutured_gauss << endl;

		report << endl;

		if(BLAS3::equal(ctx, M_orig, A))
		{
			report << "Result CORRECT" << std::endl;
			report << endl;
		}
		else
		{
			report << "Result NOT OK" << std::endl;
			report << endl;
			pass = false;
		}
	}

	report << endl;

	return pass;
}

int main(int argc, char **argv)
{
	const char *fileName = "";

	bool pass = true;
	bool validate_results = false;
	bool parallel = false;
	bool method = 0;
	bool compare = false;
	bool only_D = false;
	bool reduce_C = false;

	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ 'p', "-p", "Use parallel version", TYPE_NONE, &parallel },
		{ 'm', "-m", "Use the new method", TYPE_NONE, &method },
		{ 'c', "-c", "Compared methods", TYPE_NONE, &compare },
		{ 'd', "-d", "Compute only the Gauss elimination of the matrix", TYPE_NONE, &only_D },
		{ 'o', "-o", "Reduce only C", TYPE_NONE, &reduce_C },
		{ '\0' }
	};

	parseArguments(argc, argv, args, "", 0);

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDepth(3);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(PROGRESS_REPORT).setMaxDepth(3);

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Faugère-Lachartre Multiline Version");

	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(fileName);
	Modular<modulus_type> R(modulus);

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A;
		NS::loadF4Matrix__low_memory_syscall_no_checks(R, A, fileName);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");
	report << endl;

	if(reduce_C && compare)
	{
		commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(6);
		SparseMatrix<Ring::Element> C (A.rowdim(), A.coldim());
		Context<Ring> ctx(R);
		BLAS3::copy(ctx, A, C);

		pass = testFaugereLachartre_invert_steps(R, A, only_D, parallel, validate_results);
		pass &= testFaugereLachartre_invert_steps_separated_C_D(R, C, only_D, validate_results);

		commentator.stop(MSG_STATUS (pass));
		return pass ? 0 : -1;
	}

	if(reduce_C)
	{
		commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(6);
		pass = testFaugereLachartre_invert_steps_separated_C_D(R, A, only_D, validate_results);
		return pass ? 0 : -1;
	}

	if(compare)
	{
		commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
		SparseMatrix<Ring::Element> C (A.rowdim(), A.coldim());
		Context<Ring> ctx(R);
		BLAS3::copy(ctx, A, C);

		pass = testFaugereLachartre_invert_steps(R, A, only_D, parallel, validate_results);
		pass &= testFaugereLachartre_normal(R, C, only_D, validate_results);

		commentator.stop(MSG_STATUS (pass));
		return pass ? 0 : -1;
	}

	if(method)
		pass = testFaugereLachartre_invert_steps(R, A, only_D, parallel, validate_results);
	else
		pass = testFaugereLachartre_normal(R, A, only_D, validate_results);

	commentator.stop(MSG_STATUS (pass));

	return pass ? 0 : -1;
}

