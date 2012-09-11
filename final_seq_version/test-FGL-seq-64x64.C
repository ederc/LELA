/*
 * test-FGL-seq.C
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */

#include "consts-macros.h"
#include "types.h"
#include "matrix-utils.h"
#include "level3-ops.h"
#include "parallel/echelon.h"

#include "indexer.h"


#include "lela/matrix/sparse.h"
#include "lela/util/commentator.h"
#include "../util/support.h"

#include "structured-gauss-lib.h"

#include "../only-D/indexer.h"
#include "../only-D/matrix-util.h"
#include "../only-D/matrix-op.h"
//
//#include "matrix-op-m.C"


using namespace LELA;
using namespace std;

template <typename Ring>
bool testFaugereLachartre_old_method(const Ring& R,
		SparseMatrix<typename Ring::Element>& A, bool validate_results,
		bool only_D, bool free_memory_on_the_go,
		bool horizontal)
{
	Context<Ring> ctx (R);
	SequentialIndexer<typename Ring::Element, IndexType> outer_indexer;

	/************************************************************************************************/
	Indexer<uint32> simple_indexer;
	/************************************************************************************************/
	
	uint32 rank;
	bool reduced = false;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Index type " << __IndexType__ << endl;

	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > sub_A, sub_B, sub_C, sub_D;
	SparseMultilineMatrix<typename Ring::Element> sub_D_multiline;

	SparseMatrix<typename Ring::Element> M_orig;
	
	if(validate_results)
	{
		M_orig = SparseMatrix<typename Ring::Element> (A.rowdim(), A.coldim());
		BLAS3::copy(ctx, A, M_orig);
	}

	/************************************************************************************************/
	commentator.start("SIMPLE INDEXER");
	commentator.start("[Indexer] constructing indexes");
		simple_indexer.processMatrix(A);
	commentator.stop(MSG_DONE);
   	report << "Pivots found: " << simple_indexer.Npiv << endl << endl;

	SparseMatrix<typename Ring::Element>	sub_A_simple (simple_indexer.Npiv, simple_indexer.Npiv),
					sub_B_simple (simple_indexer.Npiv, A.coldim () - simple_indexer.Npiv),
					sub_C_simple (A.rowdim () - simple_indexer.Npiv, simple_indexer.Npiv),
					sub_D_simple (A.rowdim () - simple_indexer.Npiv, A.coldim () - simple_indexer.Npiv);

	commentator.start("[Indexer] constructing sub matrices");
		simple_indexer.constructSubMatrices(A, sub_A_simple, sub_B_simple, sub_C_simple, sub_D_simple, false);
	commentator.stop(MSG_DONE);
	commentator.stop("SIMPLE INDEXER");
	/************************************************************************************************/
	
commentator.start("FG_LACHARTRE", "FG_LACHARTRE");
commentator.start("ROUND 1", "ROUND 1");

	commentator.start("[Bloc] construting submatrices");
		outer_indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, free_memory_on_the_go);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("[construting submatrices]"); report << endl;
	report << "Pivots found: " << outer_indexer.Npiv << endl << endl;


	MatrixUtils::dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_B.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D.pbm");
	
	MatrixUtils::dumpMatrixAsPbmImage(sub_A_simple, "sub_A_simple.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_B_simple, "sub_B_simple.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_C_simple, "sub_C_simple.pbm");
	MatrixUtils::dumpMatrixAsPbmImage(sub_D_simple, "sub_D_simple.pbm");

	SHOW_MATRIX_INFO_BLOC(sub_A);
	SHOW_MATRIX_INFO_BLOC(sub_B);
	SHOW_MATRIX_INFO_BLOC(sub_C);
	SHOW_MATRIX_INFO_BLOC(sub_D);

	/************************************************************************************************/
	if(MatrixUtils::equal(R, sub_A, sub_A_simple))
		report << "sub_A OK" << endl;
	else
		report << "sub_A NOT OK" << endl;
	
	if(MatrixUtils::equal(R, sub_B, sub_B_simple))
		report << "sub_B OK" << endl;
	else
		report << "sub_B NOT OK" << endl;
	
	if(MatrixUtils::equal(R, sub_C, sub_C_simple))
		report << "sub_C OK" << endl;
	else
		report << "sub_C NOT OK" << endl;
	
	if(MatrixUtils::equal(R, sub_D, sub_D_simple))
		report << "sub_D OK" << endl;
	else
		report << "sub_D NOT OK" << endl;
	
	/************************************************************************************************/
	
	
	commentator.start("[Bloc] B = A^-1 B", "[B = A^-1 B]");
		Level3Ops::reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop("[B = A^-1 B]");
	SHOW_MATRIX_INFO_BLOC(sub_B);
	sub_A.free(true);
	MatrixUtils::show_mem_usage("[B = A^-1 B]"); report << endl;

	//MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_B.pbm");

	commentator.start("[Bloc] D = D - C*B", "[D = D - C*B]");
		if(horizontal)
			Level3Ops::reduceNonPivotsByPivots_horizontal(R, sub_C, sub_B, sub_D);
		else
			Level3Ops::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
	commentator.stop("[D = D - C*B]");
	SHOW_MATRIX_INFO_BLOC(sub_D);
	sub_C.free(true);
	MatrixUtils::show_mem_usage("[D = D - C*B]"); report << endl;

	//return true;

	commentator.start("echelonize D", "[echelonize D]");
		rank = Level3Ops::echelonize(R, sub_D, sub_D_multiline, free_memory_on_the_go);
	commentator.stop("[echelonize D]");
	MatrixUtils::show_mem_usage("[echelonize D]"); report << endl;
	report << "Rank of D " << rank << endl;

	if(only_D)
	{
		SequentialIndexer<typename Ring::Element, IndexType> inner_idxr;
		inner_idxr.processMatrix(sub_D_multiline);
		outer_indexer.combineInnerIndexer(inner_idxr, true);

		MatrixUtils::show_mem_usage("[After processing]"); report << endl;

		commentator.start("[Bloc] Reconstructing matrix", "[Reconstructing matrix]");
			outer_indexer.reconstructMatrix(A, sub_B, sub_D_multiline, true);
		commentator.stop("[Reconstructing matrix]"); report << endl;
		MatrixUtils::show_mem_usage("[Reconstructing matrix]"); report << endl;

		SHOW_MATRIX_INFO_SPARSE(A);

	}
commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "------------------------------------------------" << endl;
commentator.start("ROUND 2", "ROUND 2");

	SequentialIndexer<typename Ring::Element, IndexType> inner_indexer;

	commentator.start("[MultiLineIndexer] constructing indexes");
		inner_indexer.processMatrix(sub_D_multiline);
	commentator.stop(MSG_DONE);
	report << "Pivots found: " << inner_indexer.Npiv << endl << endl;

	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > D1, D2, B1, B2;

	commentator.start("[Bloc] constructing submatrices B1, B1, D1, D2");
		inner_indexer.constructSubMatrices(sub_B, sub_D_multiline, B1, B2, D1, D2, free_memory_on_the_go);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("[construting submatrices]"); report << endl;

	commentator.start("D2 = D1^-1 x D2");
		Level3Ops::reducePivotsByPivots(R, D1, D2);
	commentator.stop(MSG_DONE);
	D1.free (true);
	MatrixUtils::show_mem_usage("[D2 = D1^-1 x D2]"); report << endl;

	commentator.start("B2 <- B2 - D2 D1");
		Level3Ops::reduceNonPivotsByPivots(R, B1, D2, B2);
	commentator.stop(MSG_DONE);
	B1.free (true);
	MatrixUtils::show_mem_usage("[D2 = D1^-1 x D2]"); report << endl;

	commentator.start("[MultiLineIndexer] Reconstructing indexes");
		outer_indexer.combineInnerIndexer(inner_indexer);
	commentator.stop(MSG_DONE);

	commentator.start("[Indexer] Reconstructing matrix");
		outer_indexer.reconstructMatrix(A, B2, D2, free_memory_on_the_go);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("[Reconstructing matrix]"); report << endl;
	reduced = true;

commentator.stop("ROUND 2");
	}
commentator.stop("FG_LACHARTRE", "FG_LACHARTRE");

	report << "True rank " << outer_indexer.Npiv + rank << endl;
	
	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "<< Computing reduced echelon form of the matrix using structured Gaussian elimination >>" << endl;

		commentator.start("Structured rref of original matrix", "STRUCTURED_RREF");
			size_t rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			//FAST
			for(uint32 i=0; i<A.rowdim ()/2; ++i)
			{
				A[i].swap(A[A.rowdim () - 1 - i]);
			}

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
bool testFaugereLachartre_new_method(const Ring& R,
		SparseMatrix<typename Ring::Element>& A, bool validate_results,
		bool only_D, bool free_memory_on_the_go)
{
	Context<Ring> ctx (R);

	SequentialIndexer<typename Ring::Element, IndexType> outer_indexer;

	size_t rank;
	bool reduced = false;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Index type " << __IndexType__ << endl;

	SparseMatrix<typename Ring::Element> M_orig;

	if(validate_results)
	{
		M_orig = SparseMatrix<typename Ring::Element> (A.rowdim(), A.coldim());
		BLAS3::copy(ctx, A, M_orig);
	}

commentator.start("FGL BLOC NEW METHOD");
commentator.start("ROUND 1");
	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > sub_A, sub_B, sub_C, sub_D;
	SparseMultilineMatrix<uint16> sub_D_multiline;

	commentator.start("[Bloc] construting submatrices");
		outer_indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, free_memory_on_the_go);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("[construting submatrices]"); report << endl;
	report << "Pivots found: " << outer_indexer.Npiv << endl << endl;


// 	MatrixUtils::dumpMatrixAsPbmImage(sub_C_multiline, "sub_C_multiline.pbm");
	
	//MatrixUtils::copy(sub_C_multiline, sub_C);
	
	SHOW_MATRIX_INFO_BLOC(sub_A);
	SHOW_MATRIX_INFO_BLOC(sub_B);
	SHOW_MATRIX_INFO_BLOC(sub_C);
	SHOW_MATRIX_INFO_BLOC(sub_D);
	
	commentator.start("[Bloc] C = MatrixOps::reduceC", "[reduceC]");
		Level3Ops::reduceC(R, sub_A, sub_C);
	commentator.stop("[reduceC]");
	//SHOW_MATRIX_INFO_BLOC(sub_C);
	MatrixUtils::show_mem_usage("[reduceC]"); report << endl;
	
	commentator.start("[Bloc] D = D - C*B", "[D = D - C*B]");
		Level3Ops::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D, false);
	commentator.stop("[D = D - C*B]");
	SHOW_MATRIX_INFO_BLOC(sub_D);
	sub_C.free(true);
	MatrixUtils::show_mem_usage("[D = D - C*B]"); report << endl;

	commentator.start("[Bloc] echelonize", "[echelonize D]");
//		if(parallel)
//			rank = Level3ParallelEchelon::echelonize__Parallel(R, sub_D, sub_D_multiline, free_memory_on_the_go, NUM_THREADS);
//		else
			rank = Level3Ops::echelonize(R, sub_D, sub_D_multiline, free_memory_on_the_go);
	commentator.stop("[echelonize D]");
	MatrixUtils::show_mem_usage("[echelonize D]"); report << endl;
	report << "Rank of D " << rank << endl;
	report << endl;

	SequentialIndexer<typename Ring::Element, IndexType> inner_idxr;
	commentator.start("[Bloc] Processing new matrix D");
		inner_idxr.processMatrix(sub_D_multiline);
	commentator.stop(MSG_DONE);
	report << "Pivots found: " << inner_idxr.Npiv << endl << endl;

	commentator.start("[Bloc] Combine inner indexer");
		outer_indexer.combineInnerIndexer(inner_idxr, true);
	commentator.stop(MSG_DONE); report << endl;

	commentator.start("[Bloc] Reconstructing matrix", "Reconstructing matrix]");
		outer_indexer.reconstructMatrix(A, sub_A, sub_B, sub_D_multiline, free_memory_on_the_go);
	commentator.stop("[Reconstructing matrix]");
	MatrixUtils::show_mem_usage("[Reconstructing matrix]"); report << endl;

	report << "True rank " << outer_indexer.Npiv + rank << endl;
	outer_indexer.freeMemory (); //not used anymore
commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "----------------------------------------------------------------------------------------" << endl;
commentator.start("ROUND 2");
		SequentialIndexer<typename Ring::Element, IndexType> idxr2;

		SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> >
				sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime;

		SparseMatrix<typename Ring::Element> dummySparse;
		SparseMultilineMatrix<typename Ring::Element> dummyMultiline;

		commentator.start("[Indexer] constructing sub matrices 2");
			idxr2.constructSubMatrices(A, sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime, free_memory_on_the_go);
		commentator.stop(MSG_DONE);
		MatrixUtils::show_mem_usage("[constructing sub matrices 2]"); report << endl;

		commentator.start("[Bloc] B1 = A1^-1 B1");
			Level3Ops::reducePivotsByPivots(R, sub_A_prime, sub_B_prime);
		commentator.stop(MSG_DONE);
		sub_A_prime.free(true);
		MatrixUtils::show_mem_usage("[B1 = A1^-1 B1]"); report << endl;

		SequentialIndexer<typename Ring::Element, IndexType> inner_dummy_idxr;
		inner_dummy_idxr.processMatrix(dummySparse);
		idxr2.combineInnerIndexer(inner_dummy_idxr, true);

		commentator.start("[Bloc] Reconstructing final matrix");
			idxr2.reconstructMatrix(A, sub_B_prime, free_memory_on_the_go);
		commentator.stop(MSG_DONE); report << endl;
		MatrixUtils::show_mem_usage("[Reconstructing final matrix]"); report << endl;

		reduced = true;

commentator.stop("ROUND 2");
	}

commentator.stop("FGL BLOC NEW METHOD");

	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------------" << endl;
		report << "<< Computing reduced echelon form of the matrix using structured Gaussian elimination >>" << endl;

		commentator.start("Structured rref of original matrix", "STRUCTURED_RREF");
			size_t rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			//FAST
			for(uint32 i=0; i<A.rowdim ()/2; ++i)
			{
				A[i].swap(A[A.rowdim () - 1 - i]);
			}

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
bool testFaugereLachartre_new_method_multiline_C(const Ring& R,
						SparseMatrix<typename Ring::Element>& A, bool validate_results, bool only_D,
						bool free_memory_on_the_go, bool horizontal,
						bool reconstruct_old)
{
	Context<Ring> ctx (R);
	SequentialIndexer<typename Ring::Element, IndexType> outer_indexer;

	size_t rank;
	bool reduced = false;

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Index type " << __IndexType__ << endl;

	SparseMatrix<typename Ring::Element> M_orig;

	if(validate_results)
	{
		M_orig = SparseMatrix<typename Ring::Element> (A.rowdim(), A.coldim());
		BLAS3::copy(ctx, A, M_orig);
	}

commentator.start("FGL BLOC NEW METHOD");
commentator.start("ROUND 1");
	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > sub_A, sub_B, sub_C, sub_D;
	SparseMultilineMatrix<uint16> sub_D_multiline, sub_C_multiline, sub_A_multiline;

	commentator.start("[Bloc] construting submatrices");
		outer_indexer.constructSubMatrices(A, sub_A_multiline, sub_B, sub_C_multiline, sub_D, free_memory_on_the_go);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("[construting submatrices]"); report << endl;
	report << "Pivots found: " << outer_indexer.Npiv << endl << endl;

	//MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
 	//MatrixUtils::copy(sub_C_multiline, sub_C);

	SHOW_MATRIX_INFO_BLOC(sub_B);
	SHOW_MATRIX_INFO_BLOC(sub_D);
	SHOW_MATRIX_INFO_MULTILINE(sub_A_multiline);
	SHOW_MATRIX_INFO_MULTILINE(sub_C_multiline);
	report << endl;

	commentator.start("[Bloc] C = MatrixOps::reduceC", "[reduceC]");
		Level3Ops::reduceC(R, sub_A_multiline, sub_C_multiline);
	commentator.stop("[reduceC]");	report << endl;

	commentator.start("Copy sub_C_multiline to sub_C_bloc");
		Level3Ops::copyMultilineMatrixToBlocMatrixRTL(sub_C_multiline, sub_C, free_memory_on_the_go);
	commentator.stop("[Copy sub_C_multiline to sub_C_bloc]");	report << endl;
	MatrixUtils::show_mem_usage("[Copy sub_C_multiline to sub_C_bloc]"); report << endl;


	if(reconstruct_old)
	{
		commentator.start("Copy sub_A_multiline to sub_A_bloc");
			Level3Ops::copyMultilineMatrixToBlocMatrixRTL(sub_A_multiline, sub_A, free_memory_on_the_go);
		commentator.stop("[Copy sub_C_multiline to sub_C_bloc]");	report << endl;
		MatrixUtils::show_mem_usage("[Copy sub_A_multiline to sub_A_bloc]"); report << endl;
	}


	commentator.start("[Bloc] D = D - C*B", "[D = D - C*B]");
		Level3Ops::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D, false);
	commentator.stop("[D = D - C*B]");
	SHOW_MATRIX_INFO_BLOC(sub_D);
	sub_C.free(true);
	MatrixUtils::show_mem_usage("[D = D - C*B]"); report << endl;

	//	commentator.start("[Bloc] echelonize", "[echelonize D]");
//		rank = Level3Ops::echelonize(R, sub_D, sub_D_multiline, free_memory_on_the_go);
//	commentator.stop("[echelonize D]");

	commentator.start("[Bloc] echelonize", "[echelonize D]");
		rank = Level3Ops::echelonize(R, sub_D, sub_D_multiline, free_memory_on_the_go);
	commentator.stop("[echelonize D]");


//	if(parallel)
//		MatrixUtils::dumpMatrixAsPbmImage(sub_D_multiline, "sub_d_m_parallel.pbm");
//	else
//		MatrixUtils::dumpMatrixAsPbmImage(sub_D_multiline, "sub_d_m_seq.pbm");

	MatrixUtils::show_mem_usage("[echelonize D]"); report << endl;
	report << "Rank of D " << rank << endl;
	report << endl;


	SequentialIndexer<typename Ring::Element, IndexType> inner_idxr;
	commentator.start("[Bloc] Processing new matrix D");
		inner_idxr.processMatrix(sub_D_multiline);
	commentator.stop(MSG_DONE);
	report << "Pivots found: " << inner_idxr.Npiv << endl << endl;


	commentator.start("[Bloc] Combine inner indexer");
		outer_indexer.combineInnerIndexer(inner_idxr, true);
	commentator.stop(MSG_DONE); report << endl;

	//TODO: IF RREF, construct new matrices directly from sub_A, sub_B, sub_D_multiline
	///and skip this step
	commentator.start("[Bloc] Reconstructing matrix", "Reconstructing matrix]");
		if(reconstruct_old)
			outer_indexer.reconstructMatrix(A, sub_A, sub_B, sub_D_multiline, free_memory_on_the_go);
		else
			outer_indexer.reconstructMatrix(A, sub_A_multiline, sub_B, sub_D_multiline, free_memory_on_the_go);
	commentator.stop("[Reconstructing matrix]");
	MatrixUtils::show_mem_usage("[Reconstructing matrix]"); report << endl;



	report << "True rank " << outer_indexer.Npiv + rank << endl;
	outer_indexer.freeMemory (); //not used anymore
commentator.stop("ROUND 1");

	if(!only_D)
	{
report << "----------------------------------------------------------------------------------------" << endl;
commentator.start("ROUND 2");
		SequentialIndexer<typename Ring::Element, IndexType> idxr2;

		SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> >
				sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime;

		SparseMatrix<typename Ring::Element> dummySparse;
		SparseMultilineMatrix<typename Ring::Element> dummyMultiline;

		commentator.start("[Indexer] constructing sub matrices 2");
			idxr2.constructSubMatrices(A, sub_A_prime, sub_B_prime, sub_C_prime, sub_D_prime, free_memory_on_the_go);
		commentator.stop(MSG_DONE);
		MatrixUtils::show_mem_usage("[constructing sub matrices 2]"); report << endl;

		commentator.start("[Bloc] B1 = A1^-1 B1");
			Level3Ops::reducePivotsByPivots(R, sub_A_prime, sub_B_prime);
		commentator.stop(MSG_DONE);
		sub_A_prime.free(true);
		MatrixUtils::show_mem_usage("[B1 = A1^-1 B1]"); report << endl;

		SequentialIndexer<typename Ring::Element, IndexType>  inner_dummy_idxr;
		inner_dummy_idxr.processMatrix(dummySparse);
		idxr2.combineInnerIndexer(inner_dummy_idxr, true);

		commentator.start("[Bloc] Reconstructing final matrix");
			idxr2.reconstructMatrix(A, sub_B_prime, free_memory_on_the_go);
		commentator.stop(MSG_DONE); report << endl;
		MatrixUtils::show_mem_usage("[Reconstructing final matrix]"); report << endl;

		reduced = true;

		//MatrixUtils::dumpMatrixAsPbmImage(A, "A-rref.pbm");

commentator.stop("ROUND 2");
	}

commentator.stop("FGL BLOC NEW METHOD");

	bool pass =true;

	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------------" << endl;
		report << "<< Computing reduced echelon form of the matrix using structured Gaussian elimination >>" << endl;

		commentator.start("Structured rref of original matrix", "STRUCTURED_RREF");
			size_t rank_strucutured_gauss = StructuredGauss::echelonize_reduced(R, M_orig);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		if(!reduced)
		{
			//FAST
			for(uint32 i=0; i<A.rowdim ()/2; ++i)
			{
				A[i].swap(A[A.rowdim () - 1 - i]);
			}

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
	bool new_method = false;
	bool validate_results = false;
	bool free_mem = false;
	bool only_D = false;
	bool use_multiline = false;
	bool horizontal = false;
	bool reconstruct_old = false;

	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		{ 'n', "-n", "Use the new method (computations on CD only)", TYPE_NONE, &new_method },
		{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ 'k', "-k", "**DO NOT** free memory as early as possible", TYPE_NONE, &free_mem},
		{ 'd', "-d", "Compute only the Gauss elimination of the matrix (no RREF form)", TYPE_NONE, &only_D },
		{ 'm', "-m", "Computation on submatrix C as a multiline matrix and not a bloc matrix", TYPE_NONE, &use_multiline },
		{ 'u', "-u", "Perform parallel computations horizontally (row majot then column)", TYPE_NONE, &horizontal },
		{ 'c', "-c", "Use old reconstrt matrix (doesn't matter)", TYPE_NONE, &reconstruct_old },
		{ '\0' }
	};

	parseArguments(argc, argv, args, "", 0);
	free_mem = !free_mem;
	only_D = !only_D;

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDepth(3);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(PROGRESS_REPORT).setMaxDepth(3);

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Faugère-Lachartre Bloc Version", "Faugère-Lachartre Bloc Version");

	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtils::loadF4Modulus(fileName);
	Modular<modulus_type> R(modulus);

	Context<Ring> ctx (R);
	MatrixUtils::show_mem_usage("starting");


	commentator.start("Loading matrix loadF4Matrix__low_memory SYS CALL");
		SparseMatrix<Ring::Element> A;
		MatrixUtils::loadF4Matrix__low_memory_syscall_no_checks(R, A, fileName);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("Loading matrix");

	report << endl;

	if(use_multiline)
		pass = testFaugereLachartre_new_method_multiline_C(R, A, validate_results, only_D, free_mem, horizontal,
				reconstruct_old);

	else
	{
		if(new_method)
			pass = testFaugereLachartre_new_method(R, A, validate_results, only_D, free_mem);
		else
			pass = testFaugereLachartre_old_method(R, A, validate_results, only_D, free_mem, horizontal);
	}

	SHOW_MATRIX_INFO_SPARSE(A);
	MatrixUtils::show_mem_usage("[In main]"); report << endl;
	commentator.stop("Faugère-Lachartre Bloc Version");

	return pass ? 0 : -1;
}

