/*
 * test-faugere-lachartre.C
 * Copyright 2012 Martani Fayssal (LIP6 / UPMC University Paris06)
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 *
 */



#include <iostream>

#include "matrix-util.h"
#include "matrix-op.h"
#include "indexer.h"
#include "structured-gauss-lib.h"

#include "../util/support.h"

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

#include "lela/algorithms/faugere-lachartre.h"

using namespace LELA;
using namespace std;

double matrixSizeInBytes(uint64 nb_elements)
{
	return (double)(nb_elements * 4 + nb_elements * 2) / 1024.0 / 1024.0;
}

template <typename Matrix>
void showMatrixSizeAndDensity(Matrix& A, const char *matrix_name)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::pair<uint64, double> size_density;

	size_density = MatrixUtil::getMatrixSizeAndDensity(A);
	report << "Nb elements of matrix " << matrix_name << ": "
		    << size_density.first
		    << "\tdensity: " << size_density.second
		    << " size: " << matrixSizeInBytes(size_density.first) << " MB" << endl;
}

bool testFaugereLachartre(const char *file_name, bool validate_results = false)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;
	
	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R (modulus);
	Context<Ring> ctx (R);
	
	Indexer<uint32> indexer;
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	SparseMatrix<Ring::Element> C (A.rowdim(), A.coldim());
	
	if(validate_results)
		BLAS3::copy(ctx, A, C);

	report << endl;
commentator.start("FAUGERE LACHARTRE", "FAUGERE_LACHARTRE");
	
///First round
	commentator.start("[Indexer] constructing indexes");
		indexer.processMatrix(A);
	commentator.stop(MSG_DONE);
    MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << indexer.Npiv << endl << endl;

	SparseMatrix<Ring::Element>	sub_A (indexer.Npiv, indexer.Npiv),
					sub_B (indexer.Npiv, A.coldim () - indexer.Npiv),
					sub_C (A.rowdim () - indexer.Npiv, indexer.Npiv),
					sub_D (A.rowdim () - indexer.Npiv, A.coldim () - indexer.Npiv);

	commentator.start("[Indexer] constructing sub matrices");
		indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] constructing sub matrices");
	
	commentator.start("B = A^-1 x B");
		MatrixOp::reducePivotsByPivots(R, sub_A, sub_B);
		//Equivalent to BLAS3::trsm (ctx, ctx.F.one (), sub_A, sub_B, UpperTriangular, false);
	commentator.stop(MSG_DONE);

	//MatrixUtil::freeMatrixMemory(sub_A);
	MatrixUtil::show_mem_usage("B = A^-1 x B");

	commentator.start("D <- D - CB");
		MatrixOp::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
		//Equivalent to BLAS3::gemm
	commentator.stop(MSG_DONE);
	//MatrixUtil::freeMatrixMemory(sub_C);
	MatrixUtil::show_mem_usage("D <- D - CB");
	
///D
	report << endl;
	size_t rank;

	commentator.start("Echelonize (D)");
		rank = MatrixOp::echelonize(R, sub_D);
		MatrixUtil::makeRowsUnitary(R, sub_D);		//D pivots are not necessarily unitary
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Echelonize (D)");

	report << "Rank of D: " << rank << endl;

///Second round
	report << endl;
	Indexer<uint32> indexer2;

	commentator.start("[Indexer] processing matrix D");
		indexer2.processMatrix(sub_D);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] processing matrix D");

	SparseMatrix<Ring::Element>	D1 (indexer2.Npiv, indexer2.Npiv),
					D2 (indexer2.Npiv, sub_D.coldim () - indexer2.Npiv),

					B1 (sub_B.rowdim (), indexer2.Npiv),
					B2 (sub_B.rowdim (), sub_B.coldim() - indexer2.Npiv);

	commentator.start("[Indexer] constructing submatrices B1, B1, D1, D2");
		indexer2.constructSubMatrices(sub_B, sub_D, B1, B2, D1, D2, true);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] constructing submatrices B1, B1, D1, D2");

	commentator.start("D2 = D1^-1 x D2");
		MatrixOp::reducePivotsByPivots(R, D1, D2);
	commentator.stop(MSG_DONE);
	//MatrixUtil::show_mem_usage("D2 = D1^-1 x D2");
	//MatrixUtil::freeMatrixMemory(D1);
	
	//TODO: This operation is very inefficient, since it uses sparse-sparse operations on almost dense matrices
	commentator.start("B2 <- B2 - D2 D1");
		MatrixOp::reduceNonPivotsByPivots(R, B1, D2, B2);
	commentator.stop(MSG_DONE);
	//MatrixUtil::freeMatrixMemory(B1);
	//MatrixUtil::show_mem_usage("B2 <- B2 - D2 D1");

	report << endl;
	commentator.start("[Indexer] Reconstructing indexes");
		indexer.combineInnerIndexer(indexer2);
	commentator.stop(MSG_DONE);

	commentator.start("[Indexer] Reconstructing matrix");
		indexer.reconstructMatrix(A, B2, D2);
	commentator.stop(MSG_DONE);

commentator.stop("FAUGERE LACHARTRE");

	showMatrixSizeAndDensity(A, "A");
	report << endl;
	
	if(validate_results)
	{
		report << "----------------------------------------------------------------------------------" << endl;
		report << "Computing reduced echelon form of the matrix using structured Gaussian elimination" << endl;
		
		commentator.start("Structured rref", "STRUCTURED_RREF");
			StructuredGauss::echelonize_reduced(R, C);
		commentator.stop(MSG_DONE, "STRUCTURED_RREF");

		report << endl;
	
		if(BLAS3::equal(ctx, A, C))
		{
			report << "Result CORRECT" << std::endl;
			report << endl;
			return true;
		}
		else
		{
			report << "Result NOT OK" << std::endl;
			report << endl;
			return false;
		}
	}
	
	return true;
}

void testLelaFaugereLachartre(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;
	
	modulus_type modulus = (modulus_type)MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R (modulus);
	Context<Ring> ctx (R);
	
	FaugereLachartre<Ring> Solver (ctx);
	size_t rank;
	Ring::Element det;

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
		MatrixUtil::invertMatrixRows(A);	//In case of general matrices, sort the rows by increasing row entry
	commentator.stop(MSG_DONE);
	
	Solver.echelonize (A, A, rank, det);
}

int main (int argc, char **argv)
{
	const char *fileName = "";
	bool validate_results = false;
	bool useLelaFG_Lachartre = false;

	bool pass = true;

	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName},
		{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ 'l', "-l", "Use Lela's faugere-lachartre implementation", TYPE_NONE, &useLelaFG_Lachartre}, 
		{ '\0' }
	};

	parseArguments (argc, argv, args, "", 0);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	commentator.start ("Faugère-Lachartre test suite", "FaugereLachartre");

	if(useLelaFG_Lachartre)
	{
		testLelaFaugereLachartre(fileName);
		pass = true;
	}else
		pass = testFaugereLachartre(fileName, validate_results);

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}
