/*
 * test-FG-bloc.C
 *
 *  Created on: 4 juil. 2012
 *      Author: martani
 */


#include "types.h"
#include "matrix-utils.h"
#include "matrix-ops.h"
#include "indexer.h"
#include "../../../util/support.h"

#include "../only-D/indexer.h"
#include "../only-D/matrix-util.h"
#include "../only-D/matrix-op.h"

#define SHOW_MATRIX_INFO(M)			\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << "(" << right << setw(8) << M.rowdim() << "x " << left << setw(8) << M.coldim() << ")\t"	\
			<< "- density: " << MatrixUtil::getMatrixSizeAndDensity(M).second << "% \t bloc dim "	\
			<< M.block_height() << " x " << M.block_width() << "  ****" << endl;


using namespace LELA;
using namespace std;

template <typename Ring>
bool testFaugereLachartre_bloc(const Ring& R, SparseMatrix<typename Ring::Element>& A) //, bool only_D, bool parallel, bool validate_results)
{
	Context<Ring> ctx (R);
	Indexer_ idx;
	Indexer<uint32> simple_indexer;

#ifdef DEFAULT_BLOC_HEIGHT
#if DEFAULT_BLOC_HEIGHT>256
	typedef uint16 IndexType;
	string type = "uint16";
#else
	typedef uint8 IndexType;
	string type = "uint8";
#endif
#else
	typedef uint16 IndexType;
	string type = "uint16";
#endif

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Index type " << type << endl;

	SparseBlocMatrix<typename Ring::Element, IndexType> sub_A, sub_B, sub_C, sub_D;

	commentator.start("[Bloc] construting submatrices");
		idx.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);

	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);

	/*typename SparseBlocMatrix<typename Ring::Element, uint16>::RowIterator i_;
	int x=0;
	for(i_ = sub_A.rowBegin(); i_ != sub_A.rowEnd(); ++i_)
	{
		report << "ROW " << x << " nb blocs " << i_->size() << " start idx" << sub_A.FirstBlocsColumIndexes[x++] << endl;
		for(uint32 j=0; j<i_->size(); ++j)
		{
			if((*i_)[j].empty ())
				report << "empty\t";
			else
				report << "full\t";

		}

		report << endl;
	}*/

	//MatrixUtils::dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");

	commentator.start("[Simple] construting submatrices");
		simple_indexer.processMatrix(A);

		SparseMatrix<typename Ring::Element> sub_A__(simple_indexer.Npiv, simple_indexer.Npiv),
										 sub_B__(simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv),
										 sub_C__(A.rowdim() - simple_indexer.Npiv, simple_indexer.Npiv),
										 sub_D__(A.rowdim() - simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv);

		simple_indexer.constructSubMatrices(A, sub_A__, sub_B__, sub_C__, sub_D__, false);
	commentator.stop(MSG_DONE);
	//MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_b1.pbm");

	//MatrixUtils::dumpMatrixAsPbmImage(sub_A__, "sub_A__.pbm");

	commentator.start("[Bloc] B = A^-1 B");
		MatrixOps::reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop(MSG_DONE);

	commentator.start("[Bloc] D = D - C*B");
		MatrixOps::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
	commentator.stop(MSG_DONE);

	commentator.start("[Simple] B = A^-1 B");
		MatrixOp::reducePivotsByPivots(R, sub_A__, sub_B__);
	commentator.stop(MSG_DONE);

	commentator.start("[Simple] D = D - C*B");
		MatrixOp::reduceNonPivotsByPivots(R, sub_C__, sub_B__, sub_D__);
	commentator.stop(MSG_DONE);

	if (MatrixUtils::equal(R, sub_B, sub_B__))
		report << "<<<<<<< B = A^-1 B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> B = A^-1 B FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_D, sub_D__))
		report << "<<<<<<< D = D - C*B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> D = D - C*B FAILURE <<<<<<<<<<" << endl;


	//MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_b_bloc.pbm");
	//MatrixUtils::dumpMatrixAsPbmImage(sub_B__, "sub_b_simple.pbm");

	/*Context<Ring> ctx (R);

	SparseMatrix<typename Ring::Element> tmp (sub_B.rowdim(), sub_B.coldim());
	MatrixUtils::copy(sub_B, tmp);
	MatrixUtils::dumpMatrixAsPbmImage(tmp, "tmp.pbm");*/


	//BLAS3::write(ctx, cout, tmp, FORMAT_SAGE);

	/*for(uint32 i=0; i<(uint32)std::ceil((double)sub_B.coldim() / sub_B.block_width()); ++i)
	{
		for(uint32 k=0; k<sub_B.block_height(); ++k)
		{
			BLAS1::write(ctx, cout, sub_B[0][i][k]);
			cout << endl;
		}
	}*/

	/********************************************************************************************/


	/********************************************************************************************/


	return true;
}







int main(int argc, char **argv)
{
	char *fileName = NULL;

	bool pass = true;


	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		//{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ '\0' }
	};

	parseArguments(argc, argv, args, "", 0);

	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(5);
	commentator.getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDepth(3);
	commentator.getMessageClass(TIMING_MEASURE).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator.getMessageClass(PROGRESS_REPORT).setMaxDepth(3);

	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Faugère-Lachartre Bloc Version");

	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtils::loadF4Modulus(fileName);
	Modular<modulus_type> R(modulus);

	Context<Ring> ctx (R);
	MatrixUtils::show_mem_usage("starting");

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A;
		MatrixUtils::loadF4Matrix__low_memory(R, A, fileName);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("Loading matrix");
	report << endl;

	pass = testFaugereLachartre_bloc(R, A);

	//SparseBlocMatrix<Ring::Element> M0(A.rowdim(), A.coldim(), SparseBlocMatrix<Ring::Element>::ArrangementDownTop_RightLeft);
	//SparseBlocMatrix<Ring::Element> M1(A.rowdim(), A.coldim(), SparseBlocMatrix<Ring::Element>::ArrangementDownTop_LeftRight);
	/*SparseBlocMatrix<Ring::Element> M2(A.rowdim(), A.coldim(), SparseBlocMatrix<Ring::Element>::ArrangementDownTop_RightLeft, true);

	SparseBlocMatrix<Ring::Element> sub_A, sub_B, sub_C, sub_D;

	//cin >> fileName;

	commentator.start("Copy SparseMatrix => SparseBlocMatrix");
		//MatrixUtils::copy(A, M0);
		//MatrixUtils::copy(A, M1);
		MatrixUtils::copy(A, M2);
	commentator.stop(MSG_DONE);

	SparseMatrix<Ring::Element> C(A.rowdim(), A.coldim());
	commentator.start("Copy SparseMatrix => SparseBlocMatrix");
		MatrixUtils::copy(R, M, C);
	commentator.stop(MSG_DONE);

	//MatrixUtils::dumpMatrixAsPbmImage(M, "M.pbm");
	//MatrixUtils::dumpMatrixAsPbmImage(A, "A.pbm");

	//MatrixUtil::dumpMatrixAsPbmImage(A, "A.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(C, "C.pbm");*/

	/*if(MatrixUtils::equal(R, A, M0))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << "RESULT WRONG";

	report << endl;

	if(MatrixUtils::equal(R, A, M1))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << "RESULT WRONG";*/

	report << endl;

	/*if(MatrixUtils::equal(R, A, M2))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << "RESULT WRONG";*/

	report << endl;

	commentator.stop(MSG_STATUS (pass));

	return pass ? 0 : -1;
}
