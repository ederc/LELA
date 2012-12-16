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

//#include "echelonize_multiline.C"

#define SHOW_MATRIX_INFO(M)			\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << "(" << right << setw(8) << M.rowdim() << "x " << left << setw(8) << M.coldim() << ")\t"	\
			<< "- density: " << MatrixUtils::getMatrixSizeAndDensity(M, true).second << " % \t bloc dim "	\
			<< M.bloc_height() << " x " << M.bloc_width() << "  ****" << endl;


using namespace LELA;
using namespace std;

template <typename Ring>
bool testFaugereLachartre_bloc(const Ring& R, SparseMatrix<typename Ring::Element>& A) //, bool only_D, bool parallel, bool validate_results)
{
	Context<Ring> ctx (R);
	Indexer_ idx;
	Indexer<uint32> simple_indexer;

#ifdef DEFAULT_BLOC_HEIGHT
#if DEFAULT_BLOC_HEIGHT>=256
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

	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > sub_A, sub_B, sub_C, sub_D;
	SparseMultilineMatrix<typename Ring::Element> sub_D_multiline;

	commentator.start("[Bloc] construting submatrices");
		idx.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);


	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);

	commentator.start("[Simple] construting submatrices");
		simple_indexer.processMatrix(A);

		SparseMatrix<typename Ring::Element> sub_A__(simple_indexer.Npiv, simple_indexer.Npiv),
										 sub_B__(simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv),
										 sub_C__(A.rowdim() - simple_indexer.Npiv, simple_indexer.Npiv),
										 sub_D__(A.rowdim() - simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv);

		simple_indexer.constructSubMatrices(A, sub_A__, sub_B__, sub_C__, sub_D__, false);
	commentator.stop(MSG_DONE);

	/*if (MatrixUtils::equal(R, sub_A, sub_A__))
		report << "<<<<<<< A SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> A FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_B, sub_B__))
		report << "<<<<<<< B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> B FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_C, sub_C__))
		report << "<<<<<<< C SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> C FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_D, sub_D__))
		report << "<<<<<<< D SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> D FAILURE <<<<<<<<<<" << endl;*/

//	MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_C__, "sub_C__.pbm");
//
//	return true;

//	uint32 nb = 0;
//	for(uint32 i=0; i<sub_D.rowBlocDim(); ++i)
//	{
//		//cout << "Row " << i << endl;
//		for(uint32 j=0; j<sub_D[i].size (); ++j)
//		{
//			//cout << "\tBloc " << j << endl;
//			for(uint32 k=0; k<sub_D[i][j].bloc_height(); ++k)
//			{
//				//cout << "\t\tline " << k << endl;
//				if(!(sub_D[i][j][k].is_sparse ()))
//					nb++;
//			}
//		}
//	}
//
//	cout << "nb hybrid vectors " << nb << endl;

	/*for(uint16 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		uint16 t=0;
		uint16 v1,v2, idx;
		for(uint16 j=0; j<sub_A[0][0][i].size(); ++j)
		{
			idx = sub_A[0][0][i].IndexData[j];
			v1 = sub_A[0][0][i].at(0, j);
			cout << "(" << idx << "," << v1 << ")\t-\t";
		}

		t = 0;
		cout << endl;

		for(uint16 j=0; j<sub_A[0][0][i].size(); ++j)
		{
			idx = sub_A[0][0][i].IndexData[j];
			v1 = sub_A[0][0][i].at(1, j);
			cout << "(" << idx << "," << v1 << ")\t-\t";
		}

		cout << endl;
	}

	cout <<" \n***************************\n";
	for(uint16 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		for(uint16 j=0; j<sub_A__[sub_A__.rowdim()-1-i].size(); ++j)
		{
			cout << "(" << sub_A__[sub_A__.rowdim()-1-i][j].first << ","
					<< sub_A__[sub_A__.rowdim()-1-i][j].second << ")\t-\t";
		}
		cout << endl;
	}*/


//	MatrixUtils::dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_B.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D.pbm");

	commentator.start("BLOC VERSION", "BLOC VERSION");
		commentator.start("[Bloc] B = A^-1 B");
			MatrixOps::reducePivotsByPivots(R, sub_A, sub_B);
		commentator.stop(MSG_DONE);

		commentator.start("[Bloc] D = D - C*B");
			MatrixOps::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
		commentator.stop(MSG_DONE);



//		commentator.start("echelonize D");
//			MatrixOps::echelonize(R, sub_D, sub_D_multiline);
//		commentator.stop(MSG_DONE);
	commentator.stop("BLOC VERSION");

	return true;

	commentator.start("SIMPLE VERSION", "SIMPLE VERSION");
		commentator.start("[Simple] B = A^-1 B");
			MatrixOp::reducePivotsByPivots(R, sub_A__, sub_B__);
		commentator.stop(MSG_DONE);

		commentator.start("[Simple] D = D - C*B");
			MatrixOp::reduceNonPivotsByPivots(R, sub_C__, sub_B__, sub_D__);
		commentator.stop(MSG_DONE);
	commentator.stop("SIMPLE VERSION");

	if (MatrixUtils::equal(R, sub_B, sub_B__))
		report << "<<<<<<< B = A^-1 B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> B = A^-1 B FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_D, sub_D__))
		report << "<<<<<<< D = D - C*B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> D = D - C*B FAILURE <<<<<<<<<<" << endl;


//	MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D_bloc.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_D__, "sub_D_simple.pbm");
//
//	SparseMatrix<typename Ring::Element> tmp (sub_D.rowdim(), sub_D.coldim());
//	MatrixUtils::copy(sub_D, tmp);
//	//MatrixUtils::dumpMatrixAsPbmImage(tmp, "tmp.pbm");
//
//	for(int i=sub_D.rowdim()-1; i>=0; --i)
//	{
//		if(!(sub_D__[i] == tmp[i]))
//		{
//			cout << "DIFFER ON " << i << endl;
//			//BLAS1::write(ctx, cout, sub_B__[i]);
//			//cout << endl << endl;
//			//BLAS1::write(ctx, cout, tmp[i]);
//		}
//	}


	//BLAS3::write(ctx, cout, tmp, FORMAT_SAGE);

	/*for(uint32 i=0; i<(uint32)std::ceil((double)sub_B.coldim() / sub_B.block_width()); ++i)
	{
		for(uint32 k=0; k<sub_B.block_height(); ++k)
		{
			BLAS1::write(ctx, cout, sub_B[0][i][k]);
			cout << endl;
		}
	}*/

	return true;
}


template <typename Ring>
bool testFaugereLachartre_bloc_CD(const Ring& R, SparseMatrix<typename Ring::Element>& A) //, bool only_D, bool parallel, bool validate_results)
{
	Context<Ring> ctx (R);
	Indexer_ idx;
	Indexer<uint32> simple_indexer;

#ifdef DEFAULT_BLOC_HEIGHT
#if DEFAULT_BLOC_HEIGHT>=256
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

	SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, IndexType> > sub_A, sub_B, sub_C, sub_D;
	SparseMultilineMatrix<uint16> sub_D_multiline;
	
	commentator.start("[Bloc] construting submatrices");
		idx.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);

	SHOW_MATRIX_INFO(sub_A);
	SHOW_MATRIX_INFO(sub_B);
	SHOW_MATRIX_INFO(sub_C);
	SHOW_MATRIX_INFO(sub_D);

//	commentator.start("[Simple] construting submatrices");
//		simple_indexer.processMatrix(A);
//
//		SparseMatrix<typename Ring::Element> sub_A__(simple_indexer.Npiv, simple_indexer.Npiv),
//										 sub_B__(simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv),
//										 sub_C__(A.rowdim() - simple_indexer.Npiv, simple_indexer.Npiv),
//										 sub_D__(A.rowdim() - simple_indexer.Npiv, A.coldim() - simple_indexer.Npiv);
//
//		simple_indexer.constructSubMatrices(A, sub_A__, sub_B__, sub_C__, sub_D__, false);
//	commentator.stop(MSG_DONE);

	/*if (MatrixUtils::equal(R, sub_A, sub_A__))
		report << "<<<<<<< A SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> A FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_B, sub_B__))
		report << "<<<<<<< B SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> B FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_C, sub_C__))
		report << "<<<<<<< C SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> C FAILURE <<<<<<<<<<" << endl;

	if (MatrixUtils::equal(R, sub_D, sub_D__))
		report << "<<<<<<< D SUCESS >>>>>>>>>>" << endl;
	else
		report << ">>>>>>> D FAILURE <<<<<<<<<<" << endl;*/


//	MatrixUtils::dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_B, "sub_B.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D.pbm");



	commentator.start("BLOC VERSION", "BLOC VERSION");
		commentator.start("[Bloc] C = MatrixOps::reduceC");
			MatrixOps::reduceC(R, sub_A, sub_C);
		commentator.stop(MSG_DONE);

		commentator.start("[Bloc] D = D - C*B");
			MatrixOps::reduceNonPivotsByPivots(R, sub_C, sub_B, sub_D);
		commentator.stop(MSG_DONE);

		//MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D1.pbm");

		commentator.start("[Bloc] echelonize_blocs_column");
			//ECHEL::echelonize_blocs_column(R, sub_D, 0);
			MatrixOps::echelonize(R, sub_D, sub_D_multiline, false);
		commentator.stop(MSG_DONE);

	commentator.stop("BLOC VERSION");

//	MatrixUtils::dumpMatrixAsPbmImage(sub_C, "sub_C.pbm");
//	MatrixUtils::dumpMatrixAsPbmImage(sub_D, "sub_D2.pbm");

	return true;
}


template <class Ring, class SparseMat>
void testMatrixUtils(const Ring& R, const SparseMat& A)
{
	Context<Ring> ctx (R);
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	typedef SparseBlocMatrix<SparseMultilineBloc<typename Ring::Element, uint16> > Matrix;

	commentator.start("TESTING ArrangementTopDown_LeftRight", "ArrangementTopDown_LeftRight");
	{
		Matrix M0(A.rowdim(), A.coldim(), Matrix::ArrangementTopDown_LeftRight);
		SparseMatrix<typename Ring::Element> C(A.rowdim(), A.coldim());

		report << " BLOC HEIGHT" << M0.bloc_height () << endl;

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(A, M0);
		commentator.stop(MSG_DONE);

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(M0, C);
		commentator.stop(MSG_DONE);

	//MatrixUtil::dumpMatrixAsPbmImage(A, "A.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(C, "C.pbm");

	report << endl;
	if(BLAS3::equal(ctx, A, C))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << ">>>>>>>>>>>>>>>>>>>>>>>>>>>RESULT WRONG<<<<<<<<<<<<<<<<<<<<<";
	report << endl;
	MatrixUtils::show_mem_usage("ArrangementTopDown_LeftRight");
	}
	commentator.stop("ArrangementTopDown_LeftRight");

	commentator.start("TESTING ArrangementDownTop_LeftRight", "ArrangementDownTop_LeftRight");
	{
		Matrix M0(A.rowdim(), A.coldim(), Matrix::ArrangementDownTop_LeftRight);
		SparseMatrix<typename Ring::Element> C(A.rowdim(), A.coldim());

		report << " BLOC HEIGHT" << M0.bloc_height () << endl;

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(A, M0);
		commentator.stop(MSG_DONE);

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(M0, C);
		commentator.stop(MSG_DONE);

	//MatrixUtil::dumpMatrixAsPbmImage(A, "A.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(C, "C.pbm");

	report << endl;
	if(BLAS3::equal(ctx, A, C))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << ">>>>>>>>>>>>>>>>>>>>>>>>>>>RESULT WRONG<<<<<<<<<<<<<<<<<<<<<";
	report << endl;
	MatrixUtils::show_mem_usage("ArrangementDownTop_LeftRight");
	}
	commentator.stop("ArrangementDownTop_LeftRight");

	commentator.start("TESTING ArrangementDownTop_RightLeft", "ArrangementDownTop_RightLeft");
	{
		Matrix M0(A.rowdim(), A.coldim(), Matrix::ArrangementDownTop_RightLeft);
		SparseMatrix<typename Ring::Element> C(A.rowdim(), A.coldim());

		report << " BLOC HEIGHT" << M0.bloc_height () << endl;

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(A, M0);
		commentator.stop(MSG_DONE);

		commentator.start("Copy SparseMatrix => SparseBlocMatrix");
			MatrixUtils::copy(M0, C);
		commentator.stop(MSG_DONE);

	//MatrixUtil::dumpMatrixAsPbmImage(A, "A.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(C, "C.pbm");

	report << endl;
	if(BLAS3::equal(ctx, A, C))
		report << "<<<<<<<<<<<<<<<<<<<<<<<<<<<RESULT OK>>>>>>>>>>>>>>>>>>>>>>>>";
	else
		report << ">>>>>>>>>>>>>>>>>>>>>>>>>>>RESULT WRONG<<<<<<<<<<<<<<<<<<<<<";
	report << endl;
	MatrixUtils::show_mem_usage("ArrangementDownTop_RightLeft");
	}
	commentator.stop("ArrangementDownTop_RightLeft");
}




int main(int argc, char **argv)
{
	char *fileName = NULL;

	bool pass = true;
	bool CD = false;

	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		//{ 'r', "-r", "Validate the results by comparing them to structured Gauss", TYPE_NONE, &validate_results },
		{ 'n', "-n", "Use the new method (computations on CD only)", TYPE_NONE, &CD },

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


	commentator.start("Loading matrix loadF4Matrix__low_memory SYS CALL");
		SparseMatrix<Ring::Element> A;
		MatrixUtils::loadF4Matrix__low_memory_syscall_no_checks(R, A, fileName);
	commentator.stop(MSG_DONE);
	MatrixUtils::show_mem_usage("Loading matrix");

	report << endl;

//	commentator.start("Loading matrix loadF4Matrix__low_memory LIBC");
//		SparseMatrix<Ring::Element> A;
//		MatrixUtils::loadF4Matrix__low_memory(R, A, fileName);
//	commentator.stop(MSG_DONE);
//	MatrixUtils::show_mem_usage("Loading matrix");
//	report << endl;
//
//
//	if(BLAS3::equal(ctx, A, AA))
//		report << "OK" << endl;
//	else
//		report << "NOT OK" << endl;

	//testMatrixUtils(R, A);
	if(CD)
		pass = testFaugereLachartre_bloc_CD(R, A);
	else
		pass = testFaugereLachartre_bloc(R, A);

	report << endl;

	commentator.stop(MSG_STATUS (pass));

	return pass ? 0 : -1;
}
