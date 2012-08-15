/*
 * level3-ops.C
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */


#ifndef LEVEL3_OPS_C_
#define LEVEL3_OPS_C_

#include "level3-ops.h"

#include "types.h"
#include "level1-ops.h"
#include "level2-ops.h"
#include "consts-macros.h"



template<typename Index>
void Level3Ops::reducePivotsByPivots(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Bloc version" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) B.coldim() / B.bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) A.rowdim() / A.bloc_height());

	//for all columns of blocs of B
	for (uint32 i = 0; i < nb_column_blocs_B; ++i)
	{
		//report << "Column B " << i << std::endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 j = 0;
				j < nb_row_blocs_A;
				++j)
		{
			const uint32 first_bloc_idx = A.FirstBlocsColumIndexes[j] / A.bloc_width();
			const uint32 last_bloc_idx = MIN(A[j].size () - 1, j);

			//report << "\tRow A " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << std::endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j
			TIMER_START_(memsetBlocToZero);
				Level1Ops::memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				Level1Ops::copySparseBlocToDenseBlocArray(R, B[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "/" << nb_column_blocs_B << "\trow\t" << j << std::ends;
#endif
			//for all the blocs in the current row of A
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in A and B " << k + first_bloc_idx << std::endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					Level2Ops::reduceBlocByRectangularBloc(R, A[j][k], B[k + first_bloc_idx][i], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(reduceBlocByTriangularBloc);
				Level2Ops::reduceBlocByTriangularBloc(R, A[j][last_bloc_idx], dense_bloc);
			TIMER_STOP_(reduceBlocByTriangularBloc);

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i], false);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);


		}
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);
	TIMER_REPORT_(reduceBlocByRectangularBloc);
	TIMER_REPORT_(reduceBlocByTriangularBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}


template<typename Index>
void Level3Ops::reduceNonPivotsByPivots(const Modular<uint16>& R,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
		SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
		bool invert_scalars)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Bloc version" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(D.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(C.rowdim(), D.rowdim());
	check_equal_or_raise_exception(C.coldim(), B.rowdim());
	check_equal_or_raise_exception(B.coldim(), D.coldim());



	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, D.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_column_blocs_D = (uint32) std::ceil((double) D.coldim() / D.bloc_width());
	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)C.rowdim() / C.bloc_height());

	//for all columns of blocs of D
	for (uint32 i = 0; i < nb_column_blocs_D; ++i)
	{
		//report << "Colum D|B\t" << i << endl;
		//for all rows of blocs in C
		for (uint32 j = 0; j < nb_row_blocs_C; ++j)
		{
			const uint32 first_bloc_idx = C.FirstBlocsColumIndexes[j] / C.bloc_width();
			const uint32 last_bloc_idx = C[j].size ();

			//report << "\tRow C\t" << j << "\tfirst bloc: " << first_bloc_idx << " - last: " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			TIMER_START_(memsetBlocToZero);
				Level1Ops::memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				Level1Ops::copySparseBlocToDenseBlocArray(R, D[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "/" << nb_column_blocs_D << "\trow\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B\t" << k + first_bloc_idx << endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					Level2Ops::reduceBlocByRectangularBloc(R, C[j][k], B[k + first_bloc_idx][i], dense_bloc, invert_scalars);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, D[j][i]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);
		}
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);
	TIMER_REPORT_(reduceBlocByRectangularBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}

template<typename Index>
void Level3Ops::reduceC(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[MatrixOps::reduceC] In spec Modular<uint16> Bloc version" << std::endl;

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, C.bloc_width() * sizeof(uint64));

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;
	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)C.rowdim() / C.bloc_height ());
	const uint32 nb_col_blocs_C = (uint32)std::ceil((double)C.coldim() / C.bloc_width ());

	//for all columns of blocs of D
	for (uint32 i = 0; i < nb_row_blocs_C; ++i)
	{
		check_equal_or_raise_exception(C.FirstBlocsColumIndexes[i], 0);

		//for all rows of blocs in C
		for (int j = nb_col_blocs_C-1; j >= 0; --j)
		{
			//report << "\tColumn C|A\t" << j << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			TIMER_START_(memsetBlocToZero);
		 	    Level1Ops::memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				Level1Ops::copySparseBlocToDenseBlocArray(R, C[i][j], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\trow\t" << i << "/" << nb_row_blocs_C << "\tcolumn\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (int k = nb_col_blocs_C-1; k>j; --k)
			{
				TIMER_START_(reduceBlocByRectangularBloc);
					Level2Ops::reduceBlocByRectangularBloc_C(R, C[i][k], A[k][j], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(reduceBlocByTriangularBloc);
				Level2Ops::reduceBlocByTriangularBloc_C(R, A[j][j], dense_bloc);
			TIMER_STOP_(reduceBlocByTriangularBloc);

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, C[i][j]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);

		}
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(reduceBlocByRectangularBloc);
	TIMER_REPORT_(reduceBlocByTriangularBloc);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);
	
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}


//Possibly yields false results, don't use
template<typename Index>
void Level3Ops::reduceC(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseMultilineMatrix<uint16>& C)
{
	typedef Modular<uint16> Ring;
	
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3Ops::reduceC] In spec Modular<uint16> Bloc version" << std::endl;

	uint32 C_coldim = C.coldim ();
	
	uint64 *tmpDenseArray1C;
		posix_memalign((void**)&tmpDenseArray1C, 16, C_coldim * sizeof(uint64));
		
	uint64 *tmpDenseArray2C;
		posix_memalign((void**)&tmpDenseArray2C, 16, C_coldim * sizeof(uint64));
	
	for(uint32 i=0; i < C.rowdim()/NB_ROWS_PER_MULTILINE; ++i)
	{
 		//report << "row " << i;
		
		if(C[i].empty ())
			continue;
		
		Level1Ops::memsetToZero(&tmpDenseArray1C, 1, C_coldim);
		Level1Ops::memsetToZero(&tmpDenseArray2C, 1, C_coldim);
		
		
// 		report << " size " << C[i].size() << std::endl;
		Level1Ops::copyMultiLineVectorToDenseArray(C[i], tmpDenseArray1C, tmpDenseArray2C, C_coldim);
		
		typename Ring::Element Cv1_col1=0, Cv2_col1=0;
		uint32 Cp1=0;
		typename Ring::Element Cv1_col2=0, Cv2_col2=0;
		uint32 Cp2=0;

		int row_idx_in_A, bloc_idx_in_A, row_idx_in_bloc_A;
		
		uint32 idx=0;
		uint32 v__;
		
		uint32 start_idx;
		if(C[i].is_sparse (C_coldim)) 
			start_idx = C[i].IndexData[0];
		else
			start_idx = 0;		//should not happen!
		
		for(uint32 j=start_idx; j<C_coldim; ++j)
		{
			//report << "elt " << j << std::endl;
			
			Cp1 = j;
			Cv1_col1 = tmpDenseArray1C[j] % R._modulus;
			Cv2_col1 = tmpDenseArray2C[j] % R._modulus;

			if(Cv1_col1 == 0 && Cv2_col1 == 0)
				continue;

			if (Cv1_col1 != 0)
				R.negin(Cv1_col1);
			if (Cv2_col1 != 0)
				R.negin(Cv2_col1);
			
			row_idx_in_A = C_coldim - 1 - Cp1;
			bloc_idx_in_A = row_idx_in_A / A.bloc_height ();
			row_idx_in_bloc_A = (row_idx_in_A % A.bloc_height ()) / NB_ROWS_PER_MULTILINE;
			
			if ((C_coldim - 1 - Cp1) % 2 == 1)
			{
				Cp2 = j+1;
				R.copy(Cv1_col2, tmpDenseArray1C[j+1] % R._modulus);
				R.copy(Cv2_col2, tmpDenseArray2C[j+1] % R._modulus);

				uint32 sz = A[bloc_idx_in_A][bloc_idx_in_A][row_idx_in_bloc_A].size ();
				const uint32 val1 = A[bloc_idx_in_A][bloc_idx_in_A][row_idx_in_bloc_A].at_unchecked(1, sz-2);

				uint32 tmp;
				tmp = Cv1_col2 + (uint32)Cv1_col1*val1;
				R.init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1*val1;
				R.init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					R.negin(Cv1_col2);
				if (Cv2_col2 != 0)
					R.negin(Cv2_col2);


				//for all bloc columns in the corresponding bloc
				for(int k = bloc_idx_in_A; k>=0; --k)
				{
					const MultiLineVector<uint16, Index> *rowA = &(A[bloc_idx_in_A][k][row_idx_in_bloc_A]);

					if(rowA->empty ())
						continue;

					uint32 base_bloc_idx = C_coldim - 1 - k * A.bloc_width ();

					const uint32 N = rowA->size();
					const Index *p_idx = (N != 0) ? rowA->IndexData.getStartingPointer () : NULL;
					const uint16 *p_val = (N != 0) ? rowA->ValuesData.getStartingPointer () : NULL;
					const uint16 *p_val2 = p_val + 1;

					register uint32 v1__, v2__;
					register uint32 idx;

					for(uint32 l=0; l < N; l++)
					{
						idx = p_idx[l];
						v1__ = p_val[l*2];
						v2__ = p_val2[l*2];

						tmpDenseArray1C[base_bloc_idx - idx] += Cv1_col1 * v1__;
						tmpDenseArray1C[base_bloc_idx - idx] += Cv1_col2 * v2__;

						tmpDenseArray2C[base_bloc_idx - idx] += Cv2_col1 * v1__;
						tmpDenseArray2C[base_bloc_idx - idx] += Cv2_col2 * v2__;
					}
				}

				tmpDenseArray1C[Cp1] = Cv1_col1;
				tmpDenseArray2C[Cp1] = Cv2_col1;

				tmpDenseArray1C[Cp2] = Cv1_col2;
				tmpDenseArray2C[Cp2] = Cv2_col2;

				j++;
			}
			else
			{
				//for all bloc columns in the corresponding bloc
				for(int k = bloc_idx_in_A; k>=0; --k)
				{
					//report << "bloc " << k << std::endl;
					const MultiLineVector<uint16, Index> *rowA = &(A[bloc_idx_in_A][k][row_idx_in_bloc_A]);

					if(rowA->empty ())
						continue;

					const uint32 N = rowA->size();
					const Index *p_idx = (N != 0) ? rowA->IndexData.getStartingPointer () : NULL;
					const uint16 *p_val = (N != 0) ? rowA->ValuesData.getStartingPointer () : NULL;
					p_val += row_idx_in_bloc_A % NB_ROWS_PER_MULTILINE;

					uint32 base_bloc_idx = C_coldim - 1 - k * A.bloc_width ();

					//report << "base_bloc_idx " << base_bloc_idx << endl;

					for(uint16 l=0; l < N; l++)
					{
						idx = p_idx[l];
						v__ = p_val[l*2];

						if(base_bloc_idx < idx)
							report << "ERROR idx" << idx << " base idx " << base_bloc_idx << std::endl;

						tmpDenseArray1C[base_bloc_idx - idx] += v__ * Cv1_col1 ;
						tmpDenseArray2C[base_bloc_idx - idx] += v__ * Cv2_col1;
					}
					
	// 				Level2Ops::SparseScalMulSub__one_row__vect_array(
	// 					Cv1_col1,
	// 					Cv2_col1,
	// 					A[bloc_idx_in_A][k][row_idx_in_A % A.bloc_height ()],
	// 					row_idx_in_A % NB_ROWS_PER_MULTILINE,
	// 					tmpDenseArray1C,
	// 					tmpDenseArray2C);
				}
				
				tmpDenseArray1C[Cp1] = Cv1_col1;
				tmpDenseArray2C[Cp1] = Cv2_col1;
			}
		}
		
		Level1Ops::copyDenseArraysToMultilineVector(R, tmpDenseArray1C, tmpDenseArray2C, C_coldim, C[i]);

	}
}


//performs a pseudo reduction of C by A (A^-1 transpose(C))
template<typename Ring>
void Level3Ops::reduceC(const Ring& R, const SparseMultilineMatrix<uint16>& A, SparseMultilineMatrix<uint16>& C)
{
	typedef SparseMultilineMatrix<uint16> Matrix;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3Ops::reduceC - Sequential] In spec Modular<uint16> Bloc version" << std::endl;

	uint32 C_coldim = C.coldim ();
	uint32 C_rowdim_multiline = C.multiline_rowdim ();

	uint64 *tmpDenseArray1C;
		posix_memalign((void**)&tmpDenseArray1C, 16, C_coldim * sizeof(uint64));

	uint64 *tmpDenseArray2C;
		posix_memalign((void**)&tmpDenseArray2C, 16, C_coldim * sizeof(uint64));


	for(uint32 i=0; i < C_rowdim_multiline; ++i)
	{
#ifdef SHOW_PROGRESS
		report << "                                                                    \r";
		report << "\trow\t" << i << " / " << C_rowdim_multiline << std::ends;
#endif

		Level1Ops::memsetToZero(&tmpDenseArray1C, 1, C_coldim);
		Level1Ops::memsetToZero(&tmpDenseArray2C, 1, C_coldim);

		Level1Ops::copyMultiLineVectorToDenseArray(C[i], tmpDenseArray1C, tmpDenseArray2C, C_coldim);

		uint32 start_idx;
		if(C[i].is_sparse (C_coldim))
			start_idx = C[i].IndexData[0];
		else
			start_idx = 0;		//should not happen!

		typename Ring::Element Cv1_col1=0, Cv2_col1=0;
		uint32 Cp1=0;
		typename Ring::Element Cv1_col2=0, Cv2_col2=0;
		uint32 Cp2=0;

		uint32 tmp=0;
		uint32 row_in_A;

		for(uint32 j=start_idx; j<C_coldim; ++j)
		{
			Cp1 = j;
			Cv1_col1 = tmpDenseArray1C[j] % R._modulus;
			Cv2_col1 = tmpDenseArray2C[j] % R._modulus;

			if(Cv1_col1 == 0 && Cv2_col1 == 0)
				continue;

			if (Cv1_col1 != 0)
				Cv1_col1 = R._modulus - Cv1_col1;
			if (Cv2_col1 != 0)
				Cv2_col1 = R._modulus - Cv2_col1;

			row_in_A = (A.rowdim() - 1 - Cp1) / NB_ROWS_PER_MULTILINE;

			if ((A.rowdim() - 1 - Cp1) % 2 == 1)
			{
				Cp2 = j+1;
				Cv1_col2 = tmpDenseArray1C[Cp2] % R._modulus;
				Cv2_col2 = tmpDenseArray2C[Cp2] % R._modulus;

				uint16 v__ = A[row_in_A].at_unchecked(1, 1);

				tmp = Cv1_col2 + (uint32)Cv1_col1 * v__;
				R.init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1 * v__;
				R.init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					Cv1_col2 = R._modulus - Cv1_col2;
				if (Cv2_col2 != 0)
					Cv2_col2 = R._modulus - Cv2_col2;

				Level2Ops::SparseScalMulSub__two_rows__vect_array(
						Cv1_col2,
						Cv2_col2,
						Cv1_col1,
						Cv2_col1,
						A[row_in_A],
						tmpDenseArray1C,
						tmpDenseArray2C);

				tmpDenseArray1C[Cp1] = Cv1_col1;
				tmpDenseArray2C[Cp1] = Cv2_col1;

				tmpDenseArray1C[Cp2] = Cv1_col2;
				tmpDenseArray2C[Cp2] = Cv2_col2;

				++j;
			}
			else
			{
					Level2Ops::SparseScalMulSub__one_row__vect_array(
							Cv1_col1,
							Cv2_col1,
							A[row_in_A],
							(A.rowdim() - 1 - Cp1) % NB_ROWS_PER_MULTILINE,
							tmpDenseArray1C,
							tmpDenseArray2C);

					tmpDenseArray1C[Cp1] = Cv1_col1;
					tmpDenseArray2C[Cp1] = Cv2_col1;
			}
		}

		Level1Ops::copyDenseArraysToMultilineVector(R, tmpDenseArray1C, tmpDenseArray2C, C_coldim, C[i]);

	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif

		free(tmpDenseArray1C);
		free(tmpDenseArray2C);
}


template<typename Index>
uint32 Level3Ops::echelonize(const Modular<uint16>& R,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
			SparseMultilineMatrix<uint16>& outMatrix, bool destruct_in_matrix, bool use_hybrid_method)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Bloc version" << std::endl;

	typedef uint16 Element;
	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;
	outMatrix = SparseMultilineMatrix<uint16> (inMatrix.rowdim (), inMatrix.coldim ());

	TIMER_DECLARE_(CopyTimer);
	TIMER_DECLARE_(EchelonizeTimer);

	//copy bloc matrix to multiline matrix
	check_equal_or_raise_exception(inMatrix.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(inMatrix.isFilledWithEmptyBlocs (), true);

	uint32 curr_row_base = 0;




	TIMER_START_(CopyTimer);
	
	commentator.start("Copying Matrix");
	//this will copy the bloc matrix to a multiline matrix by inversing the order of the rows
	for(uint32 i=0; i<inMatrix.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
		if(inMatrix[i].empty ())
		{
			curr_row_base += inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;
			continue;
		}

		for (uint32 j = 0; j < inMatrix[i].size(); ++j) //for each bloc in the row
		{
			if(inMatrix[i][j].empty ())
				continue;

			//check_equal_or_raise_exception(inMatrix.FirstBlocsColumIndexes[i], 0);
			uint32 bloc_idx = inMatrix.FirstBlocsColumIndexes[i] + (inMatrix.bloc_width () * j);

			uint32 idx;
			Element val1, val2;
			for (uint16 k = 0;
					k < inMatrix[i][j].bloc_height(); ++k) //for each row in the bloc
			{
				if(inMatrix[i][j][k].empty ())
					continue;

				if(inMatrix[i][j][k].is_sparse ())
				{
					for (uint32 p = 0; p < inMatrix[i][j][k].size (); ++p)
					{
						idx = inMatrix[i][j][k].IndexData[p];
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + idx);
						outMatrix[curr_row_base+k].ValuesData.push_back(val1);
						outMatrix[curr_row_base+k].ValuesData.push_back(val2);
					}
				}
				else
				{
					for (uint32 p = 0; p < inMatrix.bloc_width (); ++p)
					{
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						if(val1 != 0 || val2 != 0)
						{
							outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + p);
							outMatrix[curr_row_base+k].ValuesData.push_back(val1);
							outMatrix[curr_row_base+k].ValuesData.push_back(val2);
						}
					}
				}

				if(destruct_in_matrix)
					inMatrix[i][j][k].free ();
			}

			if(destruct_in_matrix)
				inMatrix[i][j].free ();
		}

		curr_row_base += inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;
	}

	if(destruct_in_matrix)
		inMatrix.free ();
	TIMER_STOP_(CopyTimer);
	commentator.stop("Copying Matrix");
	
	uint32 ret;

	TIMER_START_(EchelonizeTimer);
	if(!use_hybrid_method)
		ret = echelonize_in_sparse(R, outMatrix);
	else
	{
		ret = echelonize_in_hybrid(R, outMatrix);
	}
	TIMER_STOP_(EchelonizeTimer);

	TIMER_REPORT_(CopyTimer);
	TIMER_REPORT_(EchelonizeTimer)

	return ret;
}


template <typename Index>
uint32 Level3Ops::echelonize_in_sparse(const Modular<uint16>& R, SparseMultilineMatrix<uint16, Index>& A)
{
	typedef Modular<uint16> Ring;
	typedef uint16 Element;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Level3Ops::echelonize_in  ____ sparse ____" << std::endl;

	uint32 coldim = A.coldim ();
	uint32 npiv = 0;
	uint32 npiv_real = 0;
	uint32 N = A.rowdim()/NB_ROWS_PER_MULTILINE + A.rowdim()%NB_ROWS_PER_MULTILINE;

	uint32 *piv;
	posix_memalign((void**)&piv, 16, N * sizeof(uint32));
	Level1Ops::memsetToZero(&piv, 1, N);

	MultiLineVector<Element, Index> *rowA;

	uint64 *tmpDenseArray1;
	posix_memalign((void**)&tmpDenseArray1, 16, coldim * sizeof(uint64));

	uint64 *tmpDenseArray2;
	posix_memalign((void**)&tmpDenseArray2, 16, coldim * sizeof(uint64));


	uint32 i=0;

	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;

	uint32 x=0;

	TIMER_DECLARE_(normalizeVectorTimer);
	TIMER_DECLARE_(normalizeArrayTimer);
	TIMER_DECLARE_(HeadVectorTimer);
	TIMER_DECLARE_(HeadArrayTimer);

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);

	TIMER_DECLARE_(SparseScalMulSub__two_rows__vect_arrayTimer);

	TIMER_START_(normalizeVectorTimer);
		Level1Ops::normalizeMultiLineVector(R, A[0]);
	TIMER_STOP_(normalizeVectorTimer);

	if(A.rowdim() == 0)
		return 0;

	for(i=0; i < N ; ++i)
	{
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << "/" << N << std::ends;
#endif

        TIMER_START_(RazArrayTimer);
        	Level1Ops::memsetToZero(&tmpDenseArray1, 1, coldim);
        	Level1Ops::memsetToZero(&tmpDenseArray2, 1, coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			Level1Ops::copyMultiLineVectorToDenseArray(A[i], tmpDenseArray1, tmpDenseArray2, coldim);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

        typename Ring::Element h_a1 = R.one (), h_a2 = R.one ();

		uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
		uint32 tmp=0;

		for(uint32 j=0; j<npiv; ++j)
		{
			rowA = &(A[piv[j]]);
			if(rowA->empty())
				continue;

			TIMER_START_(HeadVectorTimer);
				head_line1 = Level1Ops::headMultiLineVector(*rowA, 0, h_a1, head_line1_idx);
				head_line2 = Level1Ops::headMultiLineVector(*rowA, 1, h_a2, head_line2_idx);
			TIMER_STOP_(HeadVectorTimer);

#if DEBUG
			if(head_line1 != -1 && head_line1 == head_line2)
				throw std::logic_error ("Wrong Mutiline format");
			if(head_line1 > head_line2 && head_line2 != -1)			//makes the row with the smallest column entry first in the multiline
												//This should not arrive, see condition above copyDenseArrayToMultilineVector
			{
				uint32 t = head_line1;
				head_line1 = head_line2;
				head_line2 = t;

				t = head_line1_idx;
				head_line1_idx = head_line2_idx;
				head_line2_idx = t;

				//TODO::SLOW exchange data of the two lines
				uint16 tmp1=0;

				for(uint32 x=0; x<rowA->size (); ++x)
				{
					tmp1 = rowA->ValuesData[0+x*2];
					rowA->ValuesData[0+x*2] = rowA->ValuesData[1+x*2];
					rowA->ValuesData[1+x*2] = tmp1;
				}
			}
#endif

			if(head_line1 != -1)
			{
				//R.copy(v1col1, tmpDenseArray1[head_line1] % R._modulus);
				//R.copy(v2col1, tmpDenseArray2[head_line1] % R._modulus);
				v1col1 = tmpDenseArray1[head_line1] % R._modulus;
				v2col1 = tmpDenseArray2[head_line1] % R._modulus;

				if(v1col1 != 0)
					v1col1 = R._modulus - v1col1;
				if(v2col1 != 0)
					v2col1 = R._modulus - v2col1;
			}
			else
			{
				v1col1 = 0;
				v2col1 = 0;
			}

			if(head_line2 != -1)
			{
				//R.copy(v1col2, tmpDenseArray1[head_line2] % R._modulus);
				//R.copy(v2col2, tmpDenseArray2[head_line2] % R._modulus);
				v1col2 = tmpDenseArray1[head_line2] % R._modulus;
				v2col2 = tmpDenseArray2[head_line2] % R._modulus;

				uint16 val = rowA->at_unchecked(0, head_line2_idx);

				tmp = v1col2 + (uint32)v1col1*val;
				R.init(v1col2, tmp);
				tmp = v2col2 + (uint32)v2col1*val;
				R.init(v2col2, tmp);

				if (v1col2 != 0)
					v1col2 = R._modulus - v1col2;
				if (v2col2 != 0)
					v2col2 = R._modulus - v2col2;
			}
			else
			{
				v1col2 = 0;
				v2col2 = 0;
			}

			TIMER_START_(SparseScalMulSub__two_rows__vect_arrayTimer);
				Level2Ops::SparseScalMulSub__two_rows__vect_array(
						v1col1,
						v2col1,
						v1col2,
						v2col2,
						*rowA,
						tmpDenseArray1,
						tmpDenseArray2);
				TIMER_STOP_(SparseScalMulSub__two_rows__vect_arrayTimer);
		}

		TIMER_START_(normalizeArrayTimer);
			head_line1 = Level1Ops::normalizeDenseArray(R, tmpDenseArray1, coldim);
			head_line2 = Level1Ops::normalizeDenseArray(R, tmpDenseArray2, coldim);
		TIMER_STOP_(normalizeArrayTimer);

        //assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

        //reduce by same MultilineVector
        if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
        {
        	//assert(head_line1 >= 0 && head_line1 < coldim);
        	if(tmpDenseArray2[head_line1] % R._modulus != 0)
			{
        		uint32 h = tmpDenseArray2[head_line1] % R._modulus;
				//R.negin(h);
				h = R._modulus - h;

				register uint32 v__;
				for(x = head_line1; x<coldim; ++x)
				{
					v__ = tmpDenseArray1[x] & 0x000000000000ffff;
					tmpDenseArray2[x] += h * v__;
				}
			}
        }

		TIMER_START_(HeadArrayTimer);
			//head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
			head_line2 = Level1Ops::headDenseArray(R, tmpDenseArray2, coldim, h_a2);		//only this can change
        TIMER_STOP_(HeadArrayTimer);

        TIMER_START_(CopyDenseArrayToSparseVectorTimer);
        	if(head_line1 == -1 ||
        			((head_line2 != -1) && (head_line1 > head_line2)))			//saves the line with the smallest column entry first
        		Level1Ops::copyDenseArraysToMultilineVector(R, tmpDenseArray2, tmpDenseArray1, coldim, A[i]);
        	else
        		Level1Ops::copyDenseArraysToMultilineVector(R, tmpDenseArray1, tmpDenseArray2, coldim, A[i]);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		TIMER_START_(normalizeVectorTimer);
			Level1Ops::normalizeMultiLineVector(R, A[i]);
		TIMER_STOP_(normalizeVectorTimer);

		if(!A[i].empty())
		{
			piv[npiv] = i;
			npiv++;
		}

		if(head_line1 != -1)
			npiv_real++;
		if(head_line2 != -1)
			npiv_real++;

	}
#ifdef SHOW_PROGRESS
                report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(normalizeVectorTimer);
	TIMER_REPORT_(normalizeArrayTimer);
	TIMER_REPORT_(HeadVectorTimer);
	TIMER_REPORT_(HeadArrayTimer);

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);

	TIMER_REPORT_(SparseScalMulSub__two_rows__vect_arrayTimer);

	free(piv);
	free(tmpDenseArray1);
	free(tmpDenseArray2);

	return npiv_real;
}



template <typename Index>
uint32 Level3Ops::echelonize_in_hybrid(const Modular<uint16>& R, SparseMultilineMatrix<uint16, Index>& A)
{
	typedef Modular<uint16> Ring;
	typedef uint16 Element;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Level3Ops::echelonize_in ____ hybrid ____" << std::endl;

	uint32 coldim = A.coldim ();
	uint32 npiv = 0;
	uint32 npiv_real = 0;
	uint32 N = A.rowdim()/NB_ROWS_PER_MULTILINE + A.rowdim()%NB_ROWS_PER_MULTILINE;

	uint32 *piv;
	posix_memalign((void**)&piv, 16, N * sizeof(uint32));
	Level1Ops::memsetToZero(&piv, 1, N);

	MultiLineVector<Element, Index> *rowA;

	uint64 *tmpDenseArray1;
	posix_memalign((void**)&tmpDenseArray1, 16, coldim * sizeof(uint64));

	uint64 *tmpDenseArray2;
	posix_memalign((void**)&tmpDenseArray2, 16, coldim * sizeof(uint64));


	uint32 i=0;

	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;

	TIMER_DECLARE_(normalizeVectorTimer);
	TIMER_DECLARE_(normalizeArrayTimer);
	TIMER_DECLARE_(HeadVectorTimer);
	TIMER_DECLARE_(HeadArrayTimer);

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);

	TIMER_DECLARE_(SparseScalMulSub__two_rows__vect_arrayTimer);
	TIMER_DECLARE_(DenseScalMulSub__two_rows__vect_array__variable_sizeTimer);

	TIMER_START_(normalizeVectorTimer);
		Level1Ops::normalizeMultiLineVector(R, A[0]);
	TIMER_STOP_(normalizeVectorTimer);

	if(A.rowdim() == 0)
		return 0;

	for(i=0; i < N ; ++i)
	{
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << "/" << N << std::ends;
#endif

        TIMER_START_(RazArrayTimer);
        	Level1Ops::memsetToZero(&tmpDenseArray1, 1, coldim);
        	Level1Ops::memsetToZero(&tmpDenseArray2, 1, coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			Level1Ops::copyMultiLineVectorToDenseArray(A[i], tmpDenseArray1, tmpDenseArray2, coldim);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

        typename Ring::Element h_a1 = R.one (), h_a2 = R.one ();

		uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
		uint32 tmp=0;

		for(uint32 j=0; j<npiv; ++j)
		{
			rowA = &(A[piv[j]]);
			if(rowA->empty ())
			{
				continue;
			}

			TIMER_START_(HeadVectorTimer);
				head_line1 = Level1Ops::headMultiLineVectorHybrid(*rowA, 0, h_a1, head_line1_idx, coldim);
				head_line2 = Level1Ops::headMultiLineVectorHybrid(*rowA, 1, h_a2, head_line2_idx, coldim);
			TIMER_STOP_(HeadVectorTimer);

#if DEBUG
			if(head_line1 != -1 && head_line1 == head_line2)
			{
				std::cout << "IN DEBUG 1\n";
				throw std::logic_error ("Wrong Mutiline format");
			}
			if(head_line1 > head_line2 && head_line2 != -1)			//makes the row with the smallest column entry first in the multiline
												//This should not arrive, see condition above copyDenseArrayToMultilineVector
			{
				uint32 t = head_line1;
				head_line1 = head_line2;
				head_line2 = t;

				t = head_line1_idx;
				head_line1_idx = head_line2_idx;
				head_line2_idx = t;

				//TODO::SLOW exchange data of the two lines
				uint16 tmp1=0;

				for(uint32 x=0; x<rowA->size (); ++x)
				{
					tmp1 = rowA->ValuesData[0+x*2];
					rowA->ValuesData[0+x*2] = rowA->ValuesData[1+x*2];
					rowA->ValuesData[1+x*2] = tmp1;
				}
			}
#endif

			if(head_line1 != -1)
			{
				//R.copy(v1col1, tmpDenseArray1[head_line1] % R._modulus);
				//R.copy(v2col1, tmpDenseArray2[head_line1] % R._modulus);
				v1col1 = tmpDenseArray1[head_line1] % R._modulus;
				v2col1 = tmpDenseArray2[head_line1] % R._modulus;

				if(v1col1 != 0)
					v1col1 = R._modulus - v1col1;
				if(v2col1 != 0)
					v2col1 = R._modulus - v2col1;
			}
			else
			{
				v1col1 = 0;
				v2col1 = 0;
			}

			if(head_line2 != -1)
			{
				//R.copy(v1col2, tmpDenseArray1[head_line2] % R._modulus);
				//R.copy(v2col2, tmpDenseArray2[head_line2] % R._modulus);
				v1col2 = tmpDenseArray1[head_line2] % R._modulus;
				v2col2 = tmpDenseArray2[head_line2] % R._modulus;

				uint16 val = rowA->at_unchecked(0, head_line2_idx);

				tmp = v1col2 + (uint32)v1col1*val;
				R.init(v1col2, tmp);
				tmp = v2col2 + (uint32)v2col1*val;
				R.init(v2col2, tmp);

				if (v1col2 != 0)
					v1col2 = R._modulus - v1col2;
				if (v2col2 != 0)
					v2col2 = R._modulus - v2col2;
			}
			else
			{
				v1col2 = 0;
				v2col2 = 0;
			}

			if(rowA->is_sparse (coldim))
			{
				TIMER_START_(SparseScalMulSub__two_rows__vect_arrayTimer);
					Level2Ops::SparseScalMulSub__two_rows__vect_array(
							v1col1,
							v2col1,
							v1col2,
							v2col2,
							*rowA,
							tmpDenseArray1,
							tmpDenseArray2);
				TIMER_STOP_(SparseScalMulSub__two_rows__vect_arrayTimer);
			}
			else
			{
				TIMER_START_(DenseScalMulSub__two_rows__vect_array__variable_sizeTimer);
					Level2Ops::DenseScalMulSub__two_rows__vect_array__variable_size(
							v1col1,
							v2col1,
							v1col2,
							v2col2,
							*rowA,
							tmpDenseArray1,
							tmpDenseArray2,
							head_line1);
				TIMER_STOP_(DenseScalMulSub__two_rows__vect_array__variable_sizeTimer);
			}
		}

		TIMER_START_(normalizeArrayTimer);
			head_line1 = Level1Ops::normalizeDenseArray(R, tmpDenseArray1, coldim);
			head_line2 = Level1Ops::normalizeDenseArray(R, tmpDenseArray2, coldim);
		TIMER_STOP_(normalizeArrayTimer);

        //assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

        //reduce by same MultilineVector
        if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
        {
        	//assert(head_line1 >= 0 && head_line1 < coldim);
        	if(tmpDenseArray2[head_line1] % R._modulus != 0)
			{
				register uint32 h = tmpDenseArray2[head_line1] % R._modulus;
				//R.negin(h);
				h = R._modulus - h;

				register uint32 v__;
				for(uint32 x = head_line1; x<coldim; ++x)
				{
					v__ = tmpDenseArray1[x] & 0x000000000000ffff;
					tmpDenseArray2[x] += h * v__;
				}
			}
        }

        TIMER_START_(HeadArrayTimer);
			//head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
			head_line2 = Level1Ops::headDenseArray(R, tmpDenseArray2, coldim, h_a2);		//only this can change
        TIMER_STOP_(HeadArrayTimer);

        TIMER_START_(CopyDenseArrayToSparseVectorTimer);
        	if(head_line1 == -1 ||
        			((head_line2 != -1) && (head_line1 > head_line2)))			//saves the line with the smallest column entry first
				Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray2, tmpDenseArray1, coldim, A[i]);
			else
				Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray1, tmpDenseArray2, coldim, A[i]);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		TIMER_START_(normalizeVectorTimer);
			Level1Ops::normalizeMultiLineVector(R, A[i]);
		TIMER_STOP_(normalizeVectorTimer);

		if(!A[i].empty ())
		{
			piv[npiv] = i;
			npiv++;
		}

		if(head_line1 != -1)
			npiv_real++;
		if(head_line2 != -1)
			npiv_real++;

	}
#ifdef SHOW_PROGRESS
                report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(normalizeVectorTimer);
	TIMER_REPORT_(normalizeArrayTimer);
	TIMER_REPORT_(HeadVectorTimer);
	TIMER_REPORT_(HeadArrayTimer);

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);

	TIMER_REPORT_(SparseScalMulSub__two_rows__vect_arrayTimer);
	TIMER_REPORT_(DenseScalMulSub__two_rows__vect_array__variable_sizeTimer);

	free(piv);
	free(tmpDenseArray1);
	free(tmpDenseArray2);

	return npiv_real;
}


//Right to left copy of multiline to bloc matrix
//NOTE: initializes B
//handle only sparse multiline matrices
template<typename Element, typename Index>
void Level3Ops::copyMultilineMatrixToBlocMatrixRTL(SparseMultilineMatrix<Element>& A,
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B, bool destruct_riginal)
{
	typedef SparseBlocMatrix<SparseMultilineBloc<Element, Index> > Matrix;
	B =  Matrix (A.rowdim(), A.coldim(), Matrix::ArrangementDownTop_RightLeft);


	const uint32 nb_blocs_per_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());



	uint32 row_bloc_idx = 0;
	uint32 idx;
	Element e1, e2;
	Index elt_idx_in_line;
	uint32 bloc_idx_in_row, line_idx_in_bloc;

	for(uint32 i = 0; i<A.multiline_rowdim(); ++i)
	{
		row_bloc_idx = i / (B.bloc_height() / NB_ROWS_PER_MULTILINE);
		line_idx_in_bloc = i % (B.bloc_height() / NB_ROWS_PER_MULTILINE);
		//cout << "i " << i << " row bloc " << row_bloc_idx << endl;

		if(line_idx_in_bloc == 0)
		{
			B[row_bloc_idx].resize (nb_blocs_per_dim);

			for(uint32 k=0; k<nb_blocs_per_dim; ++k)
				B[row_bloc_idx][k].init(B.bloc_height (), B.bloc_width ());

			B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
		}

		for(int j=A[i].size ()-1; j>=0; --j)	//sparse
		{
			idx = A[i].IndexData[j];
			e1 = A[i].at_unchecked(0, j);
			e2 = A[i].at_unchecked(1, j);

			bloc_idx_in_row = (A.coldim () - 1 - idx) / B.bloc_width();
			elt_idx_in_line = (A.coldim () - 1 - idx) % B.bloc_width();

			B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back(elt_idx_in_line);
			B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back(e1);
			B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back(e2);
		}

		if(destruct_riginal)
			A[i].free ();
	}

	if(destruct_riginal)
		A.free ();

}



/***************************************************************************************************************/
template<typename Index>
void Level3Ops::reducePivotsByPivots_horizontal(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "<<Level3Ops::reducePivotsByPivots_horizontal>> Modular<uint16> Bloc version" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) B.coldim() / B.bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) A.rowdim() / A.bloc_height());

	//for all columns of blocs of B
	//for all rows of blocs in A (starting from the end)
	for (uint32 j = 0;
		j < nb_row_blocs_A;
		++j)
	{
		const uint32 first_bloc_idx = A.FirstBlocsColumIndexes[j] / A.bloc_width();
		const uint32 last_bloc_idx = min(A[j].size () - 1, j);
			
		for (uint32 i = 0; i < nb_column_blocs_B; ++i)
		{
		//report << "Column B " << i << std::endl;

			//report << "\tRow A " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << std::endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j
			TIMER_START_(memsetBlocToZero);
				Level1Ops::memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				Level1Ops::copySparseBlocToDenseBlocArray(R, B[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\trow\t" << j << "/" << nb_row_blocs_A << "\tcolumn\t" << i << "/" << nb_column_blocs_B << std::ends;
#endif
			//for all the blocs in the current row of A
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in A and B " << k + first_bloc_idx << std::endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					Level2Ops::reduceBlocByRectangularBloc(R, A[j][k], B[k + first_bloc_idx][i], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(reduceBlocByTriangularBloc);
				Level2Ops::reduceBlocByTriangularBloc(R, A[j][last_bloc_idx], dense_bloc);
			TIMER_STOP_(reduceBlocByTriangularBloc);

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i], false);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);


		}
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);
	TIMER_REPORT_(reduceBlocByRectangularBloc);
	TIMER_REPORT_(reduceBlocByTriangularBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}





template<typename Index>
void Level3Ops::reduceNonPivotsByPivots_horizontal(const Modular<uint16>& R,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
		SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
		bool invert_scalars)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "<<Level3Ops::reduceNonPivotsByPivots_horizontal>> Modular<uint16> Bloc version" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(D.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(C.rowdim(), D.rowdim());
	check_equal_or_raise_exception(C.coldim(), B.rowdim());
	check_equal_or_raise_exception(B.coldim(), D.coldim());



	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, D.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_column_blocs_D = (uint32) std::ceil((double) D.coldim() / D.bloc_width());
	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)C.rowdim() / C.bloc_height());

	for (uint32 j = 0; j < nb_row_blocs_C; ++j)
	{
		const uint32 first_bloc_idx = C.FirstBlocsColumIndexes[j] / C.bloc_width();
		const uint32 last_bloc_idx = C[j].size ();
		
		//for all columns of blocs of D
		for (uint32 i = 0; i < nb_column_blocs_D; ++i)
		{
		//report << "Colum D|B\t" << i << endl;
		//for all rows of blocs in C
		

			//report << "\tRow C\t" << j << "\tfirst bloc: " << first_bloc_idx << " - last: " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			TIMER_START_(memsetBlocToZero);
				Level1Ops::memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				Level1Ops::copySparseBlocToDenseBlocArray(R, D[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\trow\t" << j << "/" << nb_row_blocs_C << "\tcolumn\t" << j << "/" << nb_column_blocs_D << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B\t" << k + first_bloc_idx << endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					Level2Ops::reduceBlocByRectangularBloc(R, C[j][k], B[k + first_bloc_idx][i], dense_bloc, invert_scalars);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, D[j][i]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);
		}
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);
	TIMER_REPORT_(reduceBlocByRectangularBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}







#endif /* LEVEL3_OPS_C_ */
