/*
 * matrix-ops.C
 *
 *  Created on: 6 juil. 2012
 *      Author: martani
 */

#ifndef MATRIX_OPS_C_
#define MATRIX_OPS_C_

#include <stdlib.h>
#include <malloc.h>

#include "matrix-ops.h"
#include "matrix-utils.h"
#include "types.h"

using namespace LELA;

#ifndef DEFAULT_BLOC_HEIGHT
#define DEFAULT_BLOC_HEIGHT 0
#error "must define DEFAULT_BLOC_HEIGHT"
#endif

#ifndef DEFAULT_BLOC_WIDTH
#define DEFAULT_BLOC_WIDTH 0
#error "must define DEFAULT_BLOC_WIDTH"
#endif

/**
 * Given a bloc of dense rows (an array of arrays), Zeros its memory
 */
template <typename Element>
inline void memsetBlocToZero(Element** arr, const uint32 nb_lines, const uint32 line_size)
{
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		memset(arr[i], 0, DEFAULT_BLOC_WIDTH * sizeof(Element));
	}
}

/**
 * Copy a SparseBloc to a bloc of dense rows (an array or arrays)
 */
template <typename Ring, typename Index, typename DoubleFlatElement>
inline void copySparseBlocToDenseBlocArray(const Ring& R, 
					   const ContiguousBloc<typename Ring::Element, Index>& bloc, 
					   DoubleFlatElement** arr)
{
	//Context<Ring> ctx (R);

	//for all the rows in the bloc
#pragma loop unroll
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		//BLAS1::write(ctx, cout, bloc[i]);
		//for all the elements in the row bloc[i]
		for(uint32 j=0; j<bloc.get_row_size(i); ++j)
		{
			arr[i][bloc.pos_in_row(i, j)] = bloc.value_in_row(i, j);
		}

		/*cout << endl << "*********** dense *************" << endl;
		for(uint32 j=0; j<bloc.bloc_width(); ++j)
			cout << arr[i][j] << ", ";

		cout << endl;*/
	}
}

template <typename Ring, typename DoubleFlatElement, typename Index>
inline void copyDenseBlocArrayToSparseBloc(const Ring& R, 
					   DoubleFlatElement** arr, 
					   ContiguousBloc<typename Ring::Element, Index>& bloc)
{
	typename Ring::Element e;
	
	//this shoud be a global
	//typename Ring::Element *reduced_non_zero_elts = new typename Ring::Element [bloc.bloc_width ()];
	//TODO: pre allocate memory
	bloc.clear ();
	
	//for each row in the bloc
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		//tmp = typename SparseBloc<typename Ring::Element>::Row ();

		uint32 nb_elts = 0;
		for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
		{
			ModularTraits<typename Ring::Element>::reduce (e, arr[i][j], R._modulus);

			if (!R.isZero(e))
			{
				++nb_elts;
				bloc.push_back_value (e);
				bloc.push_back_pos (j);
			}
		}

		bloc.set_row_size (i, nb_elts);
	}
	
	bloc.construct_rows_starting_offset ();
}

/*template <typename Ring, typename Row >
inline void copyDenseArrayToSparseRow(const Ring& R, uint64* arr, const uint32 line_size, Row& v)
{

	//Row tmp;

	v.clear ();
	//if(v.capacity () < 128)
		//v.reserve (128);

	//register uint16 l;
	register typename Ring::Element e;
	for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
	{
		//ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
		//l = arr[j] >> 16;
		//l = R._modulus - l;

		e = arr[j] % R._modulus;
		if(!R.isZero(e))
			v.push_back(typename Row::value_type (j, e));
	}

	//v.swap(tmp);
}*/

/*template <typename Ring>
inline void reduceDenseBlocModulo(const Ring& R, uint64** arr, const uint32 nb_lines, const uint32 line_size)
{
	typename Ring::Element e;

	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
		{
			ModularTraits<typename Ring::Element>::reduce (e, arr[i][j], R._modulus);
			arr[i][j] = (uint64)e;
		}
	}
}*/

template <typename Ring>
inline void reduceDenseArrayInRing(const Ring& R, uint64* arr)
{
	typename Ring::Element e;

	for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
	{
		ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
		arr[j] = (uint64)e;
	}
}


/**
 * acc <- acc + a.bloc[row]
 * Note: caller must ensure that acc is large enough (>= x.size())
 */
template <typename Index>
inline void axpy(const ContiguousBloc<uint16, Index>& bloc,
		 const uint16 row,
		 const uint16 a, 
		 uint64* arr) __attribute__((always_inline));

template <typename Index>
inline void axpy(const ContiguousBloc<uint16, Index>& bloc,
		 const uint16 row,
		 const uint16 a, 
		 uint64* arr)
{
	const Index sz = bloc.get_row_size(row);
//	for (uint32 i = 0; i < sz; ++i)
//	{
//		arr[bloc.pos_in_row(row, i)] += (uint32) a * bloc.value_in_row(row, i);
//	}

#define STEP 32

	register uint32 x=0;
	while(x<sz % STEP)
	{
		arr[bloc.pos_in_row(row, x)] += (uint32)a * bloc.value_in_row(row, x);
		++x;
	}

	for(x=sz % STEP; x<sz; x+=STEP)
	{
		for(uint8 t=0; t<STEP; ++t)
		{
			arr[bloc.pos_in_row(row, x+t)] += (uint32)a * bloc.value_in_row(row, x+t);
		}
			//R.axpyin(tmpDenseArray[(row_it_B+t)->first], 	 Cv, (row_it_B+t)->second);
	}

#undef STEP
}

inline void axpy(const uint64 *arr1, const uint16 a, uint64* arr2) __attribute__((always_inline));

//TODO: warning: arr1 must be already reduced modulo modulus in Ring
inline void axpy(const uint64 *arr1, const uint16 a, uint64* arr2)
{
	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; ++i)
	{
		arr2[i] += arr1[i] * a;
	}
}


TIMER_DECLARE_(reduceDenseBlocModulo);
TIMER_DECLARE_(axpyInTriangular);
TIMER_DECLARE_(axpyInRectangular);
TIMER_DECLARE_(reduceDenseArrayInRing);


template <typename Index>
inline void reduceBlocByRectangularBloc(const Modular<uint16>& R,
		const ContiguousBloc<uint16, Index>& bloc_A,
		const ContiguousBloc<uint16, Index>& bloc_B,
		uint64 **Bloc_acc)
{
	//cout << "REDUCE RECTANGULAR nb_lines " << bloc_A.size() << endl;

	if(bloc_A.empty() || bloc_B.empty())
		return;

	//for all the rows in Bloc_acc (same as for Bloc_A)
	for(uint16 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		//if(bloc_A[i].empty())
			//continue;

		//for all element in row bloc_A[i]
		for(uint16 j=0; j < bloc_A.get_row_size(i); ++j)
		{
			register uint16 Cv;
			register typename ContiguousBloc<uint16, Index>::IndexType Cp;

			R.copy(Cv, bloc_A.value_in_row(i, j));
			Cp = bloc_A.pos_in_row(i, j);
			R.negin(Cv);

			TIMER_START_(axpyInRectangular);
				axpy(bloc_B, Cp, Cv, Bloc_acc[i]);
			TIMER_STOP_(axpyInRectangular);
		}
	}
}


/**
 * Reduce the rows inside the bloc by themselves
 */
template<typename Index>
inline void reduceBlocByTriangularBloc(const Modular<uint16>& R,
									   const ContiguousBloc<uint16, Index>& bloc_A,
									   uint64** Bloc_acc)
{
	for(uint32 i=1; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		TIMER_START_(reduceDenseArrayInRing);
			reduceDenseArrayInRing(R, Bloc_acc[i-1]);
		TIMER_STOP_(reduceDenseArrayInRing);

		const Index sz = bloc_A.get_row_size(i) ;
		for(int j=0; j < (int)sz - 1; ++j)
		{
			register uint16 Av;
			register typename ContiguousBloc<uint16, Index>::IndexType Ap;

			R.copy(Av, bloc_A.value_in_row(i, j));
			Ap = bloc_A.pos_in_row(i, j);
			R.negin(Av);
			
			//lela_check(Ap < i);	//should point to a a row previous to the ith one in the bloc

			TIMER_START_(axpyInTriangular);
				axpy(Bloc_acc[Ap], Av, Bloc_acc[i]);
			TIMER_STOP_(axpyInTriangular);
		}
	}

	TIMER_START_(reduceDenseArrayInRing);
		reduceDenseArrayInRing(R, Bloc_acc[DEFAULT_BLOC_HEIGHT-1]);
	TIMER_STOP_(reduceDenseArrayInRing);

	//cout << "REDUCE TRIANGULAR DONE\n";
	
}

template <typename Index>
void MatrixOps::reducePivotsByPivots(const Modular<uint16>& R,
		const SparseBlocMatrix<ContiguousBloc<uint16, Index> >& A, 
		SparseBlocMatrix<ContiguousBloc<uint16, Index> >& B)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Bloc version" << std::endl;

	typedef SparseBlocMatrix<ContiguousBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		//dense_bloc[i] = new uint64[B.bloc_width()];
		dense_bloc[i] = (uint64 *) memalign(64, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseBlocModulo);
	TIMER_RESET_(axpyInTriangular);
	TIMER_RESET_(axpyInRectangular);



	//for all columns of blocs of B
	for (uint32 i = 0; i < (uint32) std::ceil((double) B.coldim() / B.bloc_width()); ++i)
	{
		//report << "Column B " << i << endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 j = 0;
				j < (uint32) std::ceil((double) A.rowdim() / A.bloc_height());
				++j)
		{
			uint32 first_bloc_idx, last_bloc_idx;
			first_bloc_idx = A.FirstBlocsColumIndexes[j] / A.bloc_width();
			last_bloc_idx = min(A[j].size () - 1, j);

			//report << "\tRow A " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j
			TIMER_START_(memsetBlocToZero);
				memsetBlocToZero(dense_bloc, B.bloc_height(), B.bloc_width());
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				copySparseBlocToDenseBlocArray(R, B[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "\trow\t" << j << std::ends;
#endif
			//for all the blocs in the current row of A
			for (int k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in A and B " << k + first_bloc_idx << endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					reduceBlocByRectangularBloc(R, A[j][k], B[k + first_bloc_idx][i], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			//report << "reduceBlocByRectangularBloc() DONE" << endl;

			TIMER_START_(reduceBlocByTriangularBloc);
				reduceBlocByTriangularBloc(R, A[j][last_bloc_idx], dense_bloc);
			TIMER_STOP_(reduceBlocByTriangularBloc);
			//report << "reduceBlocByTriangularBloc DONE" << endl;

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);
			//report << "copyDenseBlocArrayToSparseBloc DONE" << endl;

		}
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);

	report << endl;
	TIMER_REPORT_(reduceBlocByRectangularBloc);
	TIMER_REPORT_(axpyInRectangular);

	report << endl;
	TIMER_REPORT_(reduceBlocByTriangularBloc);
	TIMER_REPORT_(axpyInTriangular);

	report << endl;
	TIMER_REPORT_(reduceDenseArrayInRing);
}

template <typename Index>
void MatrixOps::reduceNonPivotsByPivots(const Modular<uint16>& R,
		const SparseBlocMatrix<ContiguousBloc<uint16, Index> >& C, 
		const SparseBlocMatrix<ContiguousBloc<uint16, Index> >& B,
		SparseBlocMatrix<ContiguousBloc<uint16, Index> >& D)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Bloc version" << std::endl;

	typedef Modular<uint16> Ring;
	typedef SparseBlocMatrix<ContiguousBloc<typename Ring::Element, uint16> > Matrix;

	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(D.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(C.rowdim(), D.rowdim());
	check_equal_or_raise_exception(C.coldim(), B.rowdim());
	check_equal_or_raise_exception(B.coldim(), D.coldim());



	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		//dense_bloc[i] = new uint64[D.bloc_width()];
		dense_bloc[i] = (uint64 *) memalign(64, D.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseBlocModulo);
	TIMER_RESET_(axpyInTriangular);
	TIMER_RESET_(axpyInRectangular);

	//for all columns of blocs of B
	for (uint32 i = 0; i < (uint32) std::ceil((double) D.coldim() / D.bloc_width()); ++i)
	{
		//report << "Colum D" << i << endl;
		//for all rows of blocs in C
		for(uint32 j=0; j<(uint32)std::ceil((double)C.rowdim() / C.bloc_height()); ++j)
		{
			uint32 first_bloc_idx, last_bloc_idx;
			first_bloc_idx = C.FirstBlocsColumIndexes[j] / C.bloc_width();
			last_bloc_idx = C[j].size ();

			//report << "\tRow C " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j
			TIMER_START_(memsetBlocToZero);
				memsetBlocToZero(dense_bloc, D.bloc_height(), D.bloc_width());
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				copySparseBlocToDenseBlocArray(R, D[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "\trow\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (int k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B " << k + first_bloc_idx << endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					reduceBlocByRectangularBloc(R, C[j][k], B[k + first_bloc_idx][i], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				copyDenseBlocArrayToSparseBloc(R, dense_bloc, D[j][i]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);
		}
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	TIMER_REPORT_(memsetBlocToZero);
	TIMER_REPORT_(copySparseBlocToDenseBlocArray);
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);

	report << endl;
	TIMER_REPORT_(reduceBlocByRectangularBloc);
	TIMER_REPORT_(axpyInRectangular);

	report << endl;
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}




#endif /* MATRIX_OPS_H_ */




