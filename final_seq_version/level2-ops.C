/*
 * level2-ops.C
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */



#ifndef LEVEL2_OPS_C_
#define LEVEL2_OPS_C_

#include "level2-ops.h"
#include "level1-ops.h"

#include "consts-macros.h"

#include <omp.h>

using namespace LELA;

void Level2Ops::DenseScalMulSub__one_row__array_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint64 *arr_source,
		uint64 *arr1,
		uint64 *arr2)
{
	register uint32 v__;
	//register uint32 v2__;

//	for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i++)
//	{
//		v1__ = arr_source[i] & 0x000000000000ffff;
//		//v1__ *= av1_col1;
//		arr1[i] += v1__ * av1_col1;;
//
//		//v2__ = arr_source[i] & 0x000000000000ffff;
//		//v2__ *= av2_col1;
//		arr2[i] += v1__ * av2_col1;;
//
//	}

	//XXX: USE THIS
//	for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
//	{
//		for(uint32 t=0; t<UNROLL_STEP__64; t++)
//		{
//			v1__ = arr_source[i+t] & 0x000000000000ffff;
//
//			arr1[i+t] += v1__ * av1_col1;;
//			arr2[i+t] += v1__ * av2_col1;;
//		}
//
//	}

	if(av1_col1 == 0)
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
		{
			for(uint32 t=0; t<UNROLL_STEP__64; t++)
			{
				v__ = arr_source[i+t] & 0x000000000000ffff;
				arr2[i+t] += v__ * av2_col1;;
			}
		}
	else if (av2_col1 == 0)
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
		{
			for(uint32 t=0; t<UNROLL_STEP__64; t++)
			{
				v__ = arr_source[i+t] & 0x000000000000ffff;
				arr1[i+t] += v__ * av1_col1;;
			}
		}
	else
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
		{
			for(uint32 t=0; t<UNROLL_STEP__64; t++)
			{
				v__ = arr_source[i+t] & 0x000000000000ffff;

				arr1[i+t] += v__ * av1_col1;;
				arr2[i+t] += v__ * av2_col1;;
			}
		}

}


void Level2Ops::DenseScalMulSub__two_rows__array_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const uint64 *arr_source1,
		const uint64 *arr_source2,
		uint64 *arr1,
		uint64 *arr2)
{
	if(av1_col1 == 0 && av2_col1 == 0)
	{
		DenseScalMulSub__one_row__array_array(
				av1_col2,
				av2_col2,
				arr_source2,
				arr1,
				arr2);

		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		DenseScalMulSub__one_row__array_array(
				av1_col1,
				av2_col1,
				arr_source1,
				arr1,
				arr2);

		return;
	}

	register uint32 v1__, v2__;

//	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i++)
//	{
//		v1__ = arr_source1[i] & 0x000000000000ffff;
//		v2__ = arr_source2[i] & 0x000000000000ffff;
//
//		arr1[i] += v1__ * av1_col1;
//		arr1[i] += v2__ * av1_col2;
//
//		arr2[i] += v1__ * av2_col1;
//		arr2[i] += v2__ * av2_col2;
//	}

	//XXX: USE THIS
	for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
	{
		for(uint32 t=0; t<UNROLL_STEP__64; t++)
		{
			v1__ = arr_source1[i+t] & 0x000000000000ffff;
			v2__ = arr_source2[i+t] & 0x000000000000ffff;

			arr1[i+t] += v1__ * av1_col1;
			arr1[i+t] += v2__ * av1_col2;

			arr2[i+t] += v1__ * av2_col1;
			arr2[i+t] += v2__ * av2_col2;
		}
	}
}


template <typename Index>
void Level2Ops::SparseScalMulSub__one_row__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2)
{
	const uint32 N = v.size();
	uint32 i = 0;
	/*for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		val1 = v.at_unchecked(line, i);

		arr1[idx] += (uint32) av1_col1 * val1;
		arr2[idx] += (uint32) av2_col1 * val1;
	}*/
	const Index *p_idx = (N != 0) ? v.IndexData.getStartingPointer () : NULL;
	const uint16 *p_val = (N != 0) ? v.ValuesData.getStartingPointer () : NULL;
	p_val += line;

	register uint32 v__;
	register uint32 idx;

	///XXX:USE THIS
//	for(; i < N; i++)
//	{
//		idx = p_idx[i];
//		v__ = p_val[i*2];
//
//		arr1[idx] += v__ * av1_col1 ;
//		arr2[idx] += v__ * av2_col1;
//	}

	//XXX: UNROLLING NOT EFFICIENT FOR THIS LOOP
//	for(; i < ROUND_DOWN(N, UNROLL_STEP__16); i+=UNROLL_STEP__16)
//	{
//		for(uint32 t=0; t<UNROLL_STEP__16; t++)
//		{
//			idx = p_idx[i+t];
//			v__ = p_val[(i+t)*2];
//
//			arr1[idx] += v__ * av1_col1 ;
//			arr2[idx] += v__ * av2_col1;
//		}
//	}
//
//	for(; i < N; i++)
//	{
//		idx = p_idx[i];
//		v__ = p_val[i*2];
//
//		arr1[idx] += v__ * av1_col1 ;
//		arr2[idx] += v__ * av2_col1;
//	}

	if (av1_col1 != 0 && av2_col1 != 0) //cannot both be 0
	{
		for(; i < N; i++)
		{
			idx = p_idx[i];
			v__ = p_val[i*2];

			arr1[idx] += v__ * av1_col1 ;
			arr2[idx] += v__ * av2_col1;
		}
	}
	else if (av1_col1 != 0)
	{
		for(; i < N; i++)
		{
			idx = p_idx[i];
			v__ = p_val[i*2];

			v__ *= av1_col1;
			arr1[idx] += v__;
		}
	}
	else //av2_col1 != 0
	{
		for(; i < N; i++)
		{
			idx = p_idx[i];
			v__ = p_val[i*2];

			v__ *= av2_col1;
			arr2[idx] += v__;
		}
	}
}


template <typename Index>
void Level2Ops::DenseScalMulSub__one_row__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2)
{
//	check_equal_or_raise_exception(v.is_sparse(), false);

	const uint16 *p_val = v.ValuesData.getStartingPointer ();
	p_val += line;

	register uint32 v__;

//	for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i++)
//	{
//		v__ = p_val[(i)*2];
//
//		arr1[i] += v__ * av1_col1;
//		arr2[i] += v__ * av2_col1;
//	}

	//XXX: UNROLLING EFFICIENT FOR THIS LOOP
//	for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__16)
//	{
//		for(Index t=0; t<UNROLL_STEP__16; t++)
//		{
//			v__ = p_val[(i+t)*2];
//
//			arr1[i+t] += v__ * av1_col1;
//			arr2[i+t] += v__ * av2_col1;
//		}
//	}
	if(av1_col1 != 0 && av2_col1 != 0)
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__16)
		{
			for(Index t=0; t<UNROLL_STEP__16; t++)
			{
				v__ = p_val[(i+t)*2];

				arr1[i+t] += v__ * av1_col1;
				arr2[i+t] += v__ * av2_col1;
			}
		}
	if(av1_col1 == 0)
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__16)
		{
			for(Index t=0; t<UNROLL_STEP__16; t++)
			{
				v__ = p_val[(i+t)*2];
				arr2[i+t] += v__ * av2_col1;
			}
		}
	else if(av2_col1 == 0)
		for(uint32 i=0; i<DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__16)
		{
			for(Index t=0; t<UNROLL_STEP__16; t++)
			{
				v__ = p_val[(i+t)*2];
				arr1[i+t] += v__ * av1_col1;
			}
		}
}


template <typename Index>
void Level2Ops::DenseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2)
{
//	check_equal_or_raise_exception(v.is_sparse(), false);

	if(av1_col1 == 0 && av2_col1 == 0)
	{
		Level2Ops::DenseScalMulSub__one_row__vect_array(
				av1_col2,
				av2_col2,
				v,
				1,
				arr1,
				arr2);

		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		DenseScalMulSub__one_row__vect_array(
				av1_col1,
				av2_col1,
				v,
				0,
				arr1,
				arr2);

		return;
	}

	const uint16 *p_val = v.ValuesData.getStartingPointer ();

	register uint32 v1__, v2__;

//	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i++)
//	{
//		v1__ = p_val[(i)*2];
//		v2__ = p_val2[(i)*2];
//
//		arr1[i] += v1__ * av1_col1;
//		arr1[i] += v2__ * av1_col2;
//
//		v1__ *= av2_col1;
//		v2__ *= av2_col2;
//		arr2[i] += v1__;
//		arr2[i] += v2__;
//	}

	//XXX: LOOP UNROLL with 16 EFFICIENT FOR THIS LOOP
	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__16)
	{
// 		__builtin_prefetch(p_val, 	__PREFETCH_READ, __PREFETCH_LOCALITY_NO_LOCALITY);
// 		__builtin_prefetch(p_val+32, 	__PREFETCH_READ, __PREFETCH_LOCALITY_NO_LOCALITY);
// 		__builtin_prefetch(arr1+i, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
// 		__builtin_prefetch(arr1+i+8, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
// 		__builtin_prefetch(arr1+i+16, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
// 		__builtin_prefetch(arr2+i, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
// 		__builtin_prefetch(arr2+i+8, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
// 		__builtin_prefetch(arr2+i+16, 	__PREFETCH_WRITE, __PREFETCH_LOCALITY_LOW);
		
		for(uint32 t=0; t<UNROLL_STEP__16; t++)
		{
			v1__ = p_val[(i+t)*2];
			v2__ = p_val[(i+t)*2+1];

			arr1[i+t] += v1__ * av1_col1;
			arr1[i+t] += v2__ * av1_col2;

			v1__ *= av2_col1;
			v2__ *= av2_col2;
			arr2[i+t] += v1__;
			arr2[i+t] += v2__;
		}
	}
}

template <typename Index>
void Level2Ops::SparseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2)
{
	if(av1_col1 == 0 && av2_col1 == 0)
	{
		SparseScalMulSub__one_row__vect_array(
				av1_col2,
				av2_col2,
				v,
				1,
				arr1,
				arr2);

		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		SparseScalMulSub__one_row__vect_array(
				av1_col1,
				av2_col1,
				v,
				0,
				arr1,
				arr2);

		return;
	}

	const uint32 N = v.size ();
	uint32 i = 0;
	const Index *p_idx = (N != 0) ? v.IndexData.getStartingPointer () : NULL;
	const uint16 *p_val = (N != 0) ? v.ValuesData.getStartingPointer () : NULL;
	const uint16 *p_val2 = p_val + 1;

	register uint32 v1__, v2__;
	register uint32 idx;



	for(; i < ROUND_DOWN(N, UNROLL_STEP__16); i+=UNROLL_STEP__16)
	{
		for(uint32 t=0; t<UNROLL_STEP__16; ++t)
		{
			idx = p_idx[i+t];
			v1__ = p_val[(i+t)*2];
			v2__ = p_val2[(t+i)*2];

			arr1[idx] += av1_col1 * v1__;
			arr1[idx] += av1_col2 * v2__;

			arr2[idx] += av2_col1 * v1__;
			arr2[idx] += av2_col2 * v2__;
		}
	}
	for(; i < N; i++)
	{
		idx = p_idx[i];
		v1__ = p_val[i*2];
		v2__ = p_val2[i*2];

		arr1[idx] += av1_col1 * v1__;
		arr1[idx] += av1_col2 * v2__;

		arr2[idx] += av2_col1 * v1__;
		arr2[idx] += av2_col2 * v2__;
	}

	/*for (uint32 i = 0; i < v.size(); ++i)
	{
		idx = p_idx[i];

		val1 = v.at_unchecked(0, i);
		val2 = v.at_unchecked(1, i);

		arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
		arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
	}

	return;*/
}


template <typename Index>
void Level2Ops::DenseScalMulSub__one_row__vect_array__variable_size(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2,
		const int start_offset)
{
	if(start_offset < 0)
		return;

	const uint16 *p_val = v.ValuesData.getStartingPointer ();
	p_val += line;

	register uint32 v__;
	uint32 i;

	//XXX: UNROLLING EFFICIENT FOR THIS LOOP
	if(av1_col1 != 0 && av2_col1 != 0)
	{
		for(i=start_offset; i < v.size(); i++)
		{
			v__ = p_val[i*2];

			arr1[i] += v__ * av1_col1;
			arr2[i] += v__ * av2_col1;
		}
	}
	if(av1_col1 == 0)
	{
		for(i=start_offset; i < v.size(); i++)
		{
			v__ = p_val[i*2];
			arr2[i] += v__ * av2_col1;
		}
	}
	else if(av2_col1 == 0)
	{
		for(i=start_offset; i < v.size(); i++)
		{
			v__ = p_val[i*2];

			arr1[i] += v__ * av1_col1;
		}
	}

}


template <typename Index>
void Level2Ops::DenseScalMulSub__two_rows__vect_array__variable_size(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2,
		const int start_offset)
{
	if(start_offset < 0)
	{
		return;
	}

//	check_equal_or_raise_exception(v.is_sparse(), false);

 	if(av1_col1 == 0 && av2_col1 == 0)
 	{
		Level2Ops::DenseScalMulSub__one_row__vect_array__variable_size(
				av1_col2,
				av2_col2,
				v,
				1,
				arr1,
				arr2,
				start_offset);
 		return;
 	}

 	if(av1_col2 == 0 && av2_col2 == 0)
 	{
 		DenseScalMulSub__one_row__vect_array__variable_size(
				av1_col1,
				av2_col1,
				v,
				0,
				arr1,
				arr2,
				start_offset);

 		return;
 	}

	const uint32 N = v.size();
	const uint16 *p_val = v.ValuesData.getStartingPointer ();

	register uint32 v1__, v2__;
	for (uint32 i = start_offset; i < N; i++)
	{
		v1__ = p_val[i * 2];
		v2__ = p_val[i * 2 + 1];

		arr1[i] += av1_col1 * v1__;
		arr1[i] += av1_col2 * v2__;

		arr2[i] += av2_col1 * v1__;
		arr2[i] += av2_col2 * v2__;
	}
}



template <typename Index>
void Level2Ops::reduceBlocByRectangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		const SparseMultilineBloc<uint16, Index>& bloc_B,
		uint64 **Bloc_acc,
		bool invert_scalars)
{
	typedef Modular<uint16> Ring;

	if(bloc_A.empty() || bloc_B.empty())
			return;

	for(int i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		uint8 is_sparse = 0;

		if(bloc_A[i].is_sparse ())
			is_sparse = 1;
		else
			is_sparse = 0;

		const Index N = is_sparse == 1 ? bloc_A[i].size() : DEFAULT_BLOC_WIDTH;

		for (uint32 j = 0; j < N; ++j)
		{
			const Index Ap1 = (is_sparse == 1 ? bloc_A[i].IndexData[j] : j);

			//R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
			register uint32 Av1_col1 = bloc_A[i].at_unchecked(0, j);
			//R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));
			register uint32 Av2_col1 = bloc_A[i].at_unchecked(1, j);

			if(invert_scalars)
			{
				if (Av1_col1 != 0)
					Av1_col1 = (uint32)R._modulus - Av1_col1;
				if (Av2_col1 != 0)
					Av2_col1 = (uint32)R._modulus - Av2_col1;
			}

			if (((Ap1 % 2) == 0) && (j < (uint32)(N - 1)))
			{
				const Index Ap2 = (is_sparse == 1 ? bloc_A[i].IndexData[j+1] : j+1);
				if (Ap2 == Ap1 + 1) //axpy 2 ROWS
				{
// 					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
// 					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
					register uint32 Av1_col2 = bloc_A[i].at_unchecked(0, j+1);
					register uint32 Av2_col2 = bloc_A[i].at_unchecked(1, j+1);
					
					if(invert_scalars)
					{
						if (Av1_col2 != 0)
							Av1_col2 = (uint32)R._modulus - Av1_col2;
						if (Av2_col2 != 0)
							Av2_col2 = (uint32)R._modulus - Av2_col2;
					}

					++j;

					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
						SparseScalMulSub__two_rows__vect_array(
								Av1_col1,
								Av2_col1,
								Av1_col2,
								Av2_col2,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
					else
					{
						DenseScalMulSub__two_rows__vect_array(
								Av1_col1,
								Av2_col1,
								Av1_col2,
								Av2_col2,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
				}
				else	//axpy ONE ROW
				{
					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
						SparseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
					else
					{
						DenseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
				}
			}
			else	//axpy ONE ROW
			{
					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
							SparseScalMulSub__one_row__vect_array(
									Av1_col1,
									Av2_col1,
									bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2],
									Bloc_acc[i * 2 + 1]);
					}
					else
					{
						DenseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}

			}
			//__builtin_prefetch (&(bloc_A[i].IndexData[j+1]), __PREFETCH_READ, __PREFETCH_LOCALITY_LOW);
			//__builtin_prefetch (&(bloc_A[i].ValuesData[(j+1)*NB_ROWS_PER_MULTILINE]), __PREFETCH_READ, __PREFETCH_LOCALITY_MODERATE);
		}
	}

}


/*
 * Note: All rows of bloc A are supposed to be sparse.
 */
template<typename Index>
void Level2Ops::reduceBlocByTriangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64** Bloc_acc,
		bool invert_scalars)
{
	typedef Modular<uint16> Ring;

	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		if(bloc_A[i].empty ())
			continue;

		//if(bloc_A[i].is_sparse ())
		//{
			int last_idx = -1;
			if(bloc_A[i].at_unchecked(1, bloc_A[i].size()-1) == 0)
				last_idx = (int)bloc_A[i].size()-1;
			else
				last_idx = (int)bloc_A[i].size()-2;

			typename Ring::Element Av1_col1, Av2_col1;
			uint32 Ap1;
			typename Ring::Element Av1_col2, Av2_col2;
			uint32 Ap2;

			for (int j = 0; j < last_idx; ++j)	//skip first two elements
			{
				Ap1 = (uint32)bloc_A[i].IndexData[j];
				R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
				R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));

				if(invert_scalars)
				{
					if (Av1_col1 != 0)
						R.negin(Av1_col1);
					if (Av2_col1 != 0)
						R.negin(Av2_col1);
				}

				//if(Ap1 >= i*2)
				//	throw std::runtime_error ("Index pointing to an out of range line.");

				if (Ap1 % 2 == 0 && j < last_idx - 1)
				{
					Ap2 = bloc_A[i].IndexData[j+1];
					if (Ap2 == Ap1+1)				//axpy TWO ARRAYS
					{
						R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
						R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
				
						if(invert_scalars)
						{
							if (Av1_col2 != 0)
								R.negin(Av1_col2);
							if (Av2_col2 != 0)
								R.negin(Av2_col2);
						}

						++j;

						DenseScalMulSub__two_rows__array_array(
								Av1_col1,
								Av2_col1,
								Av1_col2,
								Av2_col2,
								Bloc_acc[Ap1],
								Bloc_acc[Ap1+1],
								Bloc_acc[i*2],
								Bloc_acc[i*2+1]);
					}
					else	//axpy ONE ARRAY
					{
						DenseScalMulSub__one_row__array_array(
								Av1_col1,
								Av2_col1,
								Bloc_acc[Ap1],
								Bloc_acc[i*2],
								Bloc_acc[i*2+1]);
					}
				}
				else	//axpy ONE ARRAY
				{
					DenseScalMulSub__one_row__array_array(
							Av1_col1,
							Av2_col1,
							Bloc_acc[Ap1],
							Bloc_acc[i*2],
							Bloc_acc[i*2+1]);
				}
			}

			Level1Ops::reduceDenseArrayModulo(R, Bloc_acc[i*2]);

			if (bloc_A[i].size() > 1) //reduce lines within the same multiline
			{
				int j=bloc_A[i].size()-2;
				Ap1 = bloc_A[i].IndexData[j];
				R.copy(Av1_col1, bloc_A[i].at_unchecked(1, j));

				if (Av1_col1 != 0)
				{
					R.negin(Av1_col1);

					const Index offset1 = i*2+1;
					const Index offset2 = i*2;

					for (uint32 t = 0; t < DEFAULT_BLOC_WIDTH; ++t)
						Bloc_acc[offset1][t] += (uint32) Av1_col1 * Bloc_acc[offset2][t];
				}

			}

			Level1Ops::reduceDenseArrayModulo(R, Bloc_acc[i*2+1]);

	}  //for i
}


template<typename Index>
void Level2Ops::reduceBlocByRectangularBloc_C(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_C,
		const SparseMultilineBloc<uint16, Index>& bloc_A, uint64 **bloc_dense)
{
	typedef Modular<uint16> Ring;

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
	{
		uint8 is_sparse = 0;
		
		if (bloc_C[i].empty())
			continue;

		if(bloc_C[i].is_sparse ())
			is_sparse = 1;
		else
			is_sparse = 0;

		const Index N = is_sparse == 1 ? bloc_C[i].size() : DEFAULT_BLOC_WIDTH;

		for (int j = 0; j < N; ++j)
		{
			const Index Cp1 = (is_sparse == 1 ? bloc_C[i].IndexData[j] : j);
			register typename Ring::Element Cv1_col1 = bloc_C[i].at_unchecked(0, j);
			register typename Ring::Element Cv2_col1 = bloc_C[i].at_unchecked(1, j);

			if(Cv1_col1 == 0 && Cv2_col1 == 0)
				continue;

			if (bloc_A[Cp1 / NB_ROWS_PER_MULTILINE].empty())
				continue;

			if (Cp1 % 2 == 0 && j < N - 1)
			{
				const Index Cp2 = (is_sparse == 1 ? bloc_C[i].IndexData[j+1] : j+1);

				if(Cp1 == Cp2-1)	//succesive values
				{
					register typename Ring::Element Cv1_col2 = bloc_C[i].at_unchecked(0, j+1) % R._modulus;
					register typename Ring::Element Cv2_col2 = bloc_C[i].at_unchecked(1, j+1) % R._modulus;

					SparseScalMulSub__two_rows__vect_array(
							Cv1_col1,
							Cv2_col1,
							Cv1_col2,
							Cv2_col2,
							bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
							bloc_dense[i * 2],
							bloc_dense[i * 2 + 1]);

					++j;
				}
				else
				{
					SparseScalMulSub__one_row__vect_array(
							Cv1_col1,
							Cv2_col1,
							bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
							Cp1 % NB_ROWS_PER_MULTILINE,
							bloc_dense[i * 2],
							bloc_dense[i * 2 + 1]);
				}
			}
			else
			{
				SparseScalMulSub__one_row__vect_array(
						Cv1_col1,
						Cv2_col1,
						bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
						Cp1 % NB_ROWS_PER_MULTILINE,
						bloc_dense[i * 2],
						bloc_dense[i * 2 + 1]);
			}
		}
	}

}


template <typename Index>
void Level2Ops::reduceBlocByTriangularBloc_C(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64 **bloc_dense)
{
	typedef Modular<uint16> Ring;

	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		typename Ring::Element Av1_col1, Av2_col1,
								Av1_col2, Av2_col2;
		uint32 tmp;
		Index Ap1;
		Index sz;

		for (int j = DEFAULT_BLOC_WIDTH-1; j >= 0; --j)	//skip first two elements
		{
			Ap1 = j;
			Av1_col1 = bloc_dense[i*2][j] % R._modulus;
			Av2_col1 = bloc_dense[i*2+1][j] % R._modulus;

			if(Av1_col1 == 0 && Av2_col1 == 0)
				continue;

			if(bloc_A[Ap1/NB_ROWS_PER_MULTILINE].empty ())
			{
				//std::cout << "PROBLEM bloc_A TRIANGULAR\n";
				continue;
			}

			//neg
			if (Av1_col1 != 0)
				Av1_col1 = R._modulus - Av1_col1;
			if (Av2_col1 != 0)
				Av2_col1 = R._modulus - Av2_col1;

			if (Ap1 % 2 == 1)
			{
				Av1_col2 = bloc_dense[i*2][Ap1-1] % R._modulus;
				Av2_col2 = bloc_dense[i*2+1][Ap1-1] % R._modulus;

				//TODO: in case a multiline has only one entry of the pivot 1
				sz  = bloc_A[Ap1 / NB_ROWS_PER_MULTILINE].size ();

				const uint32 val1 = bloc_A[Ap1 / NB_ROWS_PER_MULTILINE].at_unchecked(1, sz-2);

				tmp = Av1_col2 + (uint32)Av1_col1*val1;
				R.init(Av1_col2, tmp);

				tmp = Av2_col2 + (uint32)Av2_col1*val1;
				R.init(Av2_col2, tmp);

				//neg
				if (Av1_col2 != 0)
					Av1_col2 = R._modulus - Av1_col2;
				if (Av2_col2 != 0)
					Av2_col2 = R._modulus - Av2_col2;

				SparseScalMulSub__two_rows__vect_array(
						Av1_col2,	//inversed!
						Av2_col2,	//AvX_col2 elements are multiplied by the 0th row of the Multiline
						Av1_col1,
						Av2_col1,
						bloc_A[Ap1 / NB_ROWS_PER_MULTILINE],
						bloc_dense[i * 2],
						bloc_dense[i * 2 + 1]);

				bloc_dense[i * 2][Ap1] = Av1_col1;
				bloc_dense[i * 2 + 1][Ap1] = Av2_col1;

				bloc_dense[i * 2][Ap1 - 1] = Av1_col2;
				bloc_dense[i * 2 + 1][Ap1 - 1] = Av2_col2;

				--j;
			}
			else
			{
				SparseScalMulSub__one_row__vect_array(
						Av1_col1,
						Av2_col1,
						bloc_A[Ap1 / NB_ROWS_PER_MULTILINE],
						Ap1 % NB_ROWS_PER_MULTILINE,
						bloc_dense[i * 2],
						bloc_dense[i * 2 + 1]);

				bloc_dense[i * 2][Ap1] = Av1_col1;
				bloc_dense[i * 2 + 1][Ap1] = Av2_col1;
			}
		}
	}

}




#endif /* LEVEL2_OPS_C_ */
