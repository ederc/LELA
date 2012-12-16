/*
 * matrix-ops.C
 *
 *  Created on: 6 juil. 2012
 *      Author: martani
 */

#ifndef MATRIX_OPS_C_
#define MATRIX_OPS_C_

#include <stdlib.h>
#include <sys/time.h>

#include "matrix-ops.h"
#include "matrix-utils.h"

#include "matrix-op-m.C"

using namespace LELA;
using namespace std;

#ifndef DEFAULT_BLOC_HEIGHT
#define DEFAULT_BLOC_HEIGHT 0
#error "must define DEFAULT_BLOC_HEIGHT"
#endif

#ifndef DEFAULT_BLOC_WIDTH
#define DEFAULT_BLOC_WIDTH 0
#error "must define DEFAULT_BLOC_WIDTH"
#endif

#ifndef NB_ROWS_PER_MULTILINE
#define NB_ROWS_PER_MULTILINE 0
#error "must define NB_ROWS_PER_MULTILINE"
#endif

#define UNROLL_STEP__64		16
#define UNROLL_STEP__16		16



#define __PREFETCH_WRITE	1
#define __PREFETCH_READ		0
#define __PREFETCH_LOCALITY_NO_LOCALITY	0
#define __PREFETCH_LOCALITY_LOW	1
#define __PREFETCH_LOCALITY_MODERATE	2
#define __PREFETCH_LOCALITY_HIGH	3

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

/**
 * Given a bloc of dense rows (an array of arrays), Zeros its memory
 */
template <typename Element>
static inline void memsetToZero(Element** arr)
{
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		memset(arr[i], 0, DEFAULT_BLOC_WIDTH * sizeof(Element));
	}
}

template <typename Element>
static inline void memsetToZero(Element** arr, const uint32 nb_lines, const uint32 line_size)
{
	for(uint32 i=0; i<nb_lines; ++i)
	{
		memset(arr[i], 0, line_size * sizeof(Element));
	}
}


/**
 * Copy a SparseBloc to a bloc of dense rows (an array or arrays)
 */
template<typename Ring, typename Index, typename DoubleFlatElement>
static void copySparseBlocToDenseBlocArray(const Ring& R,
		const SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		DoubleFlatElement** arr)
{
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		Index idx;
		typename Ring::Element val1, val2;

		if(bloc[i].empty ())
			continue;

		if(bloc[i].is_sparse ())
			for(uint32 j=0; j<bloc[i].size (); ++j)
			{
				idx = bloc[i].IndexData[j];
				val1 = bloc[i].at_unchecked(0, j);
				val2 = bloc[i].at_unchecked(1, j);

				arr[i*2][idx] = val1;
				arr[i*2+1][idx] = val2;
			}
		else
			for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
			{
				val1 = bloc[i].at_unchecked(0, j);
				val2 = bloc[i].at_unchecked(1, j);

				arr[i*2][j] = val1;
				arr[i*2+1][j] = val2;
			}
	}
}

template<typename Ring, typename DoubleFlatElement, typename Index>
static void copyDenseBlocArrayToSparseBloc(const Ring& R,
		DoubleFlatElement** arr,
		SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		bool reduce_in_Ring = true)
{
	typename Ring::Element e1, e2;
	MultiLineVector<typename Ring::Element, Index> tmp;

	if(reduce_in_Ring)
	{
		for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
		{
			bloc[i].clear ();
			tmp.clear ();

			for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
			{
				ModularTraits<typename Ring::Element>::reduce(e1, arr[i * 2][j], R._modulus);
				ModularTraits<typename Ring::Element>::reduce(e2, arr[i * 2 + 1][j], R._modulus);

				if (!R.isZero(e1) || !R.isZero(e2))
				{
					tmp.IndexData.push_back(j);
					tmp.ValuesData.push_back(e1);
					tmp.ValuesData.push_back(e2);
				}
			}

			if ((float)tmp.size() / (float)DEFAULT_BLOC_WIDTH < bloc.get_HYBRID_REPRESENTATION_THRESHOLD())
				bloc[i].swap(tmp);
			else
			{
				Index idx = 0;
				for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
					if (idx < tmp.size() && tmp.IndexData[idx] == j)
					{
						bloc[i].ValuesData.push_back(tmp.at_unchecked(0, idx));
						bloc[i].ValuesData.push_back(tmp.at_unchecked(1, idx));
						idx++;
					}
					else
					{
						bloc[i].ValuesData.push_back(0);
						bloc[i].ValuesData.push_back(0);
					}
			}
		}
	}
	else
	{
		for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
		{
			bloc[i].clear ();
			tmp.clear ();

			for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
			{
				ModularTraits<typename Ring::Element>::reduce(e1, arr[i * 2][j], R._modulus);
				ModularTraits<typename Ring::Element>::reduce(e2, arr[i * 2 + 1][j], R._modulus);

				if (arr[i * 2][j] != 0 || arr[i * 2+1][j] != 0)
				{
					tmp.IndexData.push_back(j);
					tmp.ValuesData.push_back((typename Ring::Element)arr[i * 2][j]);
					tmp.ValuesData.push_back((typename Ring::Element)arr[i * 2+1][j]);
				}
			}

			if ((float)tmp.size() / (float)DEFAULT_BLOC_WIDTH < bloc.get_HYBRID_REPRESENTATION_THRESHOLD())
				bloc[i].swap(tmp);
			else
			{
				Index idx = 0;
				for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
					if (idx < tmp.size() && tmp.IndexData[idx] == j)
					{
						bloc[i].ValuesData.push_back(tmp.at_unchecked(0, idx));
						bloc[i].ValuesData.push_back(tmp.at_unchecked(1, idx));
						idx++;
					}
					else
					{
						bloc[i].ValuesData.push_back(0);
						bloc[i].ValuesData.push_back(0);
					}
			}
		}
	}
}


template<typename Ring, typename DoubleFlatElement, typename Index>
static void copyDenseBlocArrayToSparseBlocRTL____(const Ring& R,
		DoubleFlatElement** arr,
		SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		bool reduce_in_Ring = true)
{
	typename Ring::Element e1, e2;
	MultiLineVector<typename Ring::Element, Index> tmp;

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
	{
		bloc[i].clear();
		tmp.clear();

		for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
		{
			ModularTraits<typename Ring::Element>::reduce(e1, arr[i * 2][j],
					R._modulus);
			ModularTraits<typename Ring::Element>::reduce(e2, arr[i * 2 + 1][j],
					R._modulus);

			if (!R.isZero(e1) || !R.isZero(e2))
			{
				tmp.IndexData.push_back(j);
				tmp.ValuesData.push_back(e1);
				tmp.ValuesData.push_back(e2);
			}
		}

		if ((float) tmp.size() / (float) DEFAULT_BLOC_WIDTH
				< bloc.get_HYBRID_REPRESENTATION_THRESHOLD())
			bloc[i].swap(tmp);
		else
		{
			Index idx = 0;
			for (uint32 j = 0; j < DEFAULT_BLOC_WIDTH; ++j)
				if (idx < tmp.size() && tmp.IndexData[idx] == j)
				{
					bloc[i].ValuesData.push_back(tmp.at_unchecked(0, idx));
					bloc[i].ValuesData.push_back(tmp.at_unchecked(1, idx));
					idx++;
				}
				else
				{
					bloc[i].ValuesData.push_back(0);
					bloc[i].ValuesData.push_back(0);
				}
		}
	}
}




template <typename Ring>
static inline void reduceDenseArrayModulo(const Ring& R, uint64* arr)
{
	typename Ring::Element e;

	for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
	{
		ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
		arr[j] = (uint64)e;
	}
}

static void DenseScalMulSub__one_row__array_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint64 *arr_source,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));

static void DenseScalMulSub__one_row__array_array(const uint32 av1_col1,
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

static void DenseScalMulSub__two_rows__array_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const uint64 *arr_source1,
		const uint64 *arr_source2,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));

static void DenseScalMulSub__two_rows__array_array(const uint32 av1_col1,
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
	/*
//	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; ++i)
//	{
//		arr1[i] += av1_col1 * arr_source1[i] + av1_col2 * arr_source2[i];
//		arr2[i] += av2_col1 * arr_source1[i] + av2_col2 * arr_source2[i];
//	}

	register uint32 v1__, v2__;

	if(av1_col1 == 0)
	{
		if(av1_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr2[i] += av2_col1 * arr_source1[i] + av2_col2 * arr_source2[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v1__ *= av2_col1;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;
					v2__ *= av2_col2;

					arr2[i+t] += v1__;
					arr2[i+t] += v2__;
				}
			}
		}
		else if(av2_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col2 * arr_source2[i];
				//arr2[i] += av2_col1 * arr_source1[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v1__ *= av2_col1;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;
					v2__ *= av1_col2;

					arr1[i+t] += v2__;
					arr2[i+t] += v1__;
				}
			}
		}
		else	//both != 0
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col2 * arr_source2[i];
				//arr2[i] += av2_col1 * arr_source1[i] + av2_col2 * arr_source2[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					arr2[i+t] += v1__ * av2_col1;
					arr2[i+t] += v2__ * av2_col2;

					arr1[i+t] += v2__ * av1_col2;
				}
			}
		}
	}
	else if(av2_col1 == 0)
	{
		if(av1_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i];
				//arr2[i] += av2_col2 * arr_source2[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					v1__ *= av1_col1;
					v2__ *= av2_col2;

					arr1[i+t] += v1__;
					arr2[i+t] += v2__;
				}
			}
		}
		else if(av2_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i] + av1_col2 * arr_source2[i];


				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					v1__ *= av1_col1;
					v2__ *= av1_col2;

					arr1[i+t] += v1__;
					arr1[i+t] += v2__;
				}
			}
		}
		else	//both != 0
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i] + av1_col2 * arr_source2[i];
				//arr2[i] += av2_col2 * arr_source2[i];

				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					v1__ *= av1_col1;
					arr1[i+t] += v1__;
					arr1[i+t] += av1_col2 * v2__;

					v2__ *= av2_col2;
					arr2[i+t] += v2__;
				}
			}
		}
	}
	else	//both != 0
	{
		if(av1_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i];
				//arr2[i] += av2_col1 * arr_source1[i] + av2_col2 * arr_source2[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					arr1[i+t] += v1__ * av1_col1;

					v1__ *= av2_col1;
					v2__ *= av2_col2;
					arr2[i+t] += v1__;
					arr2[i+t] += v2__;
				}
			}
		}
		else if(av2_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i] + av1_col2 * arr_source2[i];
				//arr2[i] += av2_col1 * arr_source1[i];
				#pragma loop unroll
				for(uint32 t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = arr_source1[i+t] & 0x000000000000ffff;
					v2__ = arr_source2[i+t] & 0x000000000000ffff;

					arr2[i+t] += v1__ * av2_col1;

					v1__ *= av1_col1;
					v2__ *= av1_col2;
					arr1[i+t] += v1__;
					arr1[i+t] += v2__;
				}
			}
		}
		else	//both != 0
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				//arr1[i] += av1_col1 * arr_source1[i] + av1_col2 * arr_source2[i];
				//arr2[i] += av2_col1 * arr_source1[i] + av2_col2 * arr_source2[i];
				#pragma loop unroll
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
	}*/

}


template <typename Index>
static void SparseScalMulSub__one_row__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));
template <typename Index>
static void SparseScalMulSub__one_row__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2)
{
	const int N = v.size();
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
void DenseScalMulSub__one_row__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint32 line,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));

template <typename Index>
void DenseScalMulSub__one_row__vect_array(const uint32 av1_col1,
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
static void DenseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));

template <typename Index>
static void DenseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
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
		DenseScalMulSub__one_row__vect_array(
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
	//const uint16 *p_val2 = p_val + 1;
	
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
	

//	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; ++i)
//	{
//		val1 = v.at_unchecked(0, i);
//		val2 = v.at_unchecked(1, i);
//
//		arr1[i] += (uint64)((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2);
//		arr2[i] += (uint64)((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2);
//	}

	/*if(av1_col1 == 0)
	{
		if(av1_col2 == 0)
		{
//			arr2[i] += (uint64)((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{

				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					v1__ *= av2_col1;
					v2__ *= av2_col2;

					arr2[i+t] += v1__;
					arr2[i+t] += v2__;
				}
			}
		}
		else if(av2_col2 == 0)
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					v1__ *= av2_col1;
					v2__ *= av1_col2;

					arr1[i+t] += v2__;
					arr2[i+t] += v1__;
				}

				//arr1[i] += (uint64)((uint32) av1_col2 * val2);
				//arr2[i] += (uint64)((uint32) av2_col1 * val1);
			}
		}
		else	//both != 0
		{
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					arr1[i+t] += v2__ * av1_col2;

					v1__ *= av2_col1;
					v2__ *= av2_col2;
					arr2[i+t] += v1__;
					arr2[i+t] += v2__;
				}

//				arr1[i] += (uint64)((uint32) av1_col2 * val2);
//				arr2[i] += (uint64)((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2);
			}
		}
	}
	else if(av2_col1 == 0)
	{
		if(av1_col2 == 0)
		{
//			arr1[i] += (uint64)((uint32) av1_col1 * val1);
//			arr2[i] += (uint64)((uint32) av2_col2 * val2);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					v1__ *= av1_col1;
					v2__ *= av2_col2;

					arr1[i+t] += v1__;
					arr2[i+t] += v2__;
				}


			}
		}
		else if(av2_col2 == 0)
		{
//			arr1[i] += (uint64)((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					v1__ *= av1_col1;
					v2__ *= av1_col2;

					arr1[i+t] += v1__;
					arr1[i+t] += v2__;
				}
			}
		}
		else	//both != 0
		{
//			arr1[i] += (uint64)((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2);
//			arr2[i] += (uint64)((uint32) av2_col2 * val2);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					arr2[i+t] += v2__ * av2_col2;

					v1__ *= av1_col1;
					v2__ *= av1_col2;
					arr1[i+t] += v1__;
					arr1[i+t] += v2__;
				}
			}
		}
	}
	else	//both != 0
	{
		if(av1_col2 == 0)
		{
//			arr1[i] += (uint64)((uint32) av1_col1 * val1);
//			arr2[i] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					arr1[i+t] += av1_col1 * v1__;

					v1__ *= av2_col1;
					v2__ *= av2_col2;
					arr2[i+t] += v1__;
					arr2[i+t] += v2__;

				}
			}
		}
		else if(av2_col2 == 0)
		{
//			arr1[i] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
//			arr2[i] += (uint64)((uint32) av2_col1 * val1);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					arr2[i+t] += av2_col1 * v1__;

					v1__ *= av1_col1;
					v2__ *= av1_col2;
					arr1[i+t] += v1__;
					arr1[i+t] += v2__;
				}
			}
		}
		else	//both != 0
		{
//			arr1[i] += (uint64)((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2);
//			arr2[i] += (uint64)((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2);
			for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; i+=UNROLL_STEP__64)
			{
				for(Index t=0; t<UNROLL_STEP__64; t++)
				{
					v1__ = v.at_unchecked(0, i+t);
					v2__ = v.at_unchecked(1, i+t);

					arr1[i+t] += v1__ * av1_col1;
					arr1[i+t] += v2__ * av1_col2;

					v1__ *= av2_col1;
					v2__ *= av2_col2;
					arr2[i+t] += v1__;
					arr2[i+t] += v2__;
				}
			}
		}
	}*/
}


template <typename Index>
static void SparseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
		const uint32 av2_col1,
		const uint32 av1_col2,
		const uint32 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2) __attribute__((noinline));
//AXPY2
template <typename Index>
static void SparseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
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


TIMER_DECLARE_(reduceDenseArrayModulo);
TIMER_DECLARE_(DenseScalMulSub__one_row__array_array);
TIMER_DECLARE_(DenseScalMulSub__two_rows__array_array);
TIMER_DECLARE_(SparseScalMulSub__one_row__vect_array);
TIMER_DECLARE_(DenseScalMulSub__one_row__vect_array);
TIMER_DECLARE_(DenseScalMulSub__two_rows__vect_array);
TIMER_DECLARE_(SparseScalMulSub__two_rows__vect_array);

TIMER_DECLARE_(reduceBlocByRectangularBlocOuter);


TIMER_DECLARE_(T2);
TIMER_DECLARE_(T3);



template <typename Index>
static void reduceBlocByRectangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		const SparseMultilineBloc<uint16, Index>& bloc_B,
		uint64 **Bloc_acc) __attribute__((noinline));

template <typename Index>
static void reduceBlocByRectangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		const SparseMultilineBloc<uint16, Index>& bloc_B,
		uint64 **Bloc_acc)
{
	TIMER_START_(reduceBlocByRectangularBlocOuter);
	typedef Modular<uint16> Ring;

	if(bloc_A.empty() || bloc_B.empty())
			return;

	uint8 is_sparse = 0;

	//for all the rows in Bloc_acc (same as for Bloc_A)
	for(int i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{

		if(bloc_A[i].is_sparse ())
			is_sparse = 1;
		else
			is_sparse = 0;

		const Index N = is_sparse == 1 ? bloc_A[i].size() : DEFAULT_BLOC_WIDTH;

		TIMER_START_(T2);
		for (uint32 j = 0; j < N; ++j)
		{
			TIMER_START_(T3);
			const Index Ap1 = is_sparse == 1 ? bloc_A[i].IndexData[j] : j;

			//R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
			register uint32 Av1_col1 = bloc_A[i].at_unchecked(0, j);
			//R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));
			register uint32 Av2_col1 = bloc_A[i].at_unchecked(1, j);

			if (Av1_col1 != 0)
				Av1_col1 = (uint32)R._modulus - Av1_col1;
			if (Av2_col1 != 0)
				Av2_col1 = (uint32)R._modulus - Av2_col1;


			if (((Ap1 %2) == 0) && (j < (uint32)(N - 1)))
			{
				const Index Ap2 = is_sparse == 1 ? bloc_A[i].IndexData[j+1] : j+1;
				if (Ap2 == Ap1 + 1) //axpy 2 ROWS
				{
// 					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
// 					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
					register uint32 Av1_col2 = bloc_A[i].at_unchecked(0, j+1);
					register uint32 Av2_col2 = bloc_A[i].at_unchecked(1, j+1);

					if (Av1_col2 != 0)
						Av1_col2 = (uint32)R._modulus - Av1_col2;
					if (Av2_col2 != 0)
						Av2_col2 = (uint32)R._modulus - Av2_col2;

					++j;

					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
						TIMER_START_(SparseScalMulSub__two_rows__vect_array);
						SparseScalMulSub__two_rows__vect_array(
								Av1_col1,
								Av2_col1,
								Av1_col2,
								Av2_col2,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(SparseScalMulSub__two_rows__vect_array);
					}
					else
					{
						TIMER_START_(DenseScalMulSub__two_rows__vect_array);
						DenseScalMulSub__two_rows__vect_array(
								Av1_col1,
								Av2_col1,
								Av1_col2,
								Av2_col2,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(DenseScalMulSub__two_rows__vect_array);
					}
				}
				else	//axpy ONE ROW
				{
					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
						TIMER_START_(SparseScalMulSub__one_row__vect_array);
						SparseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(SparseScalMulSub__one_row__vect_array);
					}
					else
					{
						TIMER_START_(DenseScalMulSub__one_row__vect_array);
						DenseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(DenseScalMulSub__one_row__vect_array);
					}
				}
			}
			else	//axpy ONE ROW
			{
					if (bloc_B[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse ())
					{
						TIMER_START_(SparseScalMulSub__one_row__vect_array);
							SparseScalMulSub__one_row__vect_array(
									Av1_col1,
									Av2_col1,
									bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2],
									Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(SparseScalMulSub__one_row__vect_array);
					}
					else
					{
						TIMER_START_(DenseScalMulSub__one_row__vect_array);
						DenseScalMulSub__one_row__vect_array(
								Av1_col1,
								Av2_col1,
								bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE,
								Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
						TIMER_STOP_(DenseScalMulSub__one_row__vect_array);
					}

			}
			TIMER_STOP_(T3);
			//////__builtin_prefetch (&(bloc_A[i].IndexData[j+1]), __PREFETCH_READ, __PREFETCH_LOCALITY_LOW);
			//////__builtin_prefetch (&(bloc_A[i].ValuesData[(j+1)*NB_ROWS_PER_MULTILINE]), __PREFETCH_READ, __PREFETCH_LOCALITY_MODERATE);
		}
		TIMER_STOP_(T2);
	}
	TIMER_STOP_(reduceBlocByRectangularBlocOuter);
	//cout << "1 LINE Sparse "<< sp1l << " DENSE " << de1l << endl;
	//cout << "2 LINE Sparse "<< sp2l << " DENSE " << de2l << endl;
}


/**
 * Reduce the rows inside the bloc by themselves
 */
template<typename Index>
static void reduceBlocByTriangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64** Bloc_acc) __attribute__((noinline));

template<typename Index>
static void reduceBlocByTriangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64** Bloc_acc)
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

				if (Av1_col1 != 0)
					R.negin(Av1_col1);
				if (Av2_col1 != 0)
					R.negin(Av2_col1);

				//if(Ap1 >= i*2)
				//	throw std::runtime_error ("Index pointing to an out of range line.");

				if (Ap1 % 2 == 0 && j < last_idx - 1)
				{
					Ap2 = bloc_A[i].IndexData[j+1];
					if (Ap2 == Ap1+1)				//axpy TWO ARRAYS
					{
						R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
						R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));

						if (Av1_col2 != 0)
							R.negin(Av1_col2);
						if (Av2_col2 != 0)
							R.negin(Av2_col2);

						++j;

						TIMER_START_(DenseScalMulSub__two_rows__array_array);
							DenseScalMulSub__two_rows__array_array(
									Av1_col1,
									Av2_col1,
									Av1_col2,
									Av2_col2,
									Bloc_acc[Ap1],
									Bloc_acc[Ap1+1],
									Bloc_acc[i*2],
									Bloc_acc[i*2+1]);
						TIMER_STOP_(DenseScalMulSub__two_rows__array_array);
					}
					else	//axpy ONE ARRAY
					{
						TIMER_START_(DenseScalMulSub__one_row__array_array);
							DenseScalMulSub__one_row__array_array(
									Av1_col1,
									Av2_col1,
									Bloc_acc[Ap1],
									Bloc_acc[i*2],
									Bloc_acc[i*2+1]);
						TIMER_STOP_(DenseScalMulSub__one_row__array_array);
					}
				}
				else	//axpy ONE ARRAY
				{
					TIMER_START_(DenseScalMulSub__one_row__array_array);
						DenseScalMulSub__one_row__array_array(
								Av1_col1,
								Av2_col1,
								Bloc_acc[Ap1],
								Bloc_acc[i*2],
								Bloc_acc[i*2+1]);
					TIMER_STOP_(DenseScalMulSub__one_row__array_array);
				}
			}

			TIMER_START_(reduceDenseArrayModulo);
				reduceDenseArrayModulo(R, Bloc_acc[i*2]);
			TIMER_STOP_(reduceDenseArrayModulo);

			if (bloc_A[i].size() > 1) //reduce lines within the same multiline
			{
				int j=bloc_A[i].size()-2;
				Ap1 = bloc_A[i].IndexData[j];
				R.copy(Av1_col1, bloc_A[i].at_unchecked(1, j));

				if (Av1_col1 != 0)
				{
					R.negin(Av1_col1);

					//TODO: make this one loop
					for (uint32 t = 0; t < DEFAULT_BLOC_WIDTH; ++t)
						Bloc_acc[i*2+1][t] += (uint32) Av1_col1 * Bloc_acc[i*2][t];
				}

			}
			TIMER_START_(reduceDenseArrayModulo);
				reduceDenseArrayModulo(R, Bloc_acc[i*2+1]);
			TIMER_STOP_(reduceDenseArrayModulo);
//			else	//these lines should be the remaning % bloc height! corresponding Bloc_acc should be null too
					//no need to reduce them in the ring
//				cout << "PROBLEM " << endl;
		/*}
		else	//dense
		{
			cout << "EROOR: ROWS IN A SHALL NOT BE DENSE" << endl;
			throw std::logic_error ("EROOR: ROWS IN A SHALL NOT BE DENSE");
		}*/

	}  //for i
}

template<typename Index>
void MatrixOps::reducePivotsByPivots(const Modular<uint16>& R,
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
		//dense_bloc[i] = new uint64[B.bloc_width()];
		//dense_bloc[i] = (uint64 *) memalign(16, B.bloc_width() * sizeof(uint64));
		posix_memalign((void**)&dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseArrayModulo);
	TIMER_RESET_(DenseScalMulSub__one_row__array_array);
	TIMER_RESET_(DenseScalMulSub__two_rows__array_array);
	TIMER_RESET_(SparseScalMulSub__one_row__vect_array);
	TIMER_RESET_(DenseScalMulSub__one_row__vect_array);
	TIMER_RESET_(DenseScalMulSub__two_rows__vect_array);
	TIMER_RESET_(SparseScalMulSub__two_rows__vect_array);

	TIMER_RESET_(reduceBlocByRectangularBlocOuter);


	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) B.coldim() / B.bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) A.rowdim() / A.bloc_height());

	//for all columns of blocs of B
	for (uint32 i = 0; i < nb_column_blocs_B; ++i)
	{
		//report << "Column B " << i << endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 j = 0;
				j < nb_row_blocs_A;
				++j)
		{
			const uint32 first_bloc_idx = A.FirstBlocsColumIndexes[j] / A.bloc_width();
			const uint32 last_bloc_idx = min(A[j].size () - 1, j);

			//report << "\tRow A " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j
			TIMER_START_(memsetBlocToZero);
				memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				copySparseBlocToDenseBlocArray(R, B[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "/" << (uint32) std::ceil((double) B.coldim() / B.bloc_width()) << "\trow\t" << j << std::ends;
#endif
			//for all the blocs in the current row of A
			for (uint32 k = 0; k < last_bloc_idx; ++k)
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
				copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i], false);
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
	TIMER_REPORT_(reduceBlocByRectangularBlocOuter);

	TIMER_REPORT_(T2);
	TIMER_REPORT_(T3);

	TIMER_REPORT_(SparseScalMulSub__one_row__vect_array);
	TIMER_REPORT_(DenseScalMulSub__one_row__vect_array);
	TIMER_REPORT_(DenseScalMulSub__two_rows__vect_array);
	TIMER_REPORT_(SparseScalMulSub__two_rows__vect_array);

	report << endl;
	TIMER_REPORT_(reduceBlocByTriangularBloc);
	TIMER_REPORT_(DenseScalMulSub__one_row__array_array);
	TIMER_REPORT_(DenseScalMulSub__two_rows__array_array);
	TIMER_REPORT_(reduceDenseArrayModulo);


	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}

template<typename Index>
void MatrixOps::reduceNonPivotsByPivots(const Modular<uint16>& R,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
		SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D)
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
		//dense_bloc[i] = (uint64 *) memalign(16, D.bloc_width() * sizeof(uint64));
		posix_memalign((void**)&dense_bloc[i], 16, D.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseArrayModulo);
	TIMER_RESET_(reduceBlocByRectangularBlocOuter);

	TIMER_RESET_(DenseScalMulSub__one_row__array_array);
	TIMER_RESET_(DenseScalMulSub__two_rows__array_array);
	TIMER_RESET_(SparseScalMulSub__one_row__vect_array);
	TIMER_RESET_(DenseScalMulSub__one_row__vect_array);
	TIMER_RESET_(DenseScalMulSub__two_rows__vect_array);
	TIMER_RESET_(SparseScalMulSub__two_rows__vect_array);


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
				memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				copySparseBlocToDenseBlocArray(R, D[j][i], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "/" << (uint32) std::ceil((double) D.coldim() / D.bloc_width()) << "\trow\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B\t" << k + first_bloc_idx << endl;
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
	TIMER_REPORT_(reduceBlocByRectangularBlocOuter);
	TIMER_REPORT_(SparseScalMulSub__one_row__vect_array);
	TIMER_REPORT_(DenseScalMulSub__one_row__vect_array);
	TIMER_REPORT_(DenseScalMulSub__two_rows__vect_array);
	TIMER_REPORT_(SparseScalMulSub__two_rows__vect_array);


	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}

template <typename Index>
static void copyMultilineToDenseArrays64(const MultiLineVector<uint16, Index>& v, uint64 *arr1, uint64 *arr2)
{
	register uint32 idx;
	for (uint32 i = 0; i < v.size(); ++i)
	{
		idx = v.IndexData[i];
		arr1[idx] = (uint64) v.at_unchecked(0, i);
		arr2[idx] = (uint64) v.at_unchecked(1, i);
	}
}


template<typename Index>
void reduceBlocByRectangularBloc_C(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_C,
		const SparseMultilineBloc<uint16, Index>& bloc_A, uint64 **bloc_dense)
{
//	cout << "In reduceBlocByRectangularBloc_C" << endl;
	typedef Modular<uint16> Ring;

	uint8 is_sparse = 0;
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
	{
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

					TIMER_START_(SparseScalMulSub__two_rows__vect_array);
						SparseScalMulSub__two_rows__vect_array(
								Cv1_col1,
								Cv2_col1,
								Cv1_col2,
								Cv2_col2,
								bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
								bloc_dense[i * 2],
								bloc_dense[i * 2 + 1]);
					TIMER_STOP_(SparseScalMulSub__two_rows__vect_array);

					++j;
				}
				else
				{
					TIMER_START_(SparseScalMulSub__one_row__vect_array);
						SparseScalMulSub__one_row__vect_array(
								Cv1_col1,
								Cv2_col1,
								bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
								Cp1 % NB_ROWS_PER_MULTILINE,
								bloc_dense[i * 2],
								bloc_dense[i * 2 + 1]);
					TIMER_STOP_(SparseScalMulSub__one_row__vect_array);
				}
			}
			else
			{
				TIMER_START_(SparseScalMulSub__one_row__vect_array);
					SparseScalMulSub__one_row__vect_array(
							Cv1_col1,
							Cv2_col1,
							bloc_A[Cp1 / NB_ROWS_PER_MULTILINE],
							Cp1 % NB_ROWS_PER_MULTILINE,
							bloc_dense[i * 2],
							bloc_dense[i * 2 + 1]);
				TIMER_STOP_(SparseScalMulSub__one_row__vect_array);
			}
		}
	}

}

template <typename Index>
void reduceBlocByTriangularBloc_C(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64 **bloc_dense)
{
//	cout << "In reduceBlocByTriangularBloc_C" << endl;

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
				cout << "PROBLEM bloc_A TRIANGULAR\n";
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

				TIMER_START_(SparseScalMulSub__two_rows__vect_array);
						SparseScalMulSub__two_rows__vect_array(
								Av1_col2,	//inversed!
								Av2_col2,
								Av1_col1,
								Av2_col1,
								bloc_A[Ap1 / NB_ROWS_PER_MULTILINE],
								bloc_dense[i * 2],
								bloc_dense[i * 2 + 1]);

						bloc_dense[i * 2][Ap1] = Av1_col1;
						bloc_dense[i * 2 + 1][Ap1] = Av2_col1;

						bloc_dense[i * 2][Ap1 - 1] = Av1_col2;
						bloc_dense[i * 2 + 1][Ap1 - 1] = Av2_col2;
				TIMER_STOP_(SparseScalMulSub__two_rows__vect_array);

				--j;
			}
			else
			{
				TIMER_START_(SparseScalMulSub__one_row__vect_array);
					SparseScalMulSub__one_row__vect_array(
							Av1_col1,
							Av2_col1,
							bloc_A[Ap1 / NB_ROWS_PER_MULTILINE],
							Ap1 % NB_ROWS_PER_MULTILINE,
							bloc_dense[i * 2],
							bloc_dense[i * 2 + 1]);

					bloc_dense[i * 2][Ap1] = Av1_col1;
					bloc_dense[i * 2 + 1][Ap1] = Av2_col1;
				TIMER_STOP_(SparseScalMulSub__one_row__vect_array);
			}


		}
	}

}

template<typename Index>
void MatrixOps::reduceC(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "[MatrixOps::reduceC] In spec Modular<uint16> Bloc version" << std::endl;

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, C.bloc_width() * sizeof(uint64));

	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)C.rowdim() / C.bloc_height ());
	const uint32 nb_col_blocs_C = (uint32)std::ceil((double)C.coldim() / C.bloc_width ());

	//for all columns of blocs of D
	for (uint32 i = 0; i < nb_row_blocs_C; ++i)
	{
		check_equal_or_raise_exception(C.FirstBlocsColumIndexes[i], 0);

		//report << "ROW C\t" << i << endl;

		//for all rows of blocs in C
		for (int j = nb_col_blocs_C-1; j >= 0; --j)
		{
			//report << "\tColumn C|A\t" << j << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			TIMER_START_(memsetBlocToZero);
				memsetToZero(dense_bloc);
			TIMER_STOP_(memsetBlocToZero);

			TIMER_START_(copySparseBlocToDenseBlocArray);
				copySparseBlocToDenseBlocArray(R, C[i][j], dense_bloc);
			TIMER_STOP_(copySparseBlocToDenseBlocArray);

#ifdef SHOW_PROGRESS
		report << "                                                                                    \r";
		report << "\trow\t" << i << "/" << nb_row_blocs_C << "\tcolumn\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (int k = nb_col_blocs_C-1; k>j; --k)
			{
				//report << "\t\tBloc in C and A\t" << k << endl;
				TIMER_START_(reduceBlocByRectangularBloc);
					reduceBlocByRectangularBloc_C(R, C[i][k], A[k][j], dense_bloc);
				TIMER_STOP_(reduceBlocByRectangularBloc);
			}

			TIMER_START_(reduceBlocByTriangularBloc);
				reduceBlocByTriangularBloc_C(R, A[j][j], dense_bloc);
			TIMER_STOP_(reduceBlocByTriangularBloc);

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				copyDenseBlocArrayToSparseBloc(R, dense_bloc, C[i][j]);
			TIMER_STOP_(copyDenseBlocArrayToSparseBloc);

		}
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

}



template<typename Index>
uint32 MatrixOps::echelonize(const Modular<uint16>& R, const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A)
{

//	uint16 *P[A.rowdim ()]  __attribute__((aligned(0x1000)));
//	for (uint32 i = 0; i < A.rowdim (); ++i)
//		posix_memalign((void**)&P[i], 16, A.rowdim () * sizeof(uint64));
//
//
//	const uint32 nb_row_blocs_A = (uint32)std::ceil((double)A.rowdim() / A.bloc_height ());
//	const uint32 nb_col_blocs_A = (uint32)std::ceil((double)A.coldim() / A.bloc_width ());
//
//	//for all columns of blocs of D
//	for (uint32 i = 0; i < nb_col_blocs_A; ++i)
//	{
//		//report << "COL C\t" << i << endl;
//
//		//for all rows of blocs in C
//		for (uint32 j = 0; j < nb_row_blocs_A; ++j)
//		{
//
//		}
//	}
}











































template<typename Index>
uint32 MatrixOps::echelonize(const Modular<uint16>& R, const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
		SparseMultilineMatrix<uint16>& outMatrix, bool destruct_in_matrix)
{
	typedef uint16 Element;
	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;
	outMatrix = SparseMultilineMatrix<uint16> (inMatrix.rowdim (), inMatrix.coldim ());

	//copy bloc matrix to multiline matrix
	check_equal_or_raise_exception(inMatrix.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(inMatrix.isFilledWithEmptyBlocs (), true);

	uint32 curr_row_base = 0;

	//this will copy the bloc matrix to a multiline matrix by inversing the order of the rows
	for(uint32 i=0; i<inMatrix.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
// 		cout << "i " << i << endl;
		if(inMatrix[i].empty ())
		{
			curr_row_base += inMatrix.bloc_height ();
			continue;
		}

		for (uint32 j = 0; j < inMatrix[i].size(); ++j) //for each bloc in the row
		{
// 			cout << "\tj " << j << endl;
			if(inMatrix[i][j].empty ())
				continue;

			uint32 bloc_idx = inMatrix.FirstBlocsColumIndexes[i] + (inMatrix.bloc_width () * j);

			uint32 idx;
			Element val1, val2;
			for (uint16 k = 0;
					k < inMatrix[i][j].bloc_height(); ++k) //for each row in the bloc
			{
// 				cout << "\t\tk " << k << endl;
				if(inMatrix[i][j][k].empty ())
					continue;

				if(inMatrix[i][j][k].is_sparse ())
				{
// 					cout << "sparse size " << inMatrix[i][j][k].size() << endl;
					for (uint32 p = 0; p < inMatrix[i][j][k].size (); ++p)
					{
						idx = inMatrix[i][j][k].IndexData[p];
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);
// 						cout <<"ok rowdim " << outMatrix.rowdim() << " row " << curr_row_base+k << endl;
						outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + idx);
						outMatrix[curr_row_base+k].ValuesData.push_back(val1);
						outMatrix[curr_row_base+k].ValuesData.push_back(val2);
					}
// 					cout << "sparse end" << endl;
				}
				else
				{
// 					cout << "dense" << endl;
					for (uint32 p = 0; p < inMatrix.bloc_width (); ++p)
					{
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						if(val1 != 0 && val2 != 0)
						{
							outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + p);
							outMatrix[curr_row_base+k].ValuesData.push_back(val1);
							outMatrix[curr_row_base+k].ValuesData.push_back(val2);
						}
					}
// 					cout << "dense end" << endl;
				}
			}

		}

		curr_row_base += inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;
	}

	return NS::echelonize(R, outMatrix);
	//return 0;
}



#endif /* MATRIX_OPS_H_ */




