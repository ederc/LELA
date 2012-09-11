/*
 * matrix-op.C
 *
 *  Created on: 19 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_OP_MULTILINE_
#define MATRIX_OP_MULTILINE_

#include <vector>
#include <pthread.h>
#include <assert.h>

#include "types.h"
//#include "matrix-util-m.C"
//#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

#define UNROLL_STEP 32

#ifdef DETAILED_PROFILE_TIMERS
#  define TIMER_DECLARE_(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START_(part) part##_timer.start ()
#  define TIMER_STOP_(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT_(part) \
		commentator.report (Commentator::LEVEL_ALWAYS, TIMING_MEASURE)		\
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
//#  define TIMER_REPORT_(part)
#else
#  define TIMER_DECLARE_(part)
#  define TIMER_START_(part)
#  define TIMER_STOP_(part)
#  define TIMER_REPORT_(part)
#endif //DETAILED_PROFILE_TIMERS

using namespace LELA;
using namespace std;

#define START_UNROLL_CODE				\
	while(x<xl)					\
	{						\
		idx = v.IndexData[x];			\
		val1 = v.at_unchecked(0, x);		\
		val2 = v.at_unchecked(1, x);

#define MIDDLE_UNROLL_CODE				\
	++x;						\
	}						\
							\
	for(x=xl; x<sz; x+=UNROLL_STEP)			\
	{

#define MIDDLE_UNROLL_CODE2				\
		for(uint8 t=0; t<UNROLL_STEP; ++t)		\
		{					\
			idx = v.IndexData[x+t];		\
							\
			val1 = v.at_unchecked(0, x+t);	\
			val2 = v.at_unchecked(1, x+t);

#define END_UNROLL_CODE					\
		}					\
	}

namespace NS{

static inline void razArray64(uint64 arr[], uint32 arrSize)
{
	memset(arr, 0, arrSize*sizeof(uint64));
}

template <typename Index>
static void copyMultilineToDenseArrays64(const MultiLineVector<uint16, Index>& v,
		uint64 *arr1, uint64 *arr2)
{
	register uint32 idx;
	/*for (uint32 i = 0; i < v.size(); ++i)
	{
		idx = v.IndexData[i];
		arr1[idx] = (uint64) v.at_unchecked(0, i);
		arr2[idx] = (uint64) v.at_unchecked(1, i);
	}*/

	uint8 xl = v.size() % UNROLL_STEP;
	uint32 sz = v.size();
	register uint32 x = 0;
	while (x < xl)
	{
		idx = v.IndexData[x];
		arr1[idx] = (uint64) v.at_unchecked(0, x);
		arr2[idx] = (uint64) v.at_unchecked(1, x);

		++x;
	}

	for (x = xl; x < sz; x += UNROLL_STEP)
	{

#pragma loop unroll
		for (uint8 t = 0; t < UNROLL_STEP; ++t)
		{
			idx = v.IndexData[x+t];
			arr1[idx] = (uint64) v.at_unchecked(0, x + t);
			arr2[idx] = (uint64) v.at_unchecked(1, x + t);
		}
	}
}


template<typename Ring, typename Index>
static void copyDenseArraysToMultilineVector64(const Ring& R,
		const uint64 *arr1,
		const uint64 *arr2,
		const uint32 size,
		MultiLineVector<typename Ring::Element, Index>& v,
		bool reduce)
{
	MultiLineVector<typename Ring::Element, Index> tmp(v.nb_lines());
	typename Ring::Element e1, e2;


	uint32 sz=0;

	if (reduce)
	{
		for (uint32 i = 0; i < size; ++i)
			if ((arr1[i] % R._modulus != 0) || (arr2[i] % R._modulus != 0))
				sz++;

		tmp.reserve(sz);

		for (uint32 i = 0; i < size; ++i)
		{
			if ((arr1[i] % R._modulus != 0) || (arr2[i] % R._modulus != 0))
			{
				ModularTraits<typename Ring::Element>::reduce(e1, arr1[i], R._modulus);
				ModularTraits<typename Ring::Element>::reduce(e2, arr2[i], R._modulus);

				tmp.IndexData.push_back(i);
				tmp.ValuesData.push_back(e1);
				tmp.ValuesData.push_back(e2);
			}
		}
	}
	else
	{
		for (uint32 i = 0; i < size; ++i)
			if ((arr1[i] != 0) || (arr2[i] != 0))
				sz++;

		tmp.reserve(sz);

		for (uint32 i = 0; i < size; ++i)
		{
			if ((arr1[i] != 0) || (arr2[i] != 0))
			{
				tmp.IndexData.push_back(i);
				tmp.ValuesData.push_back((uint16) arr1[i]);
				tmp.ValuesData.push_back((uint16) arr2[i]);
			}
		}
	}
	v.swap(tmp);
}

template <typename Index>
inline void axpy(const uint16 av1_col1,
		const uint16 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint16 line,
		uint64 *arr1,
		uint64 *arr2,
		uint16 start_from = 0)
{
	lela_check(line < v.nb_lines());
	const uint32 sz = v.size();

	/*for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		val1 = v.at_unchecked(line, i);

		arr1[idx] += (uint32) av1_col1 * val1;
		arr2[idx] += (uint32) av2_col1 * val1;
	}*/


	const uint8 xl = v.size() % UNROLL_STEP;
	register uint32 x = start_from;

	register uint32 idx;
	register uint16 val1;

	if(av1_col1 == 0 && av2_col1 == 0)
		return;

	if (av1_col1 != 0 && av2_col1 != 0) //cannot both be 0
	{
		while (x < xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr1[idx] += (uint32) av1_col1 * val1;
			arr2[idx] += (uint32) av2_col1 * val1;

			++x;
		}

		for (x = xl; x < sz; x += UNROLL_STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < UNROLL_STEP; ++t)
			{
				idx = v.IndexData[x + t];
				val1 = v.at_unchecked(line, x + t);
				arr1[idx] += (uint32) av1_col1 * val1;
				arr2[idx] += (uint32) av2_col1 * val1;
			}
		}
	}
	else if (av1_col1 != 0)
	{
		while (x < xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr1[idx] += (uint32) av1_col1 * val1;

			++x;
		}

		for (x = xl; x < sz; x += UNROLL_STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < UNROLL_STEP; ++t)
			{
				idx = v.IndexData[x + t];
				val1 = v.at_unchecked(line, x + t);
				arr1[idx] += (uint32) av1_col1 * val1;
			}
		}
	}
	else //av2_col1 != 0
	{
		while (x < xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr2[idx] += (uint32) av2_col1 * val1;

			++x;
		}

		for (x = xl; x < sz; x += UNROLL_STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < UNROLL_STEP; ++t)
			{
				idx = v.IndexData[x + t];
				val1 = v.at_unchecked(line, x + t);
				arr2[idx] += (uint32) av2_col1 * val1;
			}
		}
	}

}

template <typename Index>
inline void axpy2(const uint16 av1_col1,
		const uint16 av2_col1,
		const uint16 av1_col2,
		const uint16 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2,
		uint16 start_from = 0)
{
	register uint32 idx;
	register uint16 val1, val2;

	/*for (uint32 i = 0; i < v.size(); ++i)
	{
		idx = v.IndexData[i];

		val1 = v.at_unchecked(0, i);
		val2 = v.at_unchecked(1, i);

		arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
		arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
	}

	return;*/
	/*			START_UNROLL_CODE
				arr1[idx] += (uint32) av1_col1 * val1;
				arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint32) av2_col1 * val1;
				arr2[idx] += (uint32) av2_col2 * val2;
				MIDDLE_UNROLL_CODE
#pragma loop unroll
			MIDDLE_UNROLL_CODE2
					arr1[idx] += (uint32) av1_col1 * val1;
					arr1[idx] += (uint32) av1_col2 * val2;

					arr2[idx] += (uint32) av2_col1 * val1;
					arr2[idx] += (uint32) av2_col2 * val2;
				END_UNROLL_CODE
*/
	if(av1_col1 == 0 && av2_col1 == 0)
	{
		axpy(av1_col2, av2_col2, v, 1, arr1, arr2, start_from);
		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		axpy(av1_col1, av2_col1, v, 0, arr1, arr2, start_from);
		return;
	}

	if (av1_col1 == 0)
	{
		if (av1_col2 == 0)
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				//arr1[idx] += (uint32) av1_col1 * val1;
				//arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
		else if (av2_col2 == 0)
		{
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				//arr1[idx] += (uint32) av1_col1 * val1;
				arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint32) av2_col1 * val1;
				//arr2[idx] += (uint32) av2_col2 * val2;
			}
		}
		else
		{
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				//arr1[idx] += (uint32) av1_col1 * val1;
				arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
		}
	}
	else if (av2_col1 == 0)
	{
		if (av1_col2 == 0)
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint32) av1_col1 * val1;
				//arr1[idx] += (uint32) av1_col2 * val2;

				//arr2[idx] += (uint32) av2_col1 * val1;
				arr2[idx] += (uint32) av2_col2 * val2;
			}
		else if (av2_col2 == 0)
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));

				//arr2[idx] += (uint32) av2_col1 * val1;
				//arr2[idx] += (uint32) av2_col2 * val2;
			}
		else
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));

				//arr2[idx] += (uint32) av2_col1 * val1;
				arr2[idx] += (uint32) av2_col2 * val2;
			}
	}
	else // av1_col1 && av2_col1 != 0
	{
		if (av1_col2 == 0)
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint32) av1_col1 * val1;
				//arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
		else if (av2_col2 == 0)
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));

				arr2[idx] += (uint32) av2_col1 * val1;
				//arr2[idx] += (uint32) av2_col2 * val2;
			}
		else
			for (uint32 i = start_from; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
	}
}

template <typename Index>
inline void axpy2_unroll(const uint16 av1_col1,
		const uint16 av2_col1,
		const uint16 av1_col2,
		const uint16 av2_col2,
		const MultiLineVector<uint16, Index>& v,
		uint64 *arr1,
		uint64 *arr2,
		uint16 start_from = 0)
{

	if(av1_col1 == 0 && av2_col1 == 0)
	{
		axpy(av1_col2, av2_col2, v, 1, arr1, arr2, start_from);
		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		axpy(av1_col1, av2_col1, v, 0, arr1, arr2, start_from);
		return;
	}

	const uint32 sz = v.size();
	const uint8 xl = v.size() % UNROLL_STEP;
	uint32 x = start_from;

	register uint32 idx;
	register uint16 val1, val2;

	if (av1_col1 == 0)
	{
		if (av1_col2 == 0)
		{
			START_UNROLL_CODE
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			END_UNROLL_CODE
		}
		else if (av2_col2 == 0)
		{
			START_UNROLL_CODE
				arr1[idx] += (uint32) av1_col2 * val2;
				arr2[idx] += (uint32) av2_col1 * val1;
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint32) av1_col2 * val2;
				arr2[idx] += (uint32) av2_col1 * val1;
			END_UNROLL_CODE
		}
		else
		{
			START_UNROLL_CODE
				arr1[idx] += (uint32) av1_col2 * val2;
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint32) av1_col2 * val2;
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			END_UNROLL_CODE
		}
	}
	else if (av2_col1 == 0)
	{
		if (av1_col2 == 0)
		{
			START_UNROLL_CODE
				arr1[idx] += (uint32) av1_col1 * val1;
				arr2[idx] += (uint32) av2_col2 * val2;
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint32) av1_col1 * val1;
				arr2[idx] += (uint32) av2_col2 * val2;
			END_UNROLL_CODE
		}
		else if (av2_col2 == 0)
		{
			START_UNROLL_CODE
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
			END_UNROLL_CODE
		}
		else
		{
			START_UNROLL_CODE
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint32) av2_col2 * val2;
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint32) av2_col2 * val2;
			END_UNROLL_CODE
		}
	}
	else // av1_col1 && av2_col1 != 0
	{
		if (av1_col2 == 0)
		{
			START_UNROLL_CODE
				arr1[idx] += (uint32) av1_col1 * val1;
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint32) av1_col1 * val1;
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			END_UNROLL_CODE
		}
		else if (av2_col2 == 0)
		{
			START_UNROLL_CODE
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint32) av2_col1 * val1;
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint32) av2_col1 * val1;
			END_UNROLL_CODE
		}
		else
		{
			START_UNROLL_CODE
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			MIDDLE_UNROLL_CODE
			#pragma loop unroll
			MIDDLE_UNROLL_CODE2
				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			END_UNROLL_CODE
		}
	}
}



template<typename Ring>
void reducePivotsByPivots(const Ring& R, const SparseMultilineMatrix<uint16>& A,
							SparseMultilineMatrix<uint16>& B)
{
	lela_check(A.rowdim () == B.rowdim ());
	lela_check(A.rowdim () == A.coldim ());

	//typedef Modular<uint16> Ring;
	typedef SparseMultilineMatrix<uint16> Matrix;

	typename Matrix::ConstRowIterator i_A;
	typename Matrix::RowIterator i_B;

	uint32 B_coldim = B.coldim();
	uint64 *tmpDenseArray1; // __attribute__((aligned(64)));
		tmpDenseArray1 = new uint64[B_coldim];

	uint64 *tmpDenseArray2;// __attribute__((aligned(64)));
		tmpDenseArray2 = new uint64[B_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1Timer);
	TIMER_DECLARE_(Axpy2Timer);
	TIMER_DECLARE_(AxpyOuterTimer);
	TIMER_DECLARE_(AxpyLeftTimer);

#ifdef SHOW_PROGRESS
	uint32 i = A.rowdim () / A.nb_lines_per_bloc ();
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	i_A = A.rowEnd();
	i_B = B.rowEnd();

	while (i_A != A.rowBegin())
	{ //for each multiline
		--i_A;
		--i_B;

		TIMER_START_(RazArrayTimer);
			razArray64(tmpDenseArray1, B_coldim);
			razArray64(tmpDenseArray2, B_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(*i_B, tmpDenseArray1, tmpDenseArray2);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		//line = A.nb_lines_per_bloc() - 1;		//last line
		//do {		//for each line
#ifdef SHOW_PROGRESS
		--i;
		report << "                                                                    \r";
		report << "\t" << i << std::ends;
#endif

		typename Ring::Element Av1_col1, Av2_col1;
		uint32 Ap1;
		typename Ring::Element Av1_col2, Av2_col2;
		uint32 Ap2;


		TIMER_START_(AxpyOuterTimer);
		for (uint32 j = 2; j < i_A->size(); ++j) //skip first two elements
		{
			Ap1 = i_A->IndexData[j];
			R.copy(Av1_col1, i_A->at_unchecked(0, j));
			R.copy(Av2_col1, i_A->at_unchecked(1, j));

			if (Av1_col1 != 0)
				R.negin(Av1_col1);
			if (Av2_col1 != 0)
				R.negin(Av2_col1);

			if (Ap1 % 2 == 0 && j < i_A->size() - 1)
			{
				assert(j != i_A->size() - 1);
				Ap2 = i_A->IndexData[j+1];
				if (Ap2 == Ap1+1)				//have two consecutive lines
				{
					R.copy(Av1_col2, i_A->at_unchecked(0, j+1));
					R.copy(Av2_col2, i_A->at_unchecked(1, j+1));

					if (Av1_col2 != 0)
						R.negin(Av1_col2);
					if (Av2_col2 != 0)
						R.negin(Av2_col2);

					++j;

					TIMER_START_(Axpy2Timer);

					axpy2_unroll(Av1_col1, Av2_col1,
							Av1_col2,
							Av2_col2,
							B[Ap1 / A.nb_lines_per_bloc()],
							tmpDenseArray1,
							tmpDenseArray2);
					TIMER_STOP_(Axpy2Timer);
				}
				else
				{
					TIMER_START_(Axpy1Timer);

					axpy(Av1_col1, Av2_col1,
							B[Ap1 / A.nb_lines_per_bloc()],
							Ap1 % A.nb_lines_per_bloc(),
							tmpDenseArray1,
							tmpDenseArray2);
					TIMER_STOP_(Axpy1Timer);
				}
			}
			else
			{
				TIMER_START_(Axpy1Timer);

				axpy(Av1_col1, Av2_col1,
						B[Ap1 / A.nb_lines_per_bloc()],
						Ap1 % A.nb_lines_per_bloc(),
						tmpDenseArray1,
						tmpDenseArray2);
				TIMER_STOP_(Axpy1Timer);
			}
		}

		TIMER_STOP_(AxpyOuterTimer);

		TIMER_START_(AxpyLeftTimer);
		bool reduce = true;

		if (i_A->size() > 1) //reduce lines within the same multiline
		{
			//j=1
			Ap1 = i_A->IndexData[1];
			R.copy(Av1_col1, i_A->at_unchecked(0, 1));

			if (Av1_col1 != 0)
			{
				R.negin(Av1_col1);

				for (uint32 t = 0; t < B_coldim; ++t)
					tmpDenseArray2[t] %= R._modulus; //Make sure product in next loop doesn't overflow

				for (uint32 t = 0; t < B_coldim; ++t)
					tmpDenseArray1[t] += (uint32) Av1_col1 * tmpDenseArray2[t];

				for (uint32 t = 0; t < B_coldim; ++t)
					tmpDenseArray1[t] %= R._modulus;

				reduce = false;
			}
		} TIMER_STOP_(AxpyLeftTimer);

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArraysToMultilineVector64(R, tmpDenseArray1, tmpDenseArray2, B_coldim, *i_B, reduce);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		//} while(line != 0);		//end for each line
	} //for each multiline

#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
	TIMER_REPORT_(Axpy1Timer);
	TIMER_REPORT_(Axpy2Timer);
	TIMER_REPORT_(AxpyOuterTimer);
	TIMER_REPORT_(AxpyLeftTimer);

	//report << "axpy " << called_axpy << "   axpy2 " << called_axpy2 << std::endl;
}

template<typename Ring>
void reduceNonPivotsByPivots(const Ring R,
							 const SparseMultilineMatrix<uint16>& C,
							 const SparseMultilineMatrix<uint16>& B,
							 SparseMultilineMatrix<uint16>& D,
							 bool invert_scalars = true)
{
	lela_check(B.coldim () == D.coldim ());
	lela_check(C.rowdim () == D.rowdim ());

#ifdef SHOW_PROGRESS
	uint32 i=0;
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	typedef SparseMultilineMatrix<uint16> Matrix;

	typename Matrix::ConstRowIterator i_C;
	typename Matrix::RowIterator i_D;

	uint32 D_coldim = D.coldim();
	uint64 tmpDenseArray1[D_coldim];
	uint64 tmpDenseArray2[D_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1Timer);
	TIMER_DECLARE_(Axpy2Timer);
	TIMER_DECLARE_(AxpyOuterTimer);
	TIMER_DECLARE_(AxpyLeftTimer);

	int x=0;
	for(i_C = C.rowBegin (), i_D = D.rowBegin (); i_C != C.rowEnd (); ++i_C, ++i_D){

#ifdef SHOW_PROGRESS
		++i;
		report << "                                                                    \r";
		report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			razArray64(tmpDenseArray1, D_coldim);
			razArray64(tmpDenseArray2, D_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(D[x], tmpDenseArray1, tmpDenseArray2);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		typename Ring::Element Cv1_col1, Cv2_col1;
		uint32 Cp1;
		typename Ring::Element Cv1_col2, Cv2_col2;
		uint32 Cp2;

		TIMER_START_(AxpyOuterTimer);
		for (uint32 j = 0; j < i_C->size(); ++j)
		{
			Cp1 = i_C->IndexData[j];
			R.copy(Cv1_col1, i_C->at_unchecked(0, j));
			R.copy(Cv2_col1, i_C->at_unchecked(1, j));

			if(invert_scalars)	//when the new method is used, scalars of matrix C are already 
						//in their final form, no need to compute the opposite
			{
				if (Cv1_col1 != 0)
					R.negin(Cv1_col1);
				if (Cv2_col1 != 0)
					R.negin(Cv2_col1);
			}
			if (Cp1 % 2 == 0 && j < i_C->size() - 1)
			{
				assert(j != i_C->size() - 1);
				Cp2 = i_C->IndexData[j + 1];
				if (Cp2 == Cp1 + 1) //have two consecutive lines
				{
					R.copy(Cv1_col2, i_C->at_unchecked(0, j + 1));
					R.copy(Cv2_col2, i_C->at_unchecked(1, j + 1));

					if(invert_scalars)
					{
						if (Cv1_col2 != 0)
							R.negin(Cv1_col2);
						if (Cv2_col2 != 0)
							R.negin(Cv2_col2);
					}
					++j;

					TIMER_START_(Axpy2Timer);
					axpy2_unroll(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
							B[Cp1 / C.nb_lines_per_bloc()],
							tmpDenseArray1,
							tmpDenseArray2);
					TIMER_STOP_(Axpy2Timer);
				}
				else
				{
					TIMER_START_(Axpy1Timer);
						axpy(Cv1_col1, Cv2_col1, B[Cp1 / C.nb_lines_per_bloc()],
							 Cp1 % C.nb_lines_per_bloc(),
							 tmpDenseArray1,
							 tmpDenseArray2);
					TIMER_STOP_(Axpy1Timer);
				}
			}
			else
			{
				TIMER_START_(Axpy1Timer);
					axpy(Cv1_col1, Cv2_col1, B[Cp1 / C.nb_lines_per_bloc()],
						 Cp1 % C.nb_lines_per_bloc(),
						 tmpDenseArray1,
						 tmpDenseArray2);
				TIMER_STOP_(Axpy1Timer);
			}
		}
		TIMER_STOP_(AxpyOuterTimer);

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArraysToMultilineVector64(R, tmpDenseArray1, tmpDenseArray2,
				D_coldim, D[x], true);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		x++;
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif
		TIMER_REPORT_(RazArrayTimer);
		TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
		TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
		TIMER_REPORT_(Axpy1Timer);
		TIMER_REPORT_(Axpy2Timer);
		TIMER_REPORT_(AxpyOuterTimer);
		TIMER_REPORT_(AxpyLeftTimer);
}


//performs a pseudo reduction of C by A (A^-1 transpose(C))
template<typename Ring>
void reduceCD_onlyC(const Ring& R, const SparseMultilineMatrix<uint16>& A, SparseMultilineMatrix<uint16>& C)
{
	typedef SparseMultilineMatrix<uint16> Matrix;
	typename Matrix::RowIterator i_C;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[SEQUENTIAL] NB THREADS " << endl;

#ifdef SHOW_PROGRESS
	uint32 i=0;
	//std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	uint32 C_coldim = C.coldim ();
	uint64 *tmpDenseArray1C __attribute__((aligned(0x1000)));
		tmpDenseArray1C = new uint64[C_coldim];
	uint64 *tmpDenseArray2C __attribute__((aligned(0x1000)));
		tmpDenseArray2C = new uint64[C_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1CTimer);
	TIMER_DECLARE_(Axpy2CTimer);
	TIMER_DECLARE_(AxpyOuterTimer);

	//int x=0;
	for(i_C = C.rowBegin (); i_C != C.rowEnd (); ++i_C){
#ifdef SHOW_PROGRESS
		++i;
		report << "                                                                    \r";
		report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			razArray64(tmpDenseArray1C, C_coldim);
			razArray64(tmpDenseArray2C, C_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(*i_C, tmpDenseArray1C, tmpDenseArray2C);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		typename Ring::Element Cv1_col1=0, Cv2_col1=0;
		uint32 Cp1=0;
		typename Ring::Element Cv1_col2=0, Cv2_col2=0;
		uint32 Cp2=0;

		uint32 tmp=0;
		TIMER_START_(AxpyOuterTimer);
		for(uint32 j=i_C->IndexData[0]; j<C_coldim; ++j)
		{
			if(tmpDenseArray1C[j] % R._modulus == 0 && tmpDenseArray2C[j] % R._modulus == 0)
				continue;

			Cp1 = j;
			R.copy(Cv1_col1, tmpDenseArray1C[j] % R._modulus);
			R.copy(Cv2_col1, tmpDenseArray2C[j] % R._modulus);

			if (Cv1_col1 != 0)
				R.negin(Cv1_col1);
			if (Cv2_col1 != 0)
				R.negin(Cv2_col1);

			if (Cp1 % 2 == 0 && j < C_coldim - 1 )
			{
				Cp2 = j+1;
				R.copy(Cv1_col2, tmpDenseArray1C[j+1] % R._modulus);
				R.copy(Cv2_col2, tmpDenseArray2C[j+1] % R._modulus);

				uint16 val1 = A[Cp1 / C.nb_lines_per_bloc()].at_unchecked(0, 1);

				tmp = Cv1_col2 + (uint32)Cv1_col1*val1;
				R.init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1*val1;
				R.init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					R.negin(Cv1_col2);
				if (Cv2_col2 != 0)
					R.negin(Cv2_col2);

				TIMER_START_(Axpy2CTimer);
				axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
						A[Cp1 / C.nb_lines_per_bloc()],
						tmpDenseArray1C,
						tmpDenseArray2C);

				tmpDenseArray1C[Cp1] = Cv1_col1;
				tmpDenseArray2C[Cp1] = Cv2_col1;

				tmpDenseArray1C[Cp2] = Cv1_col2;
				tmpDenseArray2C[Cp2] = Cv2_col2;
				TIMER_STOP_(Axpy2CTimer);

				++j;
			}
			else
			{
				TIMER_START_(Axpy1CTimer);
					axpy(Cv1_col1, Cv2_col1,
						 A[Cp1 / C.nb_lines_per_bloc()],
						 Cp1 % C.nb_lines_per_bloc(),
						 tmpDenseArray1C,
						 tmpDenseArray2C);
					tmpDenseArray1C[Cp1] = Cv1_col1;
					tmpDenseArray2C[Cp1] = Cv2_col1;
				TIMER_STOP_(Axpy1CTimer);
			}
		}
		TIMER_STOP_(AxpyOuterTimer);

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArraysToMultilineVector64(R, tmpDenseArray1C, tmpDenseArray2C, C_coldim, *i_C, true);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif
		TIMER_REPORT_(RazArrayTimer);
		TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
		TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
		TIMER_REPORT_(Axpy1CTimer);
		TIMER_REPORT_(Axpy2CTimer);
		TIMER_REPORT_(AxpyOuterTimer);

		delete [] tmpDenseArray1C;
		delete [] tmpDenseArray2C;
}

template<typename Ring>
void reduceCD_onlyD(const Ring& R, const SparseMultilineMatrix<uint16>& C,
		const SparseMultilineMatrix<uint16>& B,
		SparseMultilineMatrix<uint16>& D)
{
	//same code
	reduceNonPivotsByPivots(R, C, B, D);
}


template<typename Ring>
void reduceCD(const Ring R, const SparseMultilineMatrix<uint16>& A,
							const SparseMultilineMatrix<uint16>& B,
							SparseMultilineMatrix<uint16>& C,
							SparseMultilineMatrix<uint16>& D)
{
	lela_check(B.coldim () == D.coldim ());
	lela_check(C.rowdim () == D.rowdim ());

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[SEQUENTIAL] NB THREADS " << endl;

#ifdef SHOW_PROGRESS
	uint32 i=0;
	//std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	typedef SparseMultilineMatrix<uint16> Matrix;

	typename Matrix::ConstRowIterator i_C;
	typename Matrix::RowIterator i_D;

	uint32 D_coldim = D.coldim();
	uint64 *tmpDenseArray1D __attribute__((aligned(0x1000)));
		tmpDenseArray1D = new uint64[D_coldim];
	uint64 *tmpDenseArray2D __attribute__((aligned(0x1000)));
		tmpDenseArray2D = new uint64[D_coldim];

	uint32 C_coldim = C.coldim ();
	uint64 *tmpDenseArray1C __attribute__((aligned(0x1000)));
		tmpDenseArray1C = new uint64[C_coldim];
	uint64 *tmpDenseArray2C __attribute__((aligned(0x1000)));
		tmpDenseArray2C = new uint64[C_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1CTimer);
	TIMER_DECLARE_(Axpy1DTimer);
	TIMER_DECLARE_(Axpy2CTimer);
	TIMER_DECLARE_(Axpy2DTimer);
	TIMER_DECLARE_(AxpyOuterTimer);

	int x=0;
	for(i_C = C.rowBegin (), i_D = D.rowBegin (); i_C != C.rowEnd (); ++i_C, ++i_D){

#ifdef SHOW_PROGRESS
		++i;
		report << "                                                                    \r";
		report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			razArray64(tmpDenseArray1D, D_coldim);
			razArray64(tmpDenseArray2D, D_coldim);
			razArray64(tmpDenseArray1C, C_coldim);
			razArray64(tmpDenseArray2C, C_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(D[x], tmpDenseArray1D, tmpDenseArray2D);
			copyMultilineToDenseArrays64(C[x], tmpDenseArray1C, tmpDenseArray2C);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		typename Ring::Element Cv1_col1=0, Cv2_col1=0;
		uint32 Cp1=0;
		typename Ring::Element Cv1_col2=0, Cv2_col2=0;
		//uint32 Cp2=0;

		uint32 tmp=0;

		TIMER_START_(AxpyOuterTimer);
		for(uint32 j=i_C->IndexData[0]; j<C_coldim; ++j)
		{
			if(tmpDenseArray1C[j] % R._modulus == 0 && tmpDenseArray2C[j] % R._modulus == 0)
				continue;

			Cp1 = j;
			R.copy(Cv1_col1, tmpDenseArray1C[j] % R._modulus);
			R.copy(Cv2_col1, tmpDenseArray2C[j] % R._modulus);

			if (Cv1_col1 != 0)
				R.negin(Cv1_col1);
			if (Cv2_col1 != 0)
				R.negin(Cv2_col1);


			if (Cp1 % 2 == 0 && j < C_coldim - 1 )
			{
				//Cp2 = j+1;
				R.copy(Cv1_col2, tmpDenseArray1C[j+1] % R._modulus);
				R.copy(Cv2_col2, tmpDenseArray2C[j+1] % R._modulus);

				uint16 val1 = A[Cp1 / C.nb_lines_per_bloc()].at_unchecked(0, 1);

				tmp = Cv1_col2 + (uint32)Cv1_col1*val1;
				R.init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1*val1;
				R.init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					R.negin(Cv1_col2);
				if (Cv2_col2 != 0)
					R.negin(Cv2_col2);

				++j;

				//if (Cv1_col2 != 0 || Cv2_col2 != 0) //have two consecutive lines
				//{
					TIMER_START_(Axpy2CTimer);
					axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
							A[Cp1 / C.nb_lines_per_bloc()],
							tmpDenseArray1C,
							tmpDenseArray2C);
					TIMER_STOP_(Axpy2CTimer);

					TIMER_START_(Axpy2DTimer);
					axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
							B[Cp1 / C.nb_lines_per_bloc()],
							tmpDenseArray1D,
							tmpDenseArray2D);
					TIMER_STOP_(Axpy2DTimer);
				/*}
				else
				{
					TIMER_START_(Axpy1CTimer);
						axpy(Cv1_col1, Cv2_col1,
						     A[Cp1 / C.nb_lines_per_bloc()],
							 Cp1 % C.nb_lines_per_bloc(),
							 tmpDenseArray1C,
							 tmpDenseArray2C);
					TIMER_STOP_(Axpy1CTimer);

					TIMER_START_(Axpy1DTimer);
						axpy(Cv1_col1, Cv2_col1,
							 B[Cp1 / C.nb_lines_per_bloc()],
							 Cp1 % C.nb_lines_per_bloc(),
							 tmpDenseArray1D,
							 tmpDenseArray2D);
					TIMER_STOP_(Axpy1DTimer);
				}*/
			}
			else
			{
				TIMER_START_(Axpy1CTimer);
					axpy(Cv1_col1, Cv2_col1,
						 A[Cp1 / C.nb_lines_per_bloc()],
						 Cp1 % C.nb_lines_per_bloc(),
						 tmpDenseArray1C,
						 tmpDenseArray2C);
				TIMER_STOP_(Axpy1CTimer);

				TIMER_START_(Axpy1DTimer);
						axpy(Cv1_col1, Cv2_col1,
							 B[Cp1 / C.nb_lines_per_bloc()],
							 Cp1 % C.nb_lines_per_bloc(),
							 tmpDenseArray1D,
							 tmpDenseArray2D);
					TIMER_STOP_(Axpy1DTimer);
			}
		}
		TIMER_STOP_(AxpyOuterTimer);

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArraysToMultilineVector64(R, tmpDenseArray1D, tmpDenseArray2D, D_coldim, D[x], true);
			//copyDenseArraysToMultilineVector64(R, tmpDenseArray1C, tmpDenseArray2C, C_coldim, C[x], true);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		x++;
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif
		TIMER_REPORT_(RazArrayTimer);
		TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
		TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
		TIMER_REPORT_(Axpy1CTimer);
		TIMER_REPORT_(Axpy2CTimer);
		TIMER_REPORT_(Axpy1DTimer);
		TIMER_REPORT_(Axpy2DTimer);
		TIMER_REPORT_(AxpyOuterTimer);

		delete [] tmpDenseArray1C;
		delete [] tmpDenseArray2C;
		delete [] tmpDenseArray1D;
		delete [] tmpDenseArray2D;
}

typedef struct ReduceCDParams_t {
	const Modular<uint16>* R;
	const SparseMultilineMatrix<uint16>* A;
	const SparseMultilineMatrix<uint16>* B;
	const SparseMultilineMatrix<uint16>* C;
	SparseMultilineMatrix<uint16>* D;
} ReduceCDParams_t;

static void* reduceCDParallel_in(void*);

static uint32 next_row_to_reduce;

#ifdef USE_MUTEX
static pthread_mutex_t mutex_lock;
#else
static pthread_spinlock_t spinlock_lock;
#endif

template<typename Ring>
void reduceCDParallel(const Ring& R, const SparseMultilineMatrix<uint16>& A,
							const SparseMultilineMatrix<uint16>& B,
							const SparseMultilineMatrix<uint16>& C,
							SparseMultilineMatrix<uint16>& D,
							int NUM_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[PARALLEL] NB THREADS " << NUM_THREADS << endl;

	next_row_to_reduce = 0;
	ReduceCDParams_t params;
	params.A = &A;
	params.B = &B;
	params.C = &C;
	params.D = &D;
	params.R = &R;

	pthread_t threads[NUM_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&mutex_lock, NULL);
#else
	pthread_spin_init(&spinlock_lock, NULL);
#endif

	for(t=0; t<NUM_THREADS; t++){
      report << "Creating thread " << t << "\n";
      rc = pthread_create(&threads[t], NULL, reduceCDParallel_in, &params);

      if (rc){		//TODO what happened when only one thread fails
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         exit(-1);
      }
   }

	/* Free attribute and wait for the other threads */
	for(t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}

		report << "COMPLETED: thread " << t << " - handled lines: " << (long)status << endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&mutex_lock);
#else
	pthread_spin_destroy(&spinlock_lock);
#endif
}

static void* reduceCDParallel_in(void* p_params)
{
#ifdef SHOW_PROGRESS
	uint32 i=0;
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	ReduceCDParams_t params = *(ReduceCDParams_t *)p_params;

	//typedef SparseMultilineMatrix<uint16> Matrix;

	//typename Matrix::ConstRowIterator i_C;
	//typename Matrix::RowIterator i_D;

	uint32 D_coldim = params.D->coldim();
	uint64 *tmpDenseArray1D;
		tmpDenseArray1D = new uint64[D_coldim];
	uint64 *tmpDenseArray2D;
		tmpDenseArray2D = new uint64[D_coldim];

	uint32 C_coldim = params.C->coldim ();
	uint64 *tmpDenseArray1C;
		tmpDenseArray1C = new uint64[C_coldim];
	uint64 *tmpDenseArray2C;
		tmpDenseArray2C = new uint64[C_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1CTimer);
	TIMER_DECLARE_(Axpy1DTimer);
	TIMER_DECLARE_(Axpy2CTimer);
	TIMER_DECLARE_(Axpy2DTimer);
	TIMER_DECLARE_(AxpyOuterTimer);
	TIMER_DECLARE_(AxpyLeftTimer);

	const uint32 nb_multiline_rows = (params.C->rowdim() / 2) + (params.C->rowdim() % 2);
	uint32 local_row_idx, nb_rows_handled=0;

#ifdef USE_MUTEX
	pthread_mutex_lock(&mutex_lock);
#else
	pthread_spin_lock(&spinlock_lock);
#endif

	while(next_row_to_reduce < nb_multiline_rows)
	{
		local_row_idx = next_row_to_reduce;
		++next_row_to_reduce;
		++nb_rows_handled;

#ifdef USE_MUTEX
	pthread_mutex_unlock(&mutex_lock);
#else
	pthread_spin_unlock(&spinlock_lock);
#endif

#ifdef SHOW_PROGRESS
		++i;
		report << "                                                                    \r";
		report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			razArray64(tmpDenseArray1D, D_coldim);
			razArray64(tmpDenseArray2D, D_coldim);
			razArray64(tmpDenseArray1C, C_coldim);
			razArray64(tmpDenseArray2C, C_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64((*params.D)[local_row_idx], tmpDenseArray1D, tmpDenseArray2D);
			copyMultilineToDenseArrays64((*params.C)[local_row_idx], tmpDenseArray1C, tmpDenseArray2C);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		uint16 Cv1_col1, Cv2_col1;
		uint32 Cp1;
		uint16 Cv1_col2, Cv2_col2;
	 	//uint32 Cp2;

		uint32 tmp;

		TIMER_START_(AxpyOuterTimer);
		for(uint32 j=(*params.C)[local_row_idx].IndexData[0]; j<C_coldim; ++j)
		{
			if(tmpDenseArray1C[j] % params.R->_modulus == 0 && tmpDenseArray2C[j] % params.R->_modulus == 0)
				continue;

			Cp1 = j;
			params.R->copy(Cv1_col1, tmpDenseArray1C[j] % params.R->_modulus);
			params.R->copy(Cv2_col1, tmpDenseArray2C[j] % params.R->_modulus);

			if (Cv1_col1 != 0)
				params.R->negin(Cv1_col1);
			if (Cv2_col1 != 0)
				params.R->negin(Cv2_col1);


			if (Cp1 % 2 == 0 && j < C_coldim - 1 )
			{
				//Cp2 = j+1;
				params.R->copy(Cv1_col2, tmpDenseArray1C[j+1] % params.R->_modulus);
				params.R->copy(Cv2_col2, tmpDenseArray2C[j+1] % params.R->_modulus);

				uint16 val1 = (*params.A)[Cp1 / params.C->nb_lines_per_bloc()].at_unchecked(0, 1);

				tmp = Cv1_col2 + (uint32)Cv1_col1*val1;
				params.R->init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1*val1;
				params.R->init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					params.R->negin(Cv1_col2);
				if (Cv2_col2 != 0)
					params.R->negin(Cv2_col2);

				++j;

				if (Cv1_col2 != 0 || Cv2_col2 != 0) //have two consecutive lines
				{
					TIMER_START_(Axpy2CTimer);
					axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
							(*params.A)[Cp1 / params.C->nb_lines_per_bloc()],
							tmpDenseArray1C,
							tmpDenseArray2C);
					TIMER_STOP_(Axpy2CTimer);

					TIMER_START_(Axpy2DTimer);
					axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
							(*params.B)[Cp1 / params.C->nb_lines_per_bloc()],
							tmpDenseArray1D,
							tmpDenseArray2D);
					TIMER_STOP_(Axpy2DTimer);
				}
				else
				{
					TIMER_START_(Axpy1CTimer);
						axpy(Cv1_col1, Cv2_col1,
						     (*params.A)[Cp1 / params.C->nb_lines_per_bloc()],
							 Cp1 % params.C->nb_lines_per_bloc(),
							 tmpDenseArray1C,
							 tmpDenseArray2C);
					TIMER_STOP_(Axpy1CTimer);

					TIMER_START_(Axpy1DTimer);
						axpy(Cv1_col1, Cv2_col1,
							 (*params.B)[Cp1 / params.C->nb_lines_per_bloc()],
							 Cp1 % params.C->nb_lines_per_bloc(),
							 tmpDenseArray1D,
							 tmpDenseArray2D);
					TIMER_STOP_(Axpy1DTimer);
				}
			}
			else
			{
				TIMER_START_(Axpy1CTimer);
					axpy(Cv1_col1, Cv2_col1,
						 (*params.A)[Cp1 / params.C->nb_lines_per_bloc()],
						 Cp1 % params.C->nb_lines_per_bloc(),
						 tmpDenseArray1C,
						 tmpDenseArray2C);
				TIMER_STOP_(Axpy1CTimer);

				TIMER_START_(Axpy1DTimer);
						axpy(Cv1_col1, Cv2_col1,
							 (*params.B)[Cp1 / params.C->nb_lines_per_bloc()],
							 Cp1 % params.C->nb_lines_per_bloc(),
							 tmpDenseArray1D,
							 tmpDenseArray2D);
					TIMER_STOP_(Axpy1DTimer);
			}
		}
		TIMER_STOP_(AxpyOuterTimer);

		copyDenseArraysToMultilineVector64(*params.R, tmpDenseArray1D, tmpDenseArray2D, D_coldim, (*params.D)[local_row_idx], true);
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif
		TIMER_REPORT_(RazArrayTimer);
		TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
		TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
		TIMER_REPORT_(Axpy1CTimer);
		TIMER_REPORT_(Axpy2CTimer);
		TIMER_REPORT_(Axpy1DTimer);
		TIMER_REPORT_(Axpy2DTimer);
		TIMER_REPORT_(AxpyOuterTimer);
		TIMER_REPORT_(AxpyLeftTimer);

		delete [] tmpDenseArray1C;
		delete [] tmpDenseArray2C;
		delete [] tmpDenseArray1D;
		delete [] tmpDenseArray2D;

		return (void*) nb_rows_handled;
}



template <typename Index>
static inline long head(const MultiLineVector<uint16, Index>& v, const uint16 line_idx, uint16& a, uint32& head_idx)
{
	uint16 val=0;
	for(uint32 i=0; i<v.size(); ++i)
	{
		if((val = v.at_unchecked(line_idx, i)) != 0)
		{
			a = val;
			head_idx = i;
			return (long)v.IndexData[i];
		}
	}

	return -1;
}

template <typename Ring>
static long head(const Ring& R, uint64 arr[], const size_t size, uint16& a)
{
	for(uint32 i=0; i<size; ++i)
	{
		if(arr[i] % R._modulus != 0)
		{
			a = arr[i] % R._modulus;
			return (long)i;
		}
		else		//reduce on the go
			arr[i] = 0;
	}

	return -1;
}

//WARNING could fail on non prime rings
template <typename Ring, typename Index>
static inline void normalize_multiline(const Ring& R, MultiLineVector<uint16, Index>& v)
{
	uint32 idx;
	uint16 h1=0, h2=0;
	head(v, 0, h1, idx);
	head(v, 1, h2, idx);

	if(v.empty ())
		return;

	if(h1 != 0)
		if(R.invin(h1) != true)
			throw "Non Invertible Value";

	if(h2 != 0)
		if(R.invin(h2) != true)
			throw "Non Invertible Value";

	//TODO: skip when both are 1 or 0
	if((h1==0 || h1==1) && (h2==0 || h2==1))
		return;

	if(h1 == 0 || h1 == 1)
	{
		for(uint32 i=0; i<v.size (); ++i)
		{
			R.mulin(v.ValuesData[1+i*2], h2);
		}
	}
	else if(h2 == 0 || h2 == 1)
	{
		for(uint32 i=0; i<v.size (); ++i)
		{
			R.mulin(v.ValuesData[0+i*2], h1);
		}
	}
	else
	{
		for(uint32 i=0; i<v.size (); ++i)
		{
			R.mulin(v.ValuesData[0+i*2], h1);
			R.mulin(v.ValuesData[1+i*2], h2);
		}
	}
}

template <typename Ring>
static void normalize_array64(const Ring& R, uint64 arr[], const size_t size)
{
	uint16 a;
	int h1 = head(R, arr, size, a);

	if(h1 == -1)
		return;

	R.invin(a);
	for(uint32 i=h1; i<size; ++i)
	{
		arr[i] *= a;
		arr[i] %= R._modulus;
	}
}

template <typename Element>
static uint32 echelonize(const Modular<Element>& R, SparseMultilineMatrix<Element>& A)
{
	uint32 coldim = A.coldim ();
	uint32 npiv = 0;
	uint32 npiv_real = 0;
	uint32 N = A.rowdim()/2 + A.rowdim()%2;
	uint32 piv[N];
	for(uint32 i=0; i<N; ++i)
		piv[i] = 0;

	typedef SparseMultilineMatrix<uint16> Matrix;
	MultiLineVector<uint16, uint32> *rowA;

#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16>" << std::endl;
#endif
	typedef Modular<uint16> Ring;
	typename Matrix::RowIterator i_A;

	uint64 *tmpDenseArray1;
		tmpDenseArray1 = new uint64[coldim];

	uint64 *tmpDenseArray2;
		tmpDenseArray2 = new uint64[coldim];


	uint32 i=0;
	uint16 h=0;

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

	TIMER_DECLARE_(AxpyTimer);
	TIMER_DECLARE_(Axpy2Timer);

	TIMER_START_(normalizeVectorTimer);
		normalize_multiline(R, A[0]);
	TIMER_STOP_(normalizeVectorTimer);

	if(A.rowdim() == 0)
		return 0;

	for(i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++i)
	{
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << std::ends;
#endif

        TIMER_START_(RazArrayTimer);
        	razArray64(tmpDenseArray1, coldim);
        	razArray64(tmpDenseArray2, coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(*i_A, tmpDenseArray1, tmpDenseArray2);
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
				head_line1 = head(*rowA, 0, h_a1, head_line1_idx);
				head_line2 = head(*rowA, 1, h_a2, head_line2_idx);
			TIMER_STOP_(HeadVectorTimer);

			if(head_line1 != -1 && head_line1 == head_line2)
				throw "Wrong Mutiline format";
			if(head_line1 > head_line2)			//makes the row with the smallest column entry first in the multiline
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

			if(head_line1 != -1)
			{
				R.copy(v1col1, tmpDenseArray1[head_line1] % R._modulus);
				R.copy(v2col1, tmpDenseArray2[head_line1] % R._modulus);

				if(v1col1 != 0)
					R.negin(v1col1);
				if(v2col1 != 0)
					R.negin(v2col1);
			}
			else
			{
				v1col1 = 0;
				v2col1 = 0;
			}

			if(head_line2 != -1)
			{
				R.copy(v1col2, tmpDenseArray1[head_line2] % R._modulus);
				R.copy(v2col2, tmpDenseArray2[head_line2] % R._modulus);
			}
			else
			{
				v1col2 = 0;
				v2col2 = 0;
			}

			//TODO: check if buggy! if line is empty for example
			uint16 val = rowA->at_unchecked(0, head_line2_idx);

			tmp = v1col2 + (uint32)v1col1*val;
			R.init(v1col2, tmp);
			tmp = v2col2 + (uint32)v2col1*val;
			R.init(v2col2, tmp);

			if (v1col2 != 0)
				R.negin(v1col2);
			if (v2col2 != 0)
				R.negin(v2col2);

			if((v1col1 == 0 && v2col1 == 0) && (v1col2 != 0 || v2col2 != 0))
			{
				TIMER_START_(AxpyTimer);
				axpy(v1col2, v2col2,
					*rowA,
					1,	//reduce by second line only
					tmpDenseArray1,
					tmpDenseArray2);
				TIMER_STOP_(AxpyTimer);
			}else if((v1col2 == 0 && v2col2 == 0) && (v1col1 != 0 || v2col1 != 0))
			{
				TIMER_START_(AxpyTimer);
				axpy(v1col1, v2col1,
						*rowA,
						0,	//reduce by first line only
						tmpDenseArray1,
						tmpDenseArray2);
				TIMER_STOP_(AxpyTimer);
			}
			else
			{
				TIMER_START_(Axpy2Timer);
				axpy2(v1col1, v2col1, v1col2, v2col2,
						*rowA,
						tmpDenseArray1,
						tmpDenseArray2);
				TIMER_STOP_(Axpy2Timer);
			}
		}

		TIMER_START_(normalizeArrayTimer);
			normalize_array64(R, tmpDenseArray1, coldim);
			normalize_array64(R, tmpDenseArray2, coldim);
		TIMER_STOP_(normalizeArrayTimer);

		TIMER_START_(HeadArrayTimer);
			head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
			head_line2 = head(R, tmpDenseArray2, coldim, h_a2);
        TIMER_STOP_(HeadArrayTimer);

        //assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

        //reduce by same multiline
        if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
        {
        	assert(head_line1 >= 0 && head_line1 < coldim);
        	if(tmpDenseArray2[head_line1] % R._modulus != 0)
			{
				h = tmpDenseArray2[head_line1] % R._modulus;
				R.negin(h);
			}

			for(x = head_line1; x<coldim; ++x)
			{
				tmpDenseArray2[x] += (uint32)h * tmpDenseArray1[x];
			}
        }

		///////////////////////////////////////////////////////////////////////////////////////////

        TIMER_START_(HeadArrayTimer);
			//head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
			head_line2 = head(R, tmpDenseArray2, coldim, h_a2);		//only this can change
        TIMER_STOP_(HeadArrayTimer);

        TIMER_START_(CopyDenseArrayToSparseVectorTimer);
        	if((head_line2 != -1) && (head_line1 > head_line2))			//saves the line with the smallest column entry first
        		copyDenseArraysToMultilineVector64(R, tmpDenseArray2, tmpDenseArray1, coldim, *i_A, true);
        	else
        		copyDenseArraysToMultilineVector64(R, tmpDenseArray1, tmpDenseArray2, coldim, *i_A, true);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

		TIMER_START_(normalizeVectorTimer);
			normalize_multiline(R, *i_A);
		TIMER_STOP_(normalizeVectorTimer);

		if(!i_A->empty())
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

	TIMER_REPORT_(AxpyTimer);
	TIMER_REPORT_(Axpy2Timer);

	delete [] tmpDenseArray1;
	delete [] tmpDenseArray2;

	return npiv_real;
}


template <typename Ring>
static void copy(const Ring& R,
		const SparseMultilineMatrix<typename Ring::Element>& A,
		SparseMatrix<typename Ring::Element>& B)
{
	std::vector<SparseVector<typename Ring::Element> > v (A.nb_lines_per_bloc());

	typename SparseMultilineMatrix<typename Ring::Element>::ConstRowIterator i_A = A.rowBegin();
	typename SparseMatrix<typename Ring::Element>::RowIterator i_B = B.rowBegin();

	typedef SparseMatrix<typename Ring::Element> Matrix;

	Context<Ring> ctx(R);

	lela_check(A.rowdim() == B.rowdim());

	while (i_A != A.rowEnd()) //for each multiline of A
	{
		for (uint32 j = 0; j < i_A->size(); ++j)
		{
			for (uint16 i = 0; i < A.nb_lines_per_bloc(); ++i) //for each line in the multiline
			{
				if (i_A->at(i, j) != 0)
					v[i].push_back(typename Matrix::Row::value_type(i_A->IndexData[j], i_A->at(i, j)));
			}
		}

		for (uint16 i = 0; i < A.nb_lines_per_bloc(); ++i) //for each line in the multiline
		{
			if (v[i].empty() && i_B == B.rowEnd())
				return;

			if(i_B == B.rowEnd())
				throw "Non empty lines beyond B dimension detected";

			i_B->swap(v[i]);
			v[i].clear();
			++i_B;
		}

		++i_A;
	}
}








}









#endif //MATRIX_OP_MULTILINE_
