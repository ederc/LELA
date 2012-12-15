/*
 * matrix-op.C
 *
 *  Created on: 19 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_OP_MULTILINE_
#define MATRIX_OP_MULTILINE_

#include "FG-types.h"
#include "lela/ring/modular.h"

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
	for(x=xl; x<sz; x+=STEP)			\
	{

#define MIDDLE_UNROLL_CODE2				\
		for(uint8 t=0; t<STEP; ++t)		\
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

static void copyMultilineToDenseArrays64(const MultiLineVector<uint16>& v,
		uint64 *arr1, uint64 *arr2)
{
	register uint32 idx;
	/*for (uint32 i = 0; i < v.size(); ++i)
	{
		idx = v.IndexData[i];
		arr1[idx] = (uint64) v.at_unchecked(0, i);
		arr2[idx] = (uint64) v.at_unchecked(1, i);
	}*/

	const uint8 STEP = 32;
	uint8 xl = v.size() % STEP;
	uint32 sz = v.size();
	register uint32 x = 0;
	while (x < xl)
	{
		idx = v.IndexData[x];
		arr1[idx] = (uint64) v.at_unchecked(0, x);
		arr2[idx] = (uint64) v.at_unchecked(1, x);

		++x;
	}

	for (x = xl; x < sz; x += STEP)
	{

#pragma loop unroll
		for (uint8 t = 0; t < STEP; ++t)
		{
			idx = v.IndexData[x + t];
			arr1[idx] = (uint64) v.at_unchecked(0, x + t);
			arr2[idx] = (uint64) v.at_unchecked(1, x + t);
		}
	}
}

template<typename Ring>
static void copyDenseArraysToMultilineVector64(const Ring& R,
		const uint64 *arr1,
		const uint64 *arr2,
		const uint32 size,
		MultiLineVector<uint16>& v,
		bool reduce)
{
	MultiLineVector<uint16> tmp(v.nb_lines());
	typename Ring::Element e1, e2;


	uint32 sz=0;
	for (uint32 i = 0; i < size; ++i)
		if ((arr1[i] != 0) || (arr2[i] != 0))
			sz++;

	tmp.reserve(sz);

	if (reduce)
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
	else
		for (uint32 i = 0; i < size; ++i)
		{
			if ((arr1[i] != 0) || (arr2[i] != 0))
			{
				tmp.IndexData.push_back(i);
				tmp.ValuesData.push_back((uint16) arr1[i]);
				tmp.ValuesData.push_back((uint16) arr2[i]);
			}
		}

	v.swap(tmp);
}

inline void axpy2(const uint16 av1_col1,
		const uint16 av2_col1,
		const uint16 av1_col2,
		const uint16 av2_col2,
		const MultiLineVector<uint16>& v,
		uint64 *arr1,
		uint64 *arr2)
{
	register uint32 idx;
	register uint16 val1, val2;


	/*for(uint32 i=0; i<sz; ++i)
	{
		idx = v.IndexData[i];

		val1 = v.at_unchecked(0, i);
		val2 = v.at_unchecked(1, i);

		arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
		arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
	}*/

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

	if (av1_col1 == 0)
	{
		if (av1_col2 == 0)
			for (uint32 i = 0; i < v.size(); ++i)
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
			for (uint32 i = 0; i < v.size(); ++i)
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
			for (uint32 i = 0; i < v.size(); ++i)
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
			for (uint32 i = 0; i < v.size(); ++i)
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
			for (uint32 i = 0; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));

				//arr2[idx] += (uint32) av2_col1 * val1;
				//arr2[idx] += (uint32) av2_col2 * val2;
			}
		else
			for (uint32 i = 0; i < v.size(); ++i)
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
			for (uint32 i = 0; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint32) av1_col1 * val1;
				//arr1[idx] += (uint32) av1_col2 * val2;

				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
		else if (av2_col2 == 0)
			for (uint32 i = 0; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));

				arr2[idx] += (uint32) av2_col1 * val1;
				//arr2[idx] += (uint32) av2_col2 * val2;
			}
		else
			for (uint32 i = 0; i < v.size(); ++i)
			{
				idx = v.IndexData[i];

				val1 = v.at_unchecked(0, i);
				val2 = v.at_unchecked(1, i);

				arr1[idx] += (uint64)(((uint32) av1_col1 * val1) + (uint64)((uint32) av1_col2 * val2));
				arr2[idx] += (uint64)(((uint32) av2_col1 * val1) + (uint64)((uint32) av2_col2 * val2));
			}
	}
}

inline void axpy(const uint16 av1_col1,
		const uint16 av2_col1,
		const MultiLineVector<uint16>& v,
		const uint16 line,
		uint64 *arr1,
		uint64 *arr2)
{
	lela_check(line < v.nb_lines());
	const uint8 STEP = 16;
	const uint32 sz = v.size();

	register uint32 idx;
	register uint16 val1;

	/*for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		val1 = v.at_unchecked(line, i);

		arr1[idx] += (uint32) av1_col1 * val1;
		arr2[idx] += (uint32) av2_col1 * val1;
	}*/


	uint8 xl = v.size() % STEP;
	register uint32 x = 0;
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

		for (x = xl; x < sz; x += STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < STEP; ++t)
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

		for (x = xl; x < sz; x += STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < STEP; ++t)
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

		for (x = xl; x < sz; x += STEP)
		{

#pragma loop unroll
			for (uint8 t = 0; t < STEP; ++t)
			{
				idx = v.IndexData[x + t];
				val1 = v.at_unchecked(line, x + t);
				arr2[idx] += (uint32) av2_col1 * val1;
			}
		}
	}

}

template<typename Ring>
void reducePivotsByPivots(const Ring& R, const SparseMultilineMatrix<uint16>& A,
		SparseMultilineMatrix<uint16>& B)
{
	lela_check(A.coldim () == B.coldim ());
	lela_check(A.rowdim () == A.coldim ());

	//typedef Modular<uint16> Ring;
	typedef SparseMultilineMatrix<uint16> Matrix;

	typename Matrix::ConstRowIterator i_A;
	typename Matrix::RowIterator i_B;

	uint32 B_coldim = B.coldim();
	uint64 tmpDenseArray1[B_coldim] __attribute__((aligned(64)));
	uint64 tmpDenseArray2[B_coldim] __attribute__((aligned(64)));

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(Axpy1Timer);
	TIMER_DECLARE_(Axpy2Timer);
	TIMER_DECLARE_(AxpyOuterTimer);
	TIMER_DECLARE_(AxpyLeftTimer);

#ifdef SHOW_PROGRESS
	uint32 i=A.rowdim () / A.nb_lines_per_bloc ();
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
				Ap2 = i_A->IndexData[j + 1];
				if (Ap2 == Ap1 + 1) //have two consecutive lines
				{
					R.copy(Av1_col2, i_A->at_unchecked(0, j + 1));
					R.copy(Av2_col2, i_A->at_unchecked(1, j + 1));

					if (Av1_col2 != 0)
						R.negin(Av1_col2);
					if (Av2_col2 != 0)
						R.negin(Av2_col2);

					++j;

					TIMER_START_(Axpy2Timer);

					axpy2(Av1_col1, Av2_col1, Av1_col2, Av2_col2,
							B[Ap1 / A.nb_lines_per_bloc()], tmpDenseArray1,
							tmpDenseArray2);
					TIMER_STOP_(Axpy2Timer);
				}
				else
				{
					TIMER_START_(Axpy1Timer);

					axpy(Av1_col1, Av2_col1, B[Ap1 / A.nb_lines_per_bloc()],
							Ap1 % A.nb_lines_per_bloc(), tmpDenseArray1,
							tmpDenseArray2);
					TIMER_STOP_(Axpy1Timer);
				}
			}
			else
			{
				TIMER_START_(Axpy1Timer);

				axpy(Av1_col1, Av2_col1, B[Ap1 / A.nb_lines_per_bloc()],
						Ap1 % A.nb_lines_per_bloc(), tmpDenseArray1,
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
		copyDenseArraysToMultilineVector64(R, tmpDenseArray1, tmpDenseArray2,
				B_coldim, *i_B, reduce);
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
							 SparseMultilineMatrix<uint16>& D)
{
	lela_check(B.coldim () == D.coldim ());
	lela_check(C.rowdim () == D.coldim ());

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

			if (Cv1_col1 != 0)
				R.negin(Cv1_col1);
			if (Cv2_col1 != 0)
				R.negin(Cv2_col1);

			if (Cp1 % 2 == 0 && j < i_C->size() - 1)
			{
				assert(j != i_C->size() - 1);
				Cp2 = i_C->IndexData[j + 1];
				if (Cp2 == Cp1 + 1) //have two consecutive lines
				{
					R.copy(Cv1_col2, i_C->at_unchecked(0, j + 1));
					R.copy(Cv2_col2, i_C->at_unchecked(1, j + 1));

					if (Cv1_col2 != 0)
						R.negin(Cv1_col2);
					if (Cv2_col2 != 0)
						R.negin(Cv2_col2);

					++j;

					TIMER_START_(Axpy2Timer);
					axpy2(Cv1_col1, Cv2_col1, Cv1_col2, Cv2_col2,
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

		copyDenseArraysToMultilineVector64(R, tmpDenseArray1, tmpDenseArray2,
				D_coldim, D[x], true);
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


}









#endif //MATRIX_OP_MULTILINE_
