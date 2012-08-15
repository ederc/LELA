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

using namespace LELA;

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

#define UNROLL_STEP 32

/**
 * Given a bloc of dense rows (an array of arrays), Zeros its memory
 */
template <typename Element>
inline void memsetBlocToZero(Element** arr, const uint32 nb_lines, const uint32 line_size)
{
	//cout << "in memset" << endl;
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		memset(arr[i], 0, DEFAULT_BLOC_WIDTH * sizeof(Element));
	}
}

/**
 * Copy a SparseBloc to a bloc of dense rows (an array or arrays)
 */
template<typename Ring, typename Index, typename DoubleFlatElement>
inline void copySparseBlocToDenseBlocArray(const Ring& R,
		const SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		DoubleFlatElement** arr)
{
	//Context<Ring> ctx (R);
	//cout << "in copySparseBlocToDenseBlocArray" << endl;
	//for all the rows in the bloc
//#pragma loop unroll
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		//BLAS1::write(ctx, cout, bloc[i]);
		//for all the elements in the row bloc[i]
		Index idx;
		typename Ring::Element val1, val2;

		//cout << "line " << i << endl;

		if(bloc[i].empty ())
			continue;
		//cout << "line__ " << i << endl;

		/*cout << "DENSE\n";
		for(uint32 j=0; j<bloc[i].size (); ++j)
		{

			cout << (uint32)bloc[i].IndexData[j] << ", ";
			cout << (uint32)bloc[i].at(0, j) << ", ";
			cout << (uint32)bloc[i].at(1, j) << endl;
		}
		cout << "DENSE\n";*/

		for(uint32 j=0; j<bloc[i].size (); ++j)
		{

			idx = bloc[i].IndexData[j];
			val1 = bloc[i].at_unchecked(0, j);
			val2 = bloc[i].at_unchecked(1, j);

			/*cout << "idx " << (uint32)idx << " v1 " << (uint32)val1 << " v2 " << (uint32)val2 << endl;

			if(idx >= DEFAULT_BLOC_WIDTH)
				cout << "EROOOOOOOOOOOOOOOOOOOOOOOOOOOOR\n";*/

			arr[i*2][idx] = val1;
			arr[i*2+1][idx] = val2;
		}

		//cout << endl << "*********** dense *************" << endl;
		//for(uint32 j=0; j<bloc.bloc_width(); ++j)
			//cout << arr[i][j] << ", ";

		//cout << endl;
	}
	//cout << "end" << endl;
}

template <typename Ring, typename DoubleFlatElement, typename Index>
inline void copyDenseBlocArrayToSparseBloc(const Ring& R, DoubleFlatElement** arr, SparseMultilineBloc<typename Ring::Element, Index>& bloc)
{
	typename Ring::Element e1, e2;
	//typename SparseBloc<typename Ring::Element>::Row tmp;

	//Context<Ring> ctx (R);
//cout << "in copyDenseBlocArrayToSparseBloc\n";
	//this shoud be a global
	//typename Ring::Element *reduced_non_zero_elts = new typename Ring::Element [bloc.bloc_width ()];
	//TODO: pre allocate memory


	//for each row in the bloc
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		//tmp = typename SparseBloc<typename Ring::Element>::Row ();
		bloc[i].clear ();

		for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
		{
			//TODO: in generic method
			//R.init(e, arr[i][j])
			//if (!R.isZero(e))
			ModularTraits<typename Ring::Element>::reduce (e1, arr[i*2][j], R._modulus);
			ModularTraits<typename Ring::Element>::reduce (e2, arr[i*2+1][j], R._modulus);

			if (!R.isZero(e1) || !R.isZero(e2))
			{
				bloc[i].IndexData.push_back (j); // (typename SparseBloc<typename Ring::Element>::Row::value_type (j, e));
				bloc[i].ValuesData.push_back (e1);
				bloc[i].ValuesData.push_back (e2);
			}
		}

		//BLAS1::write(ctx, cout, bloc[i]);
		//point bloc[i] memory to the new elements
		//bloc[i].swap(tmp);
	}
//cout << "end \n";
}

template <typename Ring, typename Row >
inline void copyDenseArrayToSparseRow(const Ring& R, uint64* arr1, uint64* arr2, const uint32 line_size, Row& v, bool reduce)
{

	//Row tmp;

	v.clear ();
	if(v.size() != 0)
		cout << "EROOOOR\n";
	//if(v.capacity () < 128)
		//v.reserve (128);

	//register uint16 l;
	register typename Ring::Element e1, e2;

	if(reduce)
		for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
		{
			//ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
			//l = arr[j] >> 16;
			//l = R._modulus - l;

			ModularTraits<typename Ring::Element>::reduce (e1, arr1[j], R._modulus);
			ModularTraits<typename Ring::Element>::reduce (e2, arr2[j], R._modulus);

			if ((!R.isZero(e1)) || (!R.isZero(e2)))
			{
				v.IndexData.push_back (j); // (typename SparseBloc<typename Ring::Element>::Row::value_type (j, e));
				v.ValuesData.push_back (e1);
				v.ValuesData.push_back (e2);
			}
		}
	else
		for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
		{
			if (arr1[j] != 0 || arr2[j] != 0)
			{
				v.IndexData.push_back (j); // (typename SparseBloc<typename Ring::Element>::Row::value_type (j, e));
				v.ValuesData.push_back (arr1[j]);
				v.ValuesData.push_back (arr2[j]);
			}
		}

	//v.swap(tmp);
}

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

/*template <typename Ring>
inline void reduceDenseArrayModulo(const Ring& R, uint64* arr, const uint32 line_size, int a)
{
	typename Ring::Element e;

	for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
	{
		ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
		arr[j] = (uint64)e;
	}
}*/


/**
 * acc <- acc + a.v
 * Note: caller must ensure that acc is large enough (>= x.size())
 */
/*template <typename BlocSparseVector>
inline void axpy(const BlocSparseVector& v,
		const uint16 a, uint64* arr) __attribute__((always_inline));

template <typename BlocSparseVector>
inline void axpy(const BlocSparseVector& v,
		const uint16 a, uint64* arr)
{
	//typename BlocSparseVector::const_iterator it;
	//uint32 a_32 = a;

	//for(it = x.begin (); it != x.end (); ++it)
	//{
		//acc[it->first] += (uint64)(a_32 * it->second);
	//}

#define STEP 32

	register uint32 x=0;
	while(x<v.size () % STEP)
	{
		arr[v[x].first] += (uint32)a * v[x].second;
		++x;
	}

	for(x=v.size () % STEP; x<v.size (); x+=STEP)
	{
		for(uint8 t=0; t<STEP; t+=2)
		{
			arr[v[x+t].first] += (uint32)a * v[x+t].second;
			arr[v[x+t+1].first] += (uint32)a * v[x+t+1].second;
		}
			//R.axpyin(tmpDenseArray[(row_it_B+t)->first], 	 Cv, (row_it_B+t)->second);
	}

#undef STEP
}*/

/*template<typename Ring>
inline void axpy(const Ring& R, const uint64 *arr1,
		const typename Ring::Element a, uint64* arr2, const uint32 line_size, bool t)
{
	register typename Ring::Element e;
	register uint64 acc;

	for (uint32 i = 0; i < DEFAULT_BLOC_WIDTH; ++i)
	{
		//TODO: try convert to uint16 and do multiplication with uint32 since we know arr1 is already reduced
		acc = arr1[i] * a + arr2[i];

		ModularTraits<typename Ring::Element>::reduce(e, acc, R._modulus);
		arr2[i] = (uint64) e;
	}
}*/

template <typename Index>
inline void axpy(const uint16 av1_col1,
		const uint16 av2_col1,
		const MultiLineVector<uint16, Index>& v,
		const uint16 line,
		uint64 *arr1,
		uint64 *arr2)
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
	register uint32 x = 0;

	if(av1_col1 == 0 && av2_col1 == 0)
		return;

	if (av1_col1 != 0 && av2_col1 != 0) //cannot both be 0
	{
		register uint32 idx;
		register uint16 val1;
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

//#pragma loop unroll
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
		register uint32 idx;
		register uint16 val1;
		while (x < xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr1[idx] += (uint32) av1_col1 * val1;

			++x;
		}

		for (x = xl; x < sz; x += UNROLL_STEP)
		{

//#pragma loop unroll
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
		register uint32 idx;
		register uint16 val1;
		while (x < xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr2[idx] += (uint32) av2_col1 * val1;

			++x;
		}

		for (x = xl; x < sz; x += UNROLL_STEP)
		{

//#pragma loop unroll
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
		uint64 *arr2)
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
		axpy(av1_col2, av2_col2, v, 1, arr1, arr2);
		return;
	}

	if(av1_col2 == 0 && av2_col2 == 0)
	{
		axpy(av1_col1, av2_col1, v, 0, arr1, arr2);
		return;
	}

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


TIMER_DECLARE_(reduceDenseBlocModulo);
TIMER_DECLARE_(axpyInTriangular);
TIMER_DECLARE_(axpyInRectangular);
TIMER_DECLARE_(copyDenseArrayToSparseRow);
TIMER_DECLARE_(Axpy1);
TIMER_DECLARE_(Axpy2);


template <typename Index>
inline void reduceBlocByRectangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		const SparseMultilineBloc<uint16, Index>& bloc_B,
		uint64 **Bloc_acc)
{
	//cout << "REDUCE RECTANGULAR nb_lines " << bloc_A.size() << endl;
	typedef Modular<uint16> Ring;

	if(bloc_A.empty() || bloc_B.empty())
		return;

	//for all the rows in Bloc_acc (same as for Bloc_A)
	for(uint16 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		if(bloc_A[i].empty())
			continue;

		typename Ring::Element Av1_col1, Av2_col1;
		uint32 Ap1;
		typename Ring::Element Av1_col2, Av2_col2;
		uint32 Ap2;

		//for all element in row bloc_A[i]
		//for(uint16 j=0; j<bloc_A[i].size(); ++j)
		//{
			/*register uint16 Cv;
			register typename SparseBloc<uint16, Index>::IndexType Cp;

			R.copy(Cv, bloc_A[i][j].second);
			Cp = bloc_A[i][j].first;
			R.negin(Cv);

			TIMER_START_(axpyInRectangular);
				axpy(bloc_B[Cp], Cv, Bloc_acc[i]);
			TIMER_STOP_(axpyInRectangular);*/

		for (uint32 j = 0; j < bloc_A[i].size(); ++j)
		{
			Ap1 = bloc_A[i].IndexData[j];
			R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
			R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));

			if (Av1_col1 != 0)
				R.negin(Av1_col1);
			if (Av2_col1 != 0)
				R.negin(Av2_col1);

			if (Ap1 % 2 == 0 && j < bloc_A[i].size() - 1)
			{
				assert(j != bloc_A[i].size() - 1);
				Ap2 = bloc_A[i].IndexData[j + 1];
				if (Ap2 == Ap1 + 1) //have two consecutive lines
				{
					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j + 1));
					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j + 1));

					if (Av1_col2 != 0)
						R.negin(Av1_col2);
					if (Av2_col2 != 0)
						R.negin(Av2_col2);

					++j;

					//cout << "axpy 2 \n";
					TIMER_START_(Axpy2);
					axpy2(Av1_col1, Av2_col1, Av1_col2, Av2_col2,
							bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
							Bloc_acc[i*2],
							Bloc_acc[i*2+1]);
					TIMER_STOP_(Axpy2);
				}
				else
				{
					TIMER_START_(Axpy1);
					axpy(Av1_col1, Av2_col1, bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
							Ap1 % NB_ROWS_PER_MULTILINE,
							Bloc_acc[i*2],
							Bloc_acc[i*2+1]);
					TIMER_STOP_(Axpy1);
				}
			}
			else
			{
				TIMER_START_(Axpy1);
				axpy(Av1_col1, Av2_col1, bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
						Ap1 % NB_ROWS_PER_MULTILINE,
						Bloc_acc[i*2],
						Bloc_acc[i*2+1]);
				TIMER_STOP_(Axpy1);
			}
		}
	}
}


/**
 * Reduce the rows inside the bloc by themselves
 */
template<typename Index>
inline void reduceBlocByTriangularBloc(const Modular<uint16>& R,
		const SparseMultilineBloc<uint16, Index>& bloc_A,
		uint64** Bloc_acc,
		SparseMultilineBloc<uint16, Index>& bloc_B)
{
	typedef Modular<uint16> Ring;

	//if(bloc_A.empty())
		//return;

	if(bloc_A.size () < 2)
		return;

	//uint32 last_reduced_line = 0;

	//reduce first multiline
	//elements in first multiline of bloc_A are of the form
	// |1 x|
	// |  1|

	for(uint16 i=0; i<DEFAULT_BLOC_HEIGHT/2; ++i)
	{
		//cout << "line "<< i << endl;
		if(bloc_A[i].empty())	// can be empty only if it's the last bloc
		{
			TIMER_START_(copyDenseArrayToSparseRow);
				copyDenseArrayToSparseRow(R, Bloc_acc[i*2], Bloc_acc[i*2+1], bloc_B.bloc_width(), bloc_B[i], true);
			TIMER_STOP_(copyDenseArrayToSparseRow);
			continue;
		}

		/*TIMER_START_(copyDenseArrayToSparseRow);
			copyDenseArrayToSparseRow(R, Bloc_acc[last_reduced_line], bloc_B.bloc_width(), bloc_B[last_reduced_line]);
		TIMER_STOP_(copyDenseArrayToSparseRow);

		last_reduced_line = i;

		for(int j=0; j<bloc_A[i].size () - 1; ++j)
		{
			register uint16 Av;
			register typename SparseBloc<uint16, Index>::IndexType Ap;

			Ap = bloc_A[i][j].first;

			R.copy(Av, bloc_A[i][j].second);
			R.negin(Av);

			//lela_check(Ap < i);	//should point to a a row previous to the ith one in the bloc

			TIMER_START_(axpyInTriangular);
				axpy(bloc_B[Ap], Av, Bloc_acc[i]);
				//last_reduced_line = i;
			TIMER_STOP_(axpyInTriangular);
		}*/

		//on last multiline row, half of the row can be empty, we have to reach size()-1 elements
		int last_idx = -1;
		if(bloc_A[i].at_unchecked(1, bloc_A[i].size()-1) == 0)
			last_idx = (int)bloc_A[i].size()-1;
		else
			last_idx = (int)bloc_A[i].size()-2;

		typename Ring::Element Av1_col1, Av2_col1;
		uint32 Ap1;
		typename Ring::Element Av1_col2, Av2_col2;
		uint32 Ap2;

		for (int j = 0; j < last_idx; ++j) //skip first two elements
		{
			Ap1 = (uint32)bloc_A[i].IndexData[j];
			R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
			R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));

			if (Av1_col1 != 0)
				R.negin(Av1_col1);
			if (Av2_col1 != 0)
				R.negin(Av2_col1);

			if(Ap1 >= i*2)
				throw std::runtime_error ("Index pointing to an out of range line.");

			if (Ap1 % 2 == 0)
			{
				//assert(j != bloc_A[i].size() - 1);
				Ap2 = bloc_A[i].IndexData[j+1];
				if (Ap2 == Ap1+1)				//have two consecutive lines
				{
					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));

					if (Av1_col2 != 0)
						R.negin(Av1_col2);
					if (Av2_col2 != 0)
						R.negin(Av2_col2);

					++j;

					TIMER_START_(Axpy2);
					axpy2(Av1_col1, Av2_col1,
							Av1_col2,
							Av2_col2,
							bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
							Bloc_acc[i*2],
							Bloc_acc[i*2+1]);
					TIMER_STOP_(Axpy2);
				}
				else
				{
					TIMER_START_(Axpy1);

					axpy(Av1_col1, Av2_col1,
							bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
							Ap1 % NB_ROWS_PER_MULTILINE,
							Bloc_acc[i*2],
							Bloc_acc[i*2+1]);
					TIMER_STOP_(Axpy1);
				}
			}
			else
			{
				TIMER_START_(Axpy1);

				axpy(Av1_col1, Av2_col1,
						bloc_B[Ap1 / NB_ROWS_PER_MULTILINE],
						Ap1 % NB_ROWS_PER_MULTILINE,
						Bloc_acc[i*2],
						Bloc_acc[i*2+1]);
				TIMER_STOP_(Axpy1);
			}
		}

		bool reduce = true;

		if (bloc_A[i].size() > 1) //reduce lines within the same multiline
		{
			int j=bloc_A[i].size()-2;
			Ap1 = bloc_A[i].IndexData[j];
			R.copy(Av1_col1, bloc_A[i].at_unchecked(1, j));

			if (Av1_col1 != 0)
			{
				R.negin(Av1_col1);

				for (uint32 t = 0; t < DEFAULT_BLOC_WIDTH; ++t)
					Bloc_acc[i*2][t] %= R._modulus; //Make sure product in next loop doesn't overflow

				//TODO: make this one loop
				for (uint32 t = 0; t < DEFAULT_BLOC_WIDTH; ++t)
					Bloc_acc[i*2+1][t] += (uint32) Av1_col1 * Bloc_acc[i*2][t];

				for (uint32 t = 0; t < DEFAULT_BLOC_WIDTH; ++t)
					Bloc_acc[i*2+1][t] %= R._modulus;

				reduce = false;
			}
		}

		TIMER_START_(copyDenseArrayToSparseRow);
			//cout << "write back " << i  << endl;
			copyDenseArrayToSparseRow(R, Bloc_acc[i*2], Bloc_acc[i*2+1], bloc_B.bloc_width(), bloc_B[i], reduce);
		TIMER_STOP_(copyDenseArrayToSparseRow);
	}  //for i

	//cout << "REDUCE TRIANGULAR DONE\n";
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
		dense_bloc[i] = (uint64 *) memalign(64, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseBlocModulo);
	TIMER_RESET_(axpyInTriangular);
	TIMER_RESET_(axpyInRectangular);
	TIMER_RESET_(copyDenseArrayToSparseRow);


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
				reduceBlocByTriangularBloc(R, A[j][last_bloc_idx], dense_bloc, B[j][i]);
			TIMER_STOP_(reduceBlocByTriangularBloc);
			//report << "reduceBlocByTriangularBloc DONE" << endl;

			TIMER_START_(copyDenseBlocArrayToSparseBloc);
				//copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i]);
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
	TIMER_REPORT_(copyDenseArrayToSparseRow);

	report << endl;
	TIMER_REPORT_(reduceDenseBlocModulo);
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
		//dense_bloc[i] = new uint64[D.bloc_width()];
		dense_bloc[i] = (uint64 *) memalign(64, D.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	TIMER_RESET_(reduceDenseBlocModulo);
	TIMER_RESET_(axpyInTriangular);
	TIMER_RESET_(axpyInRectangular);
	TIMER_RESET_(copyDenseArrayToSparseRow);

	//for all columns of blocs of B
	for (uint32 i = 0; i < (uint32) std::ceil((double) D.coldim() / D.bloc_width()); ++i)
	{
		//report << "Colum D|B\t" << i << endl;
		//for all rows of blocs in C
		for(uint32 j=0; j<(uint32)std::ceil((double)C.rowdim() / C.bloc_height()); ++j)
		{
			uint32 first_bloc_idx, last_bloc_idx;
			first_bloc_idx = C.FirstBlocsColumIndexes[j] / C.bloc_width();
			last_bloc_idx = C[j].size ();

			//report << "\tRow C\t" << j << "\tfirst bloc: " << first_bloc_idx << " - last: " << last_bloc_idx << endl;
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
	TIMER_REPORT_(axpyInRectangular);

	report << endl;
	TIMER_REPORT_(copyDenseBlocArrayToSparseBloc);

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);
}




#endif /* MATRIX_OPS_H_ */




