/*
 * echelonize_new.C
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */



namespace EchelonizeNew
{
#include "types.h"
#include "matrix-ops.h"

template <typename Element, typename Index, typename Element2>
static inline void copyMultilineToDenseArray(const MultiLineVector<Element, Index>& v,
		Element2 *arr1, Element2 *arr2)
{
	if(v.empty ())
		return;

	register uint32 idx;
	if(v.is_sparse ())
		for (uint32 i = 0; i < v.size(); ++i)
		{
			idx = v.IndexData[i];
			arr1[idx] = (Element2) v.at_unchecked(0, i);
			arr2[idx] = (Element2) v.at_unchecked(1, i);
		}
	else
		for (uint32 i = 0; i < v.size(); ++i)
		{
			arr1[i] = (Element2) v.at_unchecked(0, i);
			arr2[i] = (Element2) v.at_unchecked(1, i);
		}
}


template<typename Ring>
static void copyDenseArraysToMultilineVector(const Ring& R,
		const uint64 *arr1,
		const uint64 *arr2,
		const uint32 size,
		MultiLineVector<typename Ring::Element>& v)
{
	typename Ring::Element e1, e2;

	v.clear();

	for (uint32 i = 0; i < size; ++i)
	{
		ModularTraits<typename Ring::Element>::reduce(e1, arr1[i], R._modulus);
		ModularTraits<typename Ring::Element>::reduce(e2, arr2[i], R._modulus);

		if ((e1 != 0) || (e2 != 0))
		{
			v.IndexData.push_back(i);
			v.ValuesData.push_back(e1);
			v.ValuesData.push_back(e2);
		}
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

static inline long head_vector(const MultiLineVector<uint16>& v, const uint16 line_idx, uint16& a, uint32& head_idx)
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
static long head_dense_array(const Ring& R, uint64 arr[], const size_t size, uint16& a)
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

template <typename Ring>
static void normalize_dense_array(const Ring& R, uint64 arr[], const size_t size)
{
	uint16 a;
	int h1 = head_dense_array(R, arr, size, a);

	if(h1 == -1)
		return;

	R.invin(a);
	for(uint32 i=h1; i<size; ++i)
	{
		arr[i] *= a;
		arr[i] %= R._modulus;
	}
}

template <typename Element, typename Index>
static inline void normalize_multiline(const Modular<Element>& R, MultiLineVector<Element, Index>& v)
{
	uint32 idx;
	uint16 h1=0, h2=0;
	long a = head_vector(v, 0, h1, idx);
	a = head_vector(v, 1, h2, idx);

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


template <typename Element>
static uint32 echelonize(const Modular<Element>& R, SparseMultilineMatrix<Element>& M)
{
	uint32 coldim = M.bloc_width();
	uint32 npiv = 0;
	uint32 npiv_real = 0;
	uint32 N = M.rowdim() / NB_ROWS_PER_MULTILINE
			+ M.rowdim() % NB_ROWS_PER_MULTILINE;
	uint32 *piv;
	posix_memalign((void**)&piv, 16, N * sizeof(uint32));

	memsetToZero(&piv, 1, N);

	MultiLineVector<Element, Index> *rowA;

#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16>" << std::endl;
#endif
	typedef Modular<uint16> Ring;

	uint64 *tmpDenseArray1;
	posix_memalign((void**)&tmpDenseArray1, 16, coldim * sizeof(uint64));

	uint64 *tmpDenseArray2;
	posix_memalign((void**)&tmpDenseArray2, 16, coldim * sizeof(uint64));


	uint32 i=0;
	uint16 h=0;

	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;

	uint32 x=0;

	normalize_multiline(R, M[0]);



	for(i=0; i<N; ++i)
	{
#ifdef SHOW_PROGRESS
                //report << "                                                                    \r";
                report << "\t" << npiv << std::ends;
#endif

       	memsetToZero(&tmpDenseArray1, 1, coldim);
      	memsetToZero(&tmpDenseArray2, 1, coldim);

      	copyMultilineToDenseArray(M[i], tmpDenseArray1, tmpDenseArray2);
        typename Ring::Element h_a1 = R.one (), h_a2 = R.one ();

		uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
		uint32 tmp=0;

		for(uint32 j=0; j<npiv; ++j)
		{
			rowA = &(M[piv[j]]);
			if(rowA->empty())
			{
				Permutation[i*2][piv[j]*2] = 0;
				Permutation[i*2+1][piv[j]*2] = 0;

				Permutation[i*2][piv[j]*2+1] = 0;
				Permutation[i*2+1][piv[j]*2+1] = 0;
				continue;
			}

			head_line1 = head_vector(*rowA, 0, h_a1, head_line1_idx);
			head_line2 = head_vector(*rowA, 1, h_a2, head_line2_idx);

			if(head_line1 != -1 && head_line1 == head_line2)
				throw "Wrong Mutiline format";
			if(head_line1 > head_line2)	//makes the row with the smallest column entry first in the multiline
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

			Permutation[i*2][piv[j]*2] = v1col1;
			Permutation[i*2+1][piv[j]*2] = v2col1;

			Permutation[i*2][piv[j]*2+1] = v1col2;
			Permutation[i*2+1][piv[j]*2+1] = v2col2;

			if((v1col1 == 0 && v2col1 == 0) && (v1col2 != 0 || v2col2 != 0))
			{
				SparseScalMulSub__one_row__vect_array(
						v1col2,
						v2col2,
						*rowA,
						1,	//reduce by second line only
						tmpDenseArray1,
						tmpDenseArray2);
			}
			else if((v1col2 == 0 && v2col2 == 0) && (v1col1 != 0 || v2col1 != 0))
			{
				SparseScalMulSub__one_row__vect_array(
						v1col2,
						v2col2,
						*rowA,
						0,	//reduce by first line only
						tmpDenseArray1,
						tmpDenseArray2);
			}
			else
			{
				SparseScalMulSub__two_rows__vect_array(
						v1col1,
						v2col1,
						v1col2,
						v2col2,
						*rowA,
						tmpDenseArray1,
						tmpDenseArray2);
			}
		}

		normalize_dense_array(R, tmpDenseArray1, coldim);	//TODO combine this so it returns the position of entry
		normalize_dense_array(R, tmpDenseArray2, coldim);

		head_line1 = head_dense_array(R, tmpDenseArray1, coldim, h_a1);
		head_line2 = head_dense_array(R, tmpDenseArray2, coldim, h_a2);


        //assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

        //reduce by same multiline
        if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
        {

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

        head_line2 = head_dense_array(R, tmpDenseArray2, coldim, h_a2);		//only this can change

	if((head_line2 != -1) && (head_line1 > head_line2))			//saves the line with the smallest column entry first
	{
		copyDenseArraysToMultilineVector(R, tmpDenseArray2, tmpDenseArray1, coldim, M[i]);
	}
	else
	{
		copyDenseArraysToMultilineVector(R, tmpDenseArray1, tmpDenseArray2, coldim, M[i]);
	}

	normalize_multiline(R, M[i]);

	if(!M[i].empty ())
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

	return npiv_real;
}



}
