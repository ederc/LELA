/*
 * level1-ops.C
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */



#ifndef LEVEL1_OPS_C_
#define LEVEL1_OPS_C_

//! Note: all DoubleFlat elements are supposed to be uint64

#include "level1-ops.h"


using namespace LELA;

/**
 * Given a bloc of dense rows (an array of arrays), Zeros its memory
 */
template <typename Element>
void Level1Ops::memsetToZero(Element** arr)
{
	for(uint32 i=0; i<DEFAULT_BLOC_HEIGHT; ++i)
	{
		memset(arr[i], 0, DEFAULT_BLOC_WIDTH * sizeof(Element));
	}
}


template <typename Element>
void Level1Ops::memsetToZero(Element** arr, const uint32 nb_lines, const uint32 line_size)
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
void Level1Ops::copySparseBlocToDenseBlocArray(const Ring& R,
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
void Level1Ops::copyDenseBlocArrayToSparseBloc(const Ring& R,
		DoubleFlatElement** arr,
		SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		bool reduce_in_Ring)
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


template <typename Ring, typename DoubleFlatElement>
void Level1Ops::reduceDenseArrayModulo(const Ring& R, DoubleFlatElement* arr)
{
	typename Ring::Element e;

	for(uint32 j=0; j<DEFAULT_BLOC_WIDTH; ++j)
	{
		ModularTraits<typename Ring::Element>::reduce (e, arr[j], R._modulus);
		arr[j] = (uint64)e;
	}
}


template <typename Element, typename Index>
long Level1Ops::headMultiLineVector(const MultiLineVector<Element, Index>& v,
			const uint16 line_idx, uint16& a, uint32& head_idx)
{
	uint16 val=0;

	for (uint32 i = 0; i < v.size(); ++i)
	{
		if ((val = v.at_unchecked(line_idx, i)) != 0)
		{
			a = val;
			head_idx = i;
			return (long) v.IndexData[i];
		}
	}

	return -1;
}


template <typename Element, typename Index>
long Level1Ops::headMultiLineVectorHybrid(const MultiLineVector<Element, Index>& v,
			const uint16 line_idx, uint16& a, uint32& head_idx, const size_t SIZE_DENSE_VECTOR)
{
	uint16 val=0;

	if (v.is_sparse(SIZE_DENSE_VECTOR))
		return headMultiLineVector(v, line_idx, a, head_idx);
	else
		for (uint32 i = 0; i < SIZE_DENSE_VECTOR; ++i)
		{
			if ((val = v.at_unchecked(line_idx, i)) != 0)
			{
				a = val;
				head_idx = i;
				return (long) i;
			}
		}

	return -1;
}


//WARNING could fail on non prime rings
template <typename Ring, typename Index>
void Level1Ops::normalizeMultiLineVector(const Ring& R, MultiLineVector<typename Ring::Element, Index>& v)
{
	uint32 idx;
	uint16 h1=0, h2=0;

	if(v.empty ())
		return;

	headMultiLineVector(v, 0, h1, idx);
	headMultiLineVector(v, 1, h2, idx);

	if(h1 != 0)
		if(R.invin(h1) != true)
			throw std::logic_error ("Non Invertible Value");

	if(h2 != 0)
		if(R.invin(h2) != true)
			throw std::logic_error ("Non Invertible Value");

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


template <typename Element, typename Index, typename Element2>
void Level1Ops::copyMultiLineVectorToDenseArray(const MultiLineVector<Element, Index>& v,
		Element2 *arr1, Element2 *arr2, size_t arr_size)
{
	if(v.empty ())
		return;

	register uint32 idx;

	if (v.is_sparse(arr_size))
		for (uint32 i = 0; i < v.size(); ++i)
		{
// 			if(i >= v.IndexData._index_vector.size ())
// 				std::cout << "ERRRO\n";
			idx = v.IndexData[i];
			arr1[idx] = (Element2) v.at_unchecked(0, i);
			arr2[idx] = (Element2) v.at_unchecked(1, i);
		}
	else
		for (uint32 i = 0; i < arr_size; ++i)
		{
			arr1[i] = (Element2) v.at_unchecked(0, i);
			arr2[i] = (Element2) v.at_unchecked(1, i);
		}
}


template<typename Ring, typename DoubleFlatElement>
long Level1Ops::headDenseArray(const Ring& R, DoubleFlatElement arr[],
		const size_t size, typename Ring::Element& a)
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


//return the index of the entry of the dense array
template <typename Ring, typename DoubleFlatElement>
long Level1Ops::normalizeDenseArray(const Ring& R, DoubleFlatElement arr[], const size_t size)
{
	uint16 a;
	int h1 = headDenseArray(R, arr, size, a);

	if(h1 == -1)	//all elements are 0
		return h1;

	R.invin(a);
	for(uint32 i=h1; i<size; ++i)
	{
		arr[i] *= a;
		arr[i] %= R._modulus;
	}

	return h1;
}


template<typename Ring, typename Index, typename DoubleFlatElement>
void Level1Ops::copyDenseArraysToMultilineVector(const Ring& R,
		const DoubleFlatElement *arr1,
		const DoubleFlatElement *arr2,
		const uint32 size,
		MultiLineVector<typename Ring::Element, Index>& v)
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


template<typename Ring, typename Index, typename DoubleFlatElement>
void Level1Ops::copyDenseArraysToMultilineVectorHybrid(const Ring& R,
		const DoubleFlatElement *arr1,
		const DoubleFlatElement *arr2,
		const uint32 size,
		MultiLineVector<typename Ring::Element, Index>& v)
{
	MultiLineVector<typename Ring::Element, Index> tmp;
	typename Ring::Element e1, e2;

	v.clear();

	for (uint32 i = 0; i < size; ++i)
	{
		ModularTraits<typename Ring::Element>::reduce(e1, arr1[i], R._modulus);
		ModularTraits<typename Ring::Element>::reduce(e2, arr2[i], R._modulus);

		if ((e1 != 0) || (e2 != 0))
		{
			tmp.IndexData.push_back(i);
			tmp.ValuesData.push_back(e1);
			tmp.ValuesData.push_back(e2);
		}
	}

	if ((float)tmp.size () / (float)size < v.get_HYBRID_REPRESENTATION_THRESHOLD ())
		v.swap(tmp);
	else
	{
		Index idx = 0;
		for (uint32 j = 0; j < size; ++j)
			if (idx < tmp.size() && tmp.IndexData[idx] == j)
			{
				v.ValuesData.push_back(tmp.at_unchecked(0, idx));
				v.ValuesData.push_back(tmp.at_unchecked(1, idx));
				idx++;
			}
			else
			{
				v.ValuesData.push_back(0);
				v.ValuesData.push_back(0);
			}
	}
}





#endif /* LEVEL1_OPS_C_ */
