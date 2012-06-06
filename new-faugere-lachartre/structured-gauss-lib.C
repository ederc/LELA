#ifndef STRUCTURED_GAUSS_LIB_C_
#define STRUCTURED_GAUSS_LIB_C_

/*
 * structured-gauss.C
 *
 *  Created on: 24 mai 2012
 *      Author: martani
 *
 *  Inplace Gaussian elimination
 */

#include "structured-gauss-lib.h"
//#include "matrix-util.h"

using namespace LELA;
using namespace std;

template <typename Vector>
inline bool StructuredGauss::isRowNull(Vector& v)
{
	return v.empty ();
}

template <typename Vector>
inline void StructuredGauss::swapRows(Vector&v, Vector& w)
{
	std::swap(v, w);
}

//reduce row so that the first entry is equal to unity
template <typename Ring, typename Vector>
void StructuredGauss::normalizeRow(const Ring& R, Vector& x)
{
	typename Vector::iterator i_x = x.begin ();
	typename Ring::Element inv;
	Context<Ring> ctx (R);

	if(x.empty ())
		return;

	assert(R.inv(inv, i_x->second) == true);		//should be invertible
	BLAS1::scal(ctx, inv, x);
}

template <typename Array, typename Matrix>
void StructuredGauss::sortRows(Matrix& A, Array pivots)
{
	uint32 row_index;
	uint t, real_j;

	uint32 perm[A.rowdim ()];
	for (uint i = 0; i < A.rowdim (); ++i)
			perm[i] = i;

	for (uint i=0, j=0; i < A.coldim (); ++i)
	{
		if(pivots[i] != 0)
		{
			row_index = pivots[i]-1;
			std::swap(A[perm[row_index]], A[j]);

			//find [who:which row] is in index j, we don not update the entry at index j
			//but the entry of the row that is in index j now
			//TODO don't search for cycles, store the original index on the go
			real_j = perm[j];
			while(j != perm[real_j])
				real_j = perm[real_j];

			t = perm[row_index];
			perm[row_index] = j;
			perm[real_j] = t;

			j++;
		}
	}
}

template <typename Ring>
inline void StructuredGauss::razArray_in(const Ring& R, typename Ring::Element arr[], uint32 arrSize)
{
	memset(arr, 0, arrSize*sizeof(typename Ring::Element));
}

//no bounds checking on the array - the caller must check for the bounds
template <typename Ring, typename Vector>
inline void StructuredGauss::copySparseVectorToDenseArray_in(Ring R, const Vector& v, typename Ring::Element array[])
{
	typename Vector::const_iterator i_v = v.begin (), i_end = v.end ();

	while(i_v != i_end)
	{
		R.copy(array[i_v->first], i_v->second);
		++i_v;
	}
}

template <typename Ring, typename Iterator>
inline void StructuredGauss::copySparseVectorToDenseArray_in(Ring R, const Iterator& v_start, const Iterator& v_end, typename Ring::Element array[])
{
	Iterator i_v = v_start;

	while(i_v != v_end)
	{
		R.copy(array[i_v->first], i_v->second);
		++i_v;
	}
}

template <typename Ring, typename Vector>
inline void StructuredGauss::copyDenseArrayToSparseVector_in(Ring& R, typename Ring::Element array[], uint32 arraySize, Vector& v)
{
	uint32 nb_elts = 0;
	Vector tmp;

	for (uint32 i = 0; i < arraySize; ++i)
		if(!R.isZero(array[i]))
			++nb_elts;

	tmp.reserve (nb_elts);

	for (uint32 i = 0; i < arraySize; ++i){
		if(!R.isZero(array[i]))
			tmp.push_back (typename Vector::value_type (i, array[i]));
	}

	std::swap(tmp, v);
}

template<typename Context, typename Element, typename Vector>
inline typename Vector::value_type::first_type StructuredGauss::head_generic(Context& ctx, Element& a, const Vector& v)
{
	return BLAS1::head(ctx, a, v);
}

template <typename Ring, typename Matrix>
size_t  StructuredGauss::echelonize_reduced(const Ring& R, Matrix& A)
{

#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif

	Context<Ring> ctx(R);

	uint32 coldim = A.coldim ();
	uint32 npiv = 0;		//number of pivots found so far
	uint32 N = A.rowdim();
	uint32 piv[N];
	uint32 pivot_rows_index[coldim];
	
	typename Matrix::Row::const_iterator it, end_iterator;
	typename Matrix::RowIterator i_A;
	typename Ring::Element a;
	typename Ring::Element tmpDenseArray[coldim];
	
	for (uint i = 0; i < coldim; ++i) {
		pivot_rows_index[i] = 0;
	}

	//make sure first row is not null
	for (uint i=0; i < N; ++i)
	{
		if(!isRowNull(A[i]))
		{
			if(i!= 0)
				swapRows(A[i], A[0]);

			normalizeRow(R, A[0]);
			piv[npiv] = 0;
			npiv++;

			pivot_rows_index[head_generic(ctx, a, A[0])] = 0+1;	//indexes are +1 shifted, 0 means no row is at this index
			break;
		}
	}

	commentator.start("First reduction", __FUNCTION__);

	uint32 new_entry, curr_head;
	uint32 j, i=1;

	i_A = A.rowBegin (); ++i_A;
	for(; i_A != A.rowEnd (); ++i_A, ++i)
	{

#ifdef SHOW_PROGRESS
		report << "                                                                    \r";
		report << "\t" << npiv;
#endif
		curr_head = head_generic(ctx, a, *i_A);

		razArray_in(R, tmpDenseArray, coldim);
		copySparseVectorToDenseArray_in(R, i_A->begin (), i_A->end (), tmpDenseArray);
		
		for (j=0; j < npiv; ++j)
		{
			register typename Ring::Element h_a = R.zero ();

			if(head_generic(ctx, h_a, A[piv[j]]) < curr_head)
				continue;

			it = A[piv[j]].begin ();
			end_iterator =  A[piv[j]].end ();

			if(it != end_iterator)
			{
				if(R.isZero(tmpDenseArray[it->first]))
					continue;
				else
				{
					R.copy(h_a, tmpDenseArray[it->first]);
					R.negin(h_a);
				}
			}

#ifdef PRAGMA_UNROLL
			register uint32 x=0, xl = A[piv[j]].size () % 32;
			while(x<xl)
			{
				R.axpyin(tmpDenseArray[it->first], h_a, it->second);
				++x;
				++it;
			}

			for(x=xl; x<A[piv[j]].size (); x+=32)
			{

//#pragma unroll(32)
				for(char t=0; t<32; ++t)
				{
					R.axpyin(tmpDenseArray[(it+t)->first], h_a, (it+t)->second);
				}

				it += 32;
			}
#else

			while (it != end_iterator)
			{
				R.axpyin(tmpDenseArray[it->first], h_a, it->second);
				++it;
			}
#endif

		}

		copyDenseArrayToSparseVector_in(R, tmpDenseArray, coldim, *i_A);
		normalizeRow(R, *i_A);
		
		new_entry = head_generic(ctx, a, *i_A);
		if(!isRowNull( *i_A ) && pivot_rows_index[new_entry]==0)
		{
			piv[npiv] = i;
			npiv++;
			pivot_rows_index[new_entry] = i+1;	//new_entry must not be -1 since the row is not empty
		}

	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    " << std::endl;
#endif

	commentator.stop(MSG_DONE);
	commentator.start("Second reduction", __FUNCTION__);
	
	for (uint i = 0; i < npiv; ++i)
	{

#ifdef SHOW_PROGRESS
		report << "                                                                    \r";
		report << "\t" << i;
#endif

		razArray_in(R, tmpDenseArray, coldim);
		copySparseVectorToDenseArray_in(R, A[piv[i]], tmpDenseArray);

		//there is no reason to start from the begining, because this row had already been eliminated by all
		//the pivots discovered before it in the previous loop
		for (uint j = i+1; j < npiv; ++j) //npiv; ++j) { //TODO change this
		{
			it = A[piv[j]].begin ();			//this vector is garanteed not to be null since it's in the pivot list
			end_iterator =  A[piv[j]].end ();

			register typename Ring::Element h_a = R.zero ();

			if(R.isZero(tmpDenseArray[it->first]))
				continue;
			else
			{
				R.copy(h_a, tmpDenseArray[it->first]);
				R.negin(h_a);
			}
#ifdef PRAGMA_UNROLL
			register uint32 x=0, xl = A[piv[j]].size () % 32;

			while(x<xl)
			{
				R.axpyin(tmpDenseArray[it->first], h_a, it->second);
				++x;
				++it;
			}

			for(x=xl; x<A[piv[j]].size (); x+=32)
			{

//#pragma unroll(32)
				for(char t=0; t<32; ++t)
				{
					R.axpyin(tmpDenseArray[(it+t)->first], h_a, (it+t)->second);
				}

				it += 32;
			}

#else
			while (it != end_iterator)
			{
				R.axpyin(tmpDenseArray[it->first], h_a, it->second);
				++it;
			}
#endif

		}

		copyDenseArrayToSparseVector_in(R, tmpDenseArray, coldim, A[piv[i]]);

	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    " << std::endl;
#endif

	commentator.stop(MSG_DONE);

	commentator.start("Sorting pivots", __FUNCTION__);
		sortRows(A, pivot_rows_index);
	commentator.stop(MSG_DONE);

	return npiv;
}


#endif //STRUCTURED_GAUSS_LIB_
