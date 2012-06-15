/*
 * matrix-op.C
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 */


#ifndef MATRIX_OP_C_
#define MATRIX_OP_C_

#include <iostream>
#include <ctime>
#include <cmath>
#include <assert.h>

#include "matrix-op.h"
#include "matrix-util.h"

template <typename Ring>
inline void MatrixOp::razArray(const Ring& R, typename Ring::Element arr[], uint32 arrSize)
{
	memset(arr, 0, arrSize*sizeof(typename Ring::Element));
}

template <typename Ring, typename Iterator>
inline void MatrixOp::copySparseVectorToDenseArray(Ring R, const Iterator& v_start, const Iterator& v_end, typename Ring::Element array[])
{
	Iterator i_v = v_start;

	while(i_v != v_end)
	{
		array[i_v->first] = i_v->second;
		//R.copy(array[i_v->first], i_v->second);	//This should be used to make sure modulo op is made
		++i_v;
	}
}

template <typename Ring, typename Vector>
inline void MatrixOp::copyDenseArrayToSparseVector(Ring& R, typename Ring::Element array[], uint32 arraySize, Vector& v)
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

inline void razArray64(uint64 arr[], uint32 arrSize)
{
	memset(arr, 0, arrSize*sizeof(uint64));
}


template <typename Iterator>
inline void copySparseVectorToDenseArray64(const Iterator& v_start, const Iterator& v_end, uint64 array[])
{
	Iterator i_v = v_start;

	while(i_v != v_end)
	{
		array[i_v->first] = (uint64)(i_v->second);
		++i_v;
	}
}

template <typename Ring, typename Vector>
inline void copyDenseArrayToSparseVector64(Ring& R, uint64 array[], uint32 arraySize, Vector& v)
{
	uint32 nb_elts = 0;
	Vector tmp;

	typename Ring::Element e;

	for (uint32 i = 0; i < arraySize; ++i)
		if(array[i] != 0)
			++nb_elts;

	tmp.reserve (nb_elts);

	for (uint32 i = 0; i < arraySize; ++i){

		//if(!R.isZero(e))
		if(array[i] % R._modulus != 0)
		{
			ModularTraits<typename Ring::Element>::reduce (e, array[i], R._modulus);
			tmp.push_back (typename Vector::value_type (i, e));
		}
	}

	std::swap(tmp, v);
}

//reduce vector x so that the first entry is equal to unity
template <typename Ring, typename Vector>
void MatrixOp::normalizeRow(const Ring& R, Vector& x)
{
	typename Vector::iterator i_x = x.begin ();
	typename Ring::Element inv;

	Context<Ring> ctx (R);

	if(x.empty ())
		return;

	if(R.inv(inv, i_x->second) != true)		//should be invertible
		throw "Non Invertible Value";

	BLAS1::scal(ctx, inv, x);				//Lela implementation is buggy over non primes P! handle it
}

template <typename Vector>
inline void axpy(const uint16 a, const Vector& v, uint64 *arr, uint8 STEP_)
{
	const uint8 STEP=32;
	uint32 sz = v.size ();
	register uint32 a32 = (uint32)a;
	//register uint32 x=0;

	//for(x=0; x<sz; ++x)
		//arr[v[x].first] += a32 * v[x].second;

	uint8 xl = v.size () % STEP;
	register uint32 x=0;
	while(x<xl)
	{
		arr[v[x].first] += a32 * v[x].second;
		++x;
		//++row_it_B;
	}

	for(x=xl; x<sz; x+=STEP)
	{

#pragma loop unroll
		for(uint8 t=0; t<STEP; ++t)
			arr[v[x+t].first] += a32 * v[x+t].second;
			//R.axpyin(tmpDenseArray[(row_it_B+t)->first], 	 Cv, (row_it_B+t)->second);

		//row_it_B += 32;
	}
}

//B <- A^-1 x B
template <typename Matrix, typename Ring>
void MatrixOp::reducePivotsByPivots(Ring& R, const Matrix& A, Matrix& B)
{
	assert(A.rowdim () == B.rowdim ());
	assert(A.coldim () == A.rowdim ());

	typename Matrix::ConstRowIterator i_A;
	typename Matrix::Row::const_iterator row_it_A, row_it_A_end;

	typename Matrix::RowIterator i_B, B_ap;
	typename Matrix::Row::const_iterator row_it_B_end;

	uint32 B_coldim = B.coldim ();

	typename Ring::Element tmpDenseArray[B_coldim];
	MatrixOp::razArray(R, tmpDenseArray, B_coldim);


	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(AxpyTimer);

#ifdef SHOW_PROGRESS
	uint32 i=A.rowdim ();
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif

	i_A = A.rowEnd ();
	i_B = B.rowEnd ();

	//skip last row
	--i_A;
	--i_B;
	if(A.rowdim () == 1)
		return;


	do {
		--i_A;
		--i_B;

#ifdef SHOW_PROGRESS
		--i;
                report << "                                                                    \r";
                report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			MatrixOp::razArray(R, tmpDenseArray, B_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			MatrixOp::copySparseVectorToDenseArray(R, i_B->begin (), i_B->end (), tmpDenseArray);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		row_it_A = i_A->begin ();
		row_it_A_end = i_A->end ();
		++row_it_A;	//skip first element

		register typename Ring::Element Av;
		register uint32 Ap;
		register typename Matrix::Row::const_iterator row_it_B;
		while(row_it_A != row_it_A_end)
		{
			Ap = row_it_A->first;
			R.copy(Av, row_it_A->second);
			R.negin(Av);

			// B[i] <- B[i] - Av * B[Ap]
			row_it_B = B[Ap].begin ();
			row_it_B_end = B[Ap].end ();

#ifdef PRAGMA_UNROLL
			TIMER_START_(AxpyTimer);
			register uint32 x=0, xl = B[Ap].size () % 32;
			while(x<xl)
			{
				R.axpyin(tmpDenseArray[row_it_B->first], Av, row_it_B->second);
				++x;
				++row_it_B;
			}

			for(x=xl; x<B[Ap].size (); x+=32)
			{

//#pragma loop unroll
				for(char t=0; t<32; ++t)
				{
					R.axpyin(tmpDenseArray[(row_it_B+t)->first], 	 Av, (row_it_B+t)->second);
				}

				row_it_B += 32;
			}
			TIMER_STOP_(AxpyTimer);

#else

			TIMER_START_(AxpyTimer);
			while(row_it_B != row_it_B_end)
			{
				R.axpyin(tmpDenseArray[row_it_B->first], Av, row_it_B->second);
				++row_it_B;
			}
			TIMER_STOP_(AxpyTimer);
#endif
			++row_it_A;
		}

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			MatrixOp::copyDenseArrayToSparseVector(R, tmpDenseArray, B_coldim, *i_B);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);
	} while(i_A != A.rowBegin ());

#ifdef SHOW_PROGRESS
        report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
	TIMER_REPORT_(AxpyTimer);
}

//B <- A^-1 x B
template <typename Matrix>
void MatrixOp::reducePivotsByPivots(Modular<uint16>& R, const Matrix& A, Matrix& B, bool t, uint8 STEP)
{
	assert(A.rowdim () == B.rowdim ());
	assert(A.coldim () == A.rowdim ());

	typedef Modular<uint16> Ring;
	
	typename Matrix::ConstRowIterator i_A;
	typename Matrix::Row::const_iterator row_it_A, row_it_A_end;

	typename Matrix::RowIterator i_B, B_ap;
	typename Matrix::Row::const_iterator row_it_B_end;

	uint32 B_coldim = B.coldim ();
	uint64 tmpDenseArray[B_coldim];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(AxpyTimer);
	TIMER_DECLARE_(AxpyOuterTimer);

#ifdef SHOW_PROGRESS
	uint32 i=A.rowdim ();
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16>" << std::endl;
#endif
	i_A = A.rowEnd ();
	i_B = B.rowEnd ();

	//skip last row
	--i_A;
	--i_B;

	if(A.rowdim () == 1)
		return;

	do {
		--i_A;
		--i_B;

#ifdef SHOW_PROGRESS
		--i;
                report << "                                                                    \r";
                report << "\t" << i << std::ends;
#endif

		TIMER_START_(RazArrayTimer);
			 razArray64(tmpDenseArray, B_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copySparseVectorToDenseArray64(i_B->begin (), i_B->end (), tmpDenseArray);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		row_it_A = i_A->begin ();
		row_it_A_end = i_A->end ();
		++row_it_A;	//skip first element

		TIMER_START_(AxpyOuterTimer);
		uint32 Ap;
		typename Ring::Element Av;
		//register uint32 Av32;
		
		while(row_it_A != row_it_A_end)
		{
			R.copy(Av, row_it_A->second);
			Ap = row_it_A->first;
			R.negin(Av);
			//Av32 = Av;

			// B[i] <- B[i] - Av * B[Ap]
			//rowB = &(B[Ap]);

			TIMER_START_(AxpyTimer);
			/*register uint32 x=0, sz = B[Ap].size ();
			for(x=0; x<sz; ++x)
				tmpDenseArray[(*rowB)[x].first] += Av32 * (*rowB)[x].second;*/

			axpy(Av, B[Ap], tmpDenseArray, STEP);
			
			TIMER_STOP_(AxpyTimer);

			++row_it_A;
		}
		TIMER_STOP_(AxpyOuterTimer);

		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArrayToSparseVector64(R, tmpDenseArray, B_coldim, *i_B);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);
		
	} while(i_A != A.rowBegin ());

#ifdef SHOW_PROGRESS
        report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
	TIMER_REPORT_(AxpyTimer);
	TIMER_REPORT_(AxpyOuterTimer);
}


// D <- D - CB
template <typename Matrix, typename Ring>
void MatrixOp::reduceNonPivotsByPivots(Ring& R, const Matrix& C, const Matrix& B, Matrix& D)
{
	assert(C.rowdim () == D.rowdim ());
	assert(C.coldim () == B.rowdim ());
	assert(B.coldim () == D.coldim ());

	typename Matrix::ConstRowIterator i_C;
	typename Matrix::RowIterator i_D;
	typename Matrix::Row::const_iterator row_it_C, row_it_C_end;
	typename Matrix::Row::const_iterator row_it_B, row_it_B_end;

	uint32 D_coldim = D.coldim ();

	typename Ring::Element tmpDenseArray[D_coldim];
	MatrixOp::razArray(R, tmpDenseArray, D_coldim);

#ifdef SHOW_PROGRESS
	uint32 i=0;
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif

	for(i_C = C.rowBegin (), i_D = D.rowBegin (); i_C != C.rowEnd (); ++i_C, ++i_D){

#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << i << std::ends;
		++i;
#endif

		MatrixOp::razArray(R, tmpDenseArray, D_coldim);
		MatrixOp::copySparseVectorToDenseArray(R, i_D->begin (), i_D->end (), tmpDenseArray);

		row_it_C = i_C->begin ();
		row_it_C_end = i_C->end ();

		register typename Ring::Element Cv;
		register uint32 Cp;

		while(row_it_C != row_it_C_end)
		{
			R.copy(Cv, row_it_C->second);
			R.negin(Cv);

			Cp = row_it_C->first;

			// D[i] <- D[i] - Cv * B[Cp]
			row_it_B = B[Cp].begin ();
			row_it_B_end = B[Cp].end ();


#ifdef PRAGMA_UNROLL
			//TIMER_START_(AxpyTimer);
			register uint32 x=0, xl = B[Cp].size () % 32;
			while(x<xl)
			{
				R.axpyin(tmpDenseArray[row_it_B->first], Cv, row_it_B->second);
				++x;
				++row_it_B;
			}

			for(x=xl; x<B[Cp].size (); x+=32)
			{

//#pragma loop unroll
				for(char t=0; t<32; ++t)
					R.axpyin(tmpDenseArray[(row_it_B+t)->first], 	 Cv, (row_it_B+t)->second);

				row_it_B += 32;
			}
			//TIMER_STOP_(AxpyTimer);

#else
			while(row_it_B != row_it_B_end)
			{
				R.axpyin(tmpDenseArray[row_it_B->first], Cv, row_it_B->second);
				++row_it_B;
			}

#endif

			++row_it_C;
		}

		MatrixOp::copyDenseArrayToSparseVector(R, tmpDenseArray, D_coldim, *i_D);

	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif

}


// D <- D - CB
template <typename Matrix>
void MatrixOp::reduceNonPivotsByPivots(Modular<uint16>& R, const Matrix& C, const Matrix& B, Matrix& D, bool t)
{
	assert(C.rowdim () == D.rowdim ());
	assert(C.coldim () == B.rowdim ());
	assert(B.coldim () == D.coldim ());

	typedef Modular<uint16> Ring;

	typename Matrix::ConstRowIterator i_C;
	typename Matrix::RowIterator i_D;
	typename Matrix::Row::const_iterator row_it_C, row_it_C_end;
	typename Matrix::Row::const_iterator row_it_B, row_it_B_end;

	typename Matrix::ConstRow *rowB;

	uint32 D_coldim = D.coldim ();
	uint64 tmpDenseArray[D_coldim];

#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	uint32 i=0;
	report << "In spec Modular<uint16>" << std::endl;
#endif

	for(i_C = C.rowBegin (), i_D = D.rowBegin (); i_C != C.rowEnd (); ++i_C, ++i_D){

#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << i << std::ends;
		++i;
#endif

		razArray64(tmpDenseArray, D_coldim);
		copySparseVectorToDenseArray64(i_D->begin (), i_D->end (), tmpDenseArray);

		row_it_C = i_C->begin ();
		row_it_C_end = i_C->end ();

		register typename Ring::Element Cv;
		register uint32 Cv32;
		uint32 Cp;

		while(row_it_C != row_it_C_end)
		{
			Cp = row_it_C->first;
			R.copy(Cv, row_it_C->second);
			R.negin(Cv);
			Cv32 = Cv;


			// D[i] <- D[i] - Cv * B[Cp]
			rowB = &(B[Cp]);

			//TIMER_START_(AxpyTimer);
			register uint32 x=0, sz = B[Cp].size ();
			for(x=0; x<sz; ++x)
				tmpDenseArray[(*rowB)[x].first] += Cv32 * (*rowB)[x].second;

			//TIMER_STOP_(AxpyTimer);

			++row_it_C;
		}

		copyDenseArrayToSparseVector64(R, tmpDenseArray, D_coldim, *i_D);
	}

#ifdef SHOW_PROGRESS
	report << "\r                                                                    \n";
#endif

}

template <typename Ring, typename Matrix>
size_t  MatrixOp::echelonize(const Ring& R, Matrix& A)
{
	uint32 coldim = A.coldim ();
	uint32 npiv = 0;
	uint32 N = A.rowdim();
	uint32 piv[N];

	typename Matrix::Row::const_iterator it, end_iterator;
	typename Matrix::RowIterator i_A;

	typename Ring::Element tmpDenseArray[coldim];
	MatrixOp::razArray(R, tmpDenseArray, coldim);

#ifdef SHOW_PROGRESS
		std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif
	//make sure first row is not null
	for (uint i=0; i < N; ++i)
	{
		if(!A[i].empty ())
		{
			if(i!= 0)
				std::swap(A[i], A[0]);

			normalizeRow(R, A[0]);
			piv[npiv] = 0;
			npiv++;

			break;
		}
	}

	uint32 curr_head;
	uint32 j, i=1;

	i_A = A.rowBegin (); ++i_A;
	for(; i_A != A.rowEnd (); ++i_A, ++i)
	{
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << std::ends;
#endif
		curr_head = i_A->front ().first;
		register typename Ring::Element h_a = R.one ();

		MatrixOp::razArray(R, tmpDenseArray, coldim);
		MatrixOp::copySparseVectorToDenseArray(R, i_A->begin (), i_A->end (), tmpDenseArray);

		for (j=0; j < npiv; ++j)
		{
			if(A[piv[j]].front ().first < curr_head)
				continue;

			it = A[piv[j]].begin ();
			end_iterator =  A[piv[j]].end ();

			if(it != end_iterator)		//nothing to reduce with this pivot, corresponding entry is 0
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

			register uint32 x=0;
			unsigned char xl = A[piv[j]].size () % 32;
			while(x<xl)
			{
				R.axpyin(tmpDenseArray[it->first], 	h_a, 	it->second);
				++x;
				++it;
			}

			for(x=xl; x<A[piv[j]].size (); x+=32)
			{

//#pragma loop unroll
				for(char t=0; t<32; ++t)
				{
					R.axpyin(tmpDenseArray[(it+t)->first], 	h_a, 	(it+t)->second);
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

		MatrixOp::copyDenseArrayToSparseVector(R, tmpDenseArray, coldim, *i_A);
		normalizeRow(R, *i_A);

		if(!i_A->empty())
		{
			piv[npiv] = i;
			npiv++;
		}
	}
#ifdef SHOW_PROGRESS
                report << "\r                                                                    \n";
#endif

	return npiv;
}

template <typename Matrix>
size_t  MatrixOp::echelonize(const Modular<uint16>& R, Matrix& A, bool t)
{
	uint32 coldim = A.coldim ();
	uint32 npiv = 0;
	uint32 N = A.rowdim();
	uint32 piv[N];

#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16>" << std::endl;
#endif
	typedef Modular<uint16> Ring;

	typename Matrix::Row::const_iterator it, end_iterator;
	typename Matrix::RowIterator i_A;
	typename Matrix::Row *rowA;

	uint64 tmpDenseArray[coldim];

	//make sure first row is not null
	for (uint i=0; i < N; ++i)
	{
		if(!A[i].empty ())
		{
			if(i!= 0)
				std::swap(A[i], A[0]);

			normalizeRow(R, A[0]);
			piv[npiv] = 0;
			npiv++;

			break;
		}
	}

	uint32 curr_head;
	uint32 j, i=1;

	i_A = A.rowBegin (); ++i_A;
	for(; i_A != A.rowEnd (); ++i_A, ++i)
	{
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << std::ends;
#endif
		curr_head = i_A->front ().first;

		razArray64(tmpDenseArray, coldim);
		copySparseVectorToDenseArray64(i_A->begin (), i_A->end (), tmpDenseArray);

		register typename Ring::Element h_a = R.one ();
		register uint32 h_a32;

		for (j=0; j < npiv; ++j)
		{
			if(A[piv[j]].front ().first < curr_head)
				continue;

			it = A[piv[j]].begin ();
			end_iterator =  A[piv[j]].end ();

			if(it != end_iterator)		//nothing to reduce with this pivot, corresponding entry is 0
			{
				if(tmpDenseArray[it->first] % R._modulus == 0)
					continue;
				else
				{
					h_a = tmpDenseArray[it->first] % R._modulus;
					R.negin(h_a);
				}
			}

			h_a32 = h_a;
			rowA = &(A[piv[j]]);
			register uint32 x=0, sz = A[piv[j]].size ();
			for(x=0; x<sz; ++x)
				tmpDenseArray[(*rowA)[x].first] += h_a32 * (*rowA)[x].second;

		}

		copyDenseArrayToSparseVector64(R, tmpDenseArray, coldim, *i_A);
		normalizeRow(R, *i_A);

		if(!i_A->empty())
		{
			piv[npiv] = i;
			npiv++;
		}
	}
#ifdef SHOW_PROGRESS
                report << "\r                                                                    \n";
#endif

	return npiv;
}


template <typename Matrix, typename Ring>
void MatrixOp::normalizeRows(Ring& R, Matrix& A)
{
	Context<Ring> ctx (R);

	typename Matrix::Row::iterator it;
	typename Matrix::RowIterator i_A;
	typename Ring::Element inv;

	for(i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A){
		it = i_A->begin ();

		R.inv(inv, it->second);
		BLAS1::scal(ctx, inv, *i_A);
	}
}


#endif	//MATRIX_OP_
