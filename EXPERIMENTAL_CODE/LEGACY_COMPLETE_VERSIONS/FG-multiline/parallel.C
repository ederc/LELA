/*
 * parallel.C
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 *  Created on: 30 juin 2012
 *      Author: martani
 */

#ifndef PARALLEL_MULTILINE_
#define PARALLEL_MULTILINE_

#include <vector>
#include <pthread.h>

#include "FG-types.h"
#include "matrix-util-m.C"
#include "matrix-op-m.C"
#include "queue.h"

#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"


using namespace LELA;


namespace NS{

	
template <typename Element>
class echelonize_params_t {
public:
	const Modular<Element> *R;
	SparseMultilineMatrix<Element> *A;
};

typedef struct {
	uint32 row_index;
	uint32 last_pivot_reduced_by;
} waiting_row_t;


uint32 *reduced_pivots;
uint32 last_reduced_pivot;
SyncQueue<uint32> queue_ready_rows;
SyncQueue<waiting_row_t> queue_not_fully_reduced_pivots;



template <typename Element>
void echelonize_parallel(const Modular<Element>& R, SparseMultilineMatrix<Element>& A)
{
	echelonize_params_t<Element> params;
	params.R = &R;
	params.A = &A;

	uint32 N = A.rowdim()/2 + A.rowdim()%2;
	reduced_pivots = new uint32[N];
	for(uint32 i=0; i<N; ++i)
		reduced_pivots[i] = 0;
}

static void* echelonize_in(void* params)
{
	uint32 coldim = params->A->coldim ();
	uint32 npiv = 0;
	uint32 npiv_real = 0;
	
	uint32 N = params->A->rowdim()/2 + params->A->rowdim()%2;
	uint32 local_last_reduced_pivot;
	
	typedef SparseMultilineMatrix<uint16> Matrix;
	
#ifdef SHOW_PROGRESS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16>" << std::endl;
#endif
	typedef Modular<uint16> Ring;
	//typename Matrix::RowIterator i_A;
	typename Matrix::Row *rowA;

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
	TIMER_DECLARE_(ExchangeTimer);
	TIMER_DECLARE_(GettingVxColx);
	TIMER_DECLARE_(ReduceBySameMultilineTimer);

	TIMER_START_(normalizeVectorTimer);
		normalize_multiline(R, A[0]);
	TIMER_STOP_(normalizeVectorTimer);

	waiting_row_t waiting_row;
	uint32 curr_row;
	
	uint32 start_reduce_idx, end_reduce_idx;
	
	//for(i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++i)
	//LOCK
	while(last_reduced_pivot < N)
	{
		local_last_reduced_pivot = last_reduced_pivot;
		
		if(queue_not_fully_reduced_pivots.dequeue(waiting_row))
		{
			start_reduce_idx = waiting_row.last_pivot_reduced_by;
			curr_row = waiting_row.row_index;
		}
		else if(queue_ready_rows.dequeue(curr_row))
		{
			start_reduce_idx = 0;
		}
		else
		{
			;
		}
		
	//UNLOCK
		
		assert(start_reduce_idx < local_last_reduced_pivot);
#ifdef SHOW_PROGRESS
                report << "                                                                    \r";
                report << "\t" << npiv << std::ends;
#endif

        TIMER_START_(RazArrayTimer);
        	razArray64(tmpDenseArray1, coldim);
        	razArray64(tmpDenseArray2, coldim);
	TIMER_STOP_(RazArrayTimer);

	TIMER_START_(CopySparseVectorToDenseArrayTimer);
		copyMultilineToDenseArrays64(params->A[curr_row], tmpDenseArray1, tmpDenseArray2);
	TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

        typename Ring::Element h_a1 = params->R->one (), h_a2 = params->R->one ();

		uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
		uint32 tmp=0;

		for(uint32 j=start_reduce_idx; j<local_last_reduced_pivot; ++j)
		{
			rowA = &((*params->A)[piv[j]]);
			TIMER_START_(HeadVectorTimer);
				head_line1 = head(*rowA, 0, h_a1, head_line1_idx);
				head_line2 = head(*rowA, 1, h_a2, head_line2_idx);
			TIMER_STOP_(HeadVectorTimer);

			if(head_line1 != -1 && head_line1 == head_line2)
				throw "Wrong Mutiline format";
			
			TIMER_START_(ExchangeTimer);
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
				//v = &(A[piv[j]]);
				uint16 tmp1=0;

				for(uint32 x=0; x<rowA->size (); ++x)
				{
					tmp1 = rowA->ValuesData[0+x*2];
					rowA->ValuesData[0+x*2] = rowA->ValuesData[1+x*2];
					rowA->ValuesData[1+x*2] = tmp1;
				}
				
				cout << "exchanging values\n";
			}
			TIMER_STOP_(ExchangeTimer);
			
			TIMER_START_(GettingVxColx);
			if(head_line1 != -1)
			{
				R.init(v1col1, tmpDenseArray1[head_line1]);
				R.init(v2col1, tmpDenseArray2[head_line1]);

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
				R.init(v1col2, tmpDenseArray1[head_line2]);
				R.init(v2col2, tmpDenseArray2[head_line2]);
			}
			else
			{
				v1col2 = 0;
				v2col2 = 0;
			}

			//TODO: check if buggy! if line is empty for example
			uint16 val = A[piv[j]].at_unchecked(0, head_line2_idx);

			tmp = v1col2 + (uint32)v1col1*val;
			R.init(v1col2, tmp);
			tmp = v2col2 + (uint32)v2col1*val;
			R.init(v2col2, tmp);

			if (v1col2 != 0)
				R.negin(v1col2);
			if (v2col2 != 0)
				R.negin(v2col2);
			TIMER_STOP_(GettingVxColx);
			
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
		normalize_array64(*params->R, tmpDenseArray1, coldim);
		normalize_array64(*params->R, tmpDenseArray2, coldim);
	TIMER_STOP_(normalizeArrayTimer);

	TIMER_START_(HeadArrayTimer);
		head_line1 = head(*params->R, tmpDenseArray1, coldim, h_a1);
		head_line2 = head(*params->R, tmpDenseArray2, coldim, h_a2);
        TIMER_STOP_(HeadArrayTimer);

	//assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

	
	if(curr_row == local_last_reduced_pivot+1)	//Completely reduce by others, reduce by self
	{
		TIMER_START_(ReduceBySameMultilineTimer);
		if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
		{
			assert(head_line1 >= 0 && head_line1 < coldim);
			if(tmpDenseArray2[head_line1] % params->R->_modulus != 0)
				{
					h = tmpDenseArray2[head_line1] % params->R->_modulus;
					params->R->negin(h);
				}

				for(x = head_line1; x<coldim; ++x)
				{
					tmpDenseArray2[x] += (uint32)h * tmpDenseArray1[x];
				}
		}
		TIMER_STOP_(ReduceBySameMultilineTimer);
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////

        TIMER_START_(CopyDenseArrayToSparseVectorTimer);
        	if((head_line2 != -1) && (head_line1 > head_line2))			//saves the line with the smallest column entry first
        		copyDenseArraysToMultilineVector64(*params->R, 
							tmpDenseArray2, 
							tmpDenseArray1, 
							coldim, 
							params->A[curr_row], 
							true);
        	else
        		copyDenseArraysToMultilineVector64(*params->R, 
							tmpDenseArray1, 
							tmpDenseArray2, 
							coldim, 
							params->A[curr_row],
							true);
	TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);

	TIMER_START_(normalizeVectorTimer);
		normalize_multiline(*params->R, params->A[curr_row]);
	TIMER_STOP_(normalizeVectorTimer);

	//if(!i_A->empty())
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

	TIMER_REPORT_(ExchangeTimer);
	TIMER_REPORT_(GettingVxColx);
	TIMER_REPORT_(ReduceBySameMultilineTimer);
	
	delete [] tmpDenseArray1;
	delete [] tmpDenseArray2;

	return npiv_real;
}

}

#end //PARALLEL_MULTILINE_