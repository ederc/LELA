/*
 * level3Parallel_echelon.C
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 *  Created on: 30 juil. 2012
 *      Author: martani (UPMC University Paris 06 / INRIA)
 *
 * ------------------------------
 * Performs parallel structured Gaussian elimination
 */


#ifndef ECHELON_C_
#define ECHELON_C_

#include <vector>
#include <sys/types.h>
#include <stdio.h>

#include "level3Parallel_echelon.h"
#include "consts-macros.h"

uint32 __echelonize_global_last_piv; //the greatest pivot available. All rows before this are already reduced
uint32 __echelonize_global_next_row_to_reduce; //the next row to reduce

#ifdef USE_MUTEX
	static pthread_mutex_t __echelonize_mutex_lock;
	static pthread_mutex_t __echelonize_waiting_list_mutex_lock;
#else
	static pthread_spinlock_t __echelonize_spinlock_lock;
	static pthread_spinlock_t __echelonize_waiting_list_spinlock_lock;
#endif


//TODO implement this as a priority heap
std::vector<Level3ParallelEchelon::waiting_row_t> waiting_list;

struct waiting_row_t_Cmp {
    bool operator() (const Level3ParallelEchelon::waiting_row_t &lhs, const Level3ParallelEchelon::waiting_row_t &rhs)
    {
        return lhs.row_idx > rhs.row_idx;
    }
};
	
bool Level3ParallelEchelon::getSmallestWaitingRow(waiting_row_t* elt)
{
	//LOCK(__echelonize);
		if(waiting_list.empty ())
		{
			//UNLOCK(__echelonize);
			return false;
		}
		
		//TODO: use a priority heap instead
		std::sort(waiting_list.begin (), waiting_list.end (), waiting_row_t_Cmp ());

		elt->row_idx = waiting_list.back ().row_idx;
		elt->last_pivot_reduced_by = waiting_list.back ().last_pivot_reduced_by;
		
		//take it out
		waiting_list.pop_back ();
	//UNLOCK(__echelonize);
	
	return true;
}

void Level3ParallelEchelon::pushRowToWaitingList(uint32 row_idx, uint32 last_pivot_reduced_by)
{
	waiting_row_t tmp;
	tmp.row_idx = row_idx;
	tmp.last_pivot_reduced_by = last_pivot_reduced_by;
	
	waiting_list.push_back (tmp);
}

template<typename Index>
uint32 Level3ParallelEchelon::echelonize__Parallel(const Modular<uint16>& R,
		SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
		SparseMultilineMatrix<uint16>& outMatrix, 
		bool destruct_in_matrix,
		int NB_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3ParallelEchelon::echelonize__Parallel] NB THREADS " << NB_THREADS << std::endl;

	typedef uint16 Element;
	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;
	outMatrix = SparseMultilineMatrix<uint16> (inMatrix.rowdim (), inMatrix.coldim ());

	check_equal_or_raise_exception(inMatrix.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(inMatrix.isFilledWithEmptyBlocs (), true);
	
	//copy bloc matrix to multiline matrix
	commentator.start("copyBlocMatrixToMultilineMatrix");
		copyBlocMatrixToMultilineMatrix(inMatrix, outMatrix, destruct_in_matrix, NB_THREADS);
	commentator.stop("copyBlocMatrixToMultilineMatrix");

	//if(inMatrix.rowdim () < 1)
	//	return inMatrix.rowdim ();

	//TODO: determine this according to the SIZE of the rows and the number of threads
	__echelonize_global_next_row_to_reduce = NB_THREADS * 2;
	__echelonize_global_last_piv = __echelonize_global_next_row_to_reduce - 1;
	

	//echelonize the few rows before starting parallel
	echelonizeRowUpTo_Sequential(R, outMatrix, 0, __echelonize_global_next_row_to_reduce - 1);

	if(outMatrix.multiline_rowdim () >= __echelonize_global_next_row_to_reduce)
	{

	echelonize_Params_t params;
	params.A = &outMatrix;
	params.R = &R;

	pthread_t threads[NB_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&__echelonize_mutex_lock, NULL);
	pthread_mutex_init(&__echelonize_waiting_list_mutex_lock, NULL);
#else
	pthread_spin_init(&__echelonize_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
	pthread_spin_init(&__echelonize_waiting_list_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
#endif


	for(t=0; t<NB_THREADS; t++){
		//report << "Creating thread " << t << "\n";
		rc = pthread_create(&threads[t], NULL, echelonize__Parallel_in, &params);

		if (rc){		//TODO what happened when only one thread fails
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			throw std::runtime_error ("Cannot create thread");
		}
	}

	/* Free attribute and wait for the other threads */
	for(t=0; t<NB_THREADS; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			throw std::runtime_error ("Cannot create thread");
		}

		//report << "COMPLETED: thread " << t << std::endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&__echelonize_mutex_lock);
	pthread_mutex_destroy(&__echelonize_waiting_list_mutex_lock);
#else
	pthread_spin_destroy(&__echelonize_spinlock_lock);
	pthread_spin_destroy(&__echelonize_waiting_list_spinlock_lock);
#endif

	}

	uint32 rank = 0;

	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;
	uint16 h_a1;

	for(uint32 i=0; i<outMatrix.multiline_rowdim (); ++i)
	{
		if(outMatrix[i].empty ())
			continue;

		head_line1 = Level1Ops::headMultiLineVectorHybrid(outMatrix[i], 0, h_a1, head_line1_idx, outMatrix.coldim ());
		head_line2 = Level1Ops::headMultiLineVectorHybrid(outMatrix[i], 1, h_a1, head_line2_idx, outMatrix.coldim ());

		if(head_line1 != -1)
			rank++;
		if(head_line2 != -1)
			rank++;
	}

	return rank;
}

void* Level3ParallelEchelon::echelonize__Parallel_in(void* p_params)
{
	/*while(true)
	//{
		//LOCK
		
		//GET LAST_GLOBAL_PIV in LAST_CURRENT_GLOBAL_PIV
		
		//IF READY_FOR_WAITING_LIST
			//GET SMALLEST ROW AT WAITING LIST
			//CALL echelonize_one_row(elt.second[==LAST_CURRENT_GLOBAL_PIV], LAST_CURRENT_GLOBAL_PIV)
		//ELSE	
			//GET_NEXT AVAIL ROW
			//CALL echelonize_one_row(0, LAST_CURRENT_GLOBAL_PIV)
			
		//UNLOCK
		
		//IF CURRENT ROW == LAST_GLOBAL_PIV + 1
			//REDUCE IT BY SAME MULTILINE
			//UPDATE LAST_GLOBAL_PIV = CURRENT ROW
			
			//SET READY_FOR_WAITING_LIST = TRUE
		//ELSE
			//SAVE TO WATING LIST: (CURRENT_ROW, LAST_CURRENT_GLOBAL_PIV)
			
		//SAVE BACK CURRENT ROW
			
	}*/
	
	echelonize_Params_t params = *(echelonize_Params_t *)p_params;

	uint32 coldim = params.A->coldim ();
	const uint32 N = params.A->multiline_rowdim ();
	
	uint64 *tmpDenseArray1;
		posix_memalign((void**)&tmpDenseArray1, 16, coldim * sizeof(uint64));
	uint64 *tmpDenseArray2;
		posix_memalign((void**)&tmpDenseArray2, 16, coldim * sizeof(uint64));
	
		
	uint32 current_row_to_reduce;
	uint32 local_last_piv;

	bool ready_for_waiting_list = false;
	bool current_row_fully_reduced = false;
	
	int nb_reduced_consecutively = 0;

	uint32 from_row;	//the pivot to start reducing from
	waiting_row_t waiting_row;
	
//	uint32 ID = pthread_self() % 100;
//	uint32 nb_reduced = 0, nb_waiting_not_fully_reduced = 0, nb_fully_reduced = 0;

	while (true)
	{
			local_last_piv = __echelonize_global_last_piv;
			if(__echelonize_global_last_piv >= N)
			{
				UNLOCK(__echelonize);
				break;
			}

			LOCK(__echelonize);
				if(!ready_for_waiting_list && __echelonize_global_next_row_to_reduce < N)
				{
					current_row_to_reduce = __echelonize_global_next_row_to_reduce;
					++__echelonize_global_next_row_to_reduce;
					from_row = 0;
				}
				else	//all rows are reduced except for those in waiting list, handle them
				{
					ready_for_waiting_list = true;
				}

				if(ready_for_waiting_list)
				{
					if (getSmallestWaitingRow(&waiting_row) == false)	//no waiting rows
					{
						if (__echelonize_global_next_row_to_reduce >= N)	//no more rows to reduce, DONE
						{
							UNLOCK(__echelonize);
							break;
						}
						else if (local_last_piv >= N)
						{
							UNLOCK(__echelonize);
							break;
						}
						else
						{
							ready_for_waiting_list = false;
							UNLOCK(__echelonize);
	//						report << "[" << ID << "] " << "waiting list empty, reducing left rows -> continue" << std::endl;
							continue;
						}
					}

					from_row = waiting_row.last_pivot_reduced_by + 1;
					current_row_to_reduce = waiting_row.row_idx;
				}
			UNLOCK(__echelonize);
		
//		report << "[" << ID << "] "
//				//<< current_row_to_reduce << " - " << local_last_piv << std::endl;
//				<< "reducing " << current_row_to_reduce << " starting from " << from_row
//				<< " to " << local_last_piv
//				<< " [Total " << params.A->multiline_rowdim () -1 << "]" << std::endl << std::endl;
//		if(current_row_to_reduce <= local_last_piv)
//			report << "ERROR" << std::endl;

		Level1Ops::memsetToZero(&tmpDenseArray1, 1, coldim);
        Level1Ops::memsetToZero(&tmpDenseArray2, 1, coldim);

		Level1Ops::copyMultiLineVectorToDenseArray((*params.A)[current_row_to_reduce], tmpDenseArray1, tmpDenseArray2, coldim);
		
//		if(from_row == 0)
//			nb_reduced++;

		echelonize_one_row(*params.R, 
					*params.A, 
					current_row_to_reduce, 
					from_row, 
					local_last_piv,
					tmpDenseArray1, 
					tmpDenseArray2);
		
		if(current_row_to_reduce == local_last_piv + 1)	//have we reduced till the row just below?
		{
			current_row_fully_reduced = true;
			//__echelonize_piv_list[current_row_to_reduce] = current_row_to_reduce;
			ready_for_waiting_list = true;

//			nb_fully_reduced++;
		}
		else
		{
			current_row_fully_reduced = false;
			ready_for_waiting_list = false;		//TODO: add a counter, for example, 
												//once X rows are reduced, consider helping with
												//the waiting rows.
			++nb_reduced_consecutively;
			if(nb_reduced_consecutively >= 5)
			{
				nb_reduced_consecutively = 0;
				ready_for_waiting_list = true;
			}

//			nb_waiting_not_fully_reduced++;
		}

		//SAVE BACK STEP
		saveBackAndReduce(*params.R,
					*params.A,
					current_row_to_reduce,
					tmpDenseArray1,
					tmpDenseArray2,
					current_row_fully_reduced);
		
		if(current_row_fully_reduced)
		{
			LOCK(__echelonize);
				++__echelonize_global_last_piv;
			UNLOCK(__echelonize);
		}
		else
		{
			LOCK(__echelonize);
				pushRowToWaitingList(current_row_to_reduce, local_last_piv);
			UNLOCK(__echelonize);
		}

	}

	free(tmpDenseArray1);
	free(tmpDenseArray2);


//	std::cout << "[" << ID << "] "
//			<< "nb_fully_reduced " << nb_fully_reduced << std::endl
//			<< " nb reduced from 0 " << nb_reduced << std::endl
//			<< " nb_waiting_not_fully_reduced " << nb_waiting_not_fully_reduced << std::endl << std::endl;
	return 0;
}





//Gived a lit of pivots in matrix A, echelonize the row at index idx_row by the list
//of the knwon pivots at the time (from first_pivot to last_pivot)
template <typename Index>
inline void Level3ParallelEchelon::echelonize_one_row(const Modular<uint16>& R,
			SparseMultilineMatrix<uint16, Index>& A,
			uint32 idx_row_to_echelonize, 
			uint32 first_pivot, 
			uint32 last_pivot,
			uint64 *tmpDenseArray1,
			uint64 *tmpDenseArray2)
{
	typedef Modular<uint16> Ring;
	typedef uint16 Element;
	
	uint32 coldim = A.coldim ();
	MultiLineVector<Element, Index> *rowA;
	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;
	
	typename Ring::Element h_a1 = R.one (), h_a2 = R.one ();

	uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
	uint32 tmp=0;
	
	for (uint32 j = first_pivot; j <= last_pivot; ++j)
	{
		rowA = &(A[j]);			//XXX: diff from original algorithm here
		if(rowA->empty ())
		{
			continue;
		}

		head_line1 = Level1Ops::headMultiLineVectorHybrid(*rowA, 0, h_a1, head_line1_idx, coldim);
		head_line2 = Level1Ops::headMultiLineVectorHybrid(*rowA, 1, h_a2, head_line2_idx, coldim);
		
		if(head_line1 != -1)
		{
			//R.copy(v1col1, tmpDenseArray1[head_line1] % R._modulus);
			//R.copy(v2col1, tmpDenseArray2[head_line1] % R._modulus);
			v1col1 = tmpDenseArray1[head_line1] % R._modulus;
			v2col1 = tmpDenseArray2[head_line1] % R._modulus;

			if(v1col1 != 0)
				v1col1 = R._modulus - v1col1;
			if(v2col1 != 0)
				v2col1 = R._modulus - v2col1;
		}
		else
		{
			v1col1 = 0;
			v2col1 = 0;
		}

		if(head_line2 != -1)
		{
			//R.copy(v1col2, tmpDenseArray1[head_line2] % R._modulus);
			//R.copy(v2col2, tmpDenseArray2[head_line2] % R._modulus);
			v1col2 = tmpDenseArray1[head_line2] % R._modulus;
			v2col2 = tmpDenseArray2[head_line2] % R._modulus;

			uint16 val = rowA->at_unchecked(0, head_line2_idx);

			tmp = v1col2 + (uint32)v1col1*val;
			R.init(v1col2, tmp);
			tmp = v2col2 + (uint32)v2col1*val;
			R.init(v2col2, tmp);

			if (v1col2 != 0)
				v1col2 = R._modulus - v1col2;
			if (v2col2 != 0)
				v2col2 = R._modulus - v2col2;
		}
		else
		{
			v1col2 = 0;
			v2col2 = 0;
		}

		if(rowA->is_sparse (coldim))
		{
			Level2Ops::SparseScalMulSub__two_rows__vect_array(
					v1col1,
					v2col1,
					v1col2,
					v2col2,
					*rowA,
					tmpDenseArray1,
					tmpDenseArray2);
		}
		else
		{
			Level2Ops::DenseScalMulSub__two_rows__vect_array__variable_size(
					v1col1,
					v2col1,
					v1col2,
					v2col2,
					*rowA,
					tmpDenseArray1,
					tmpDenseArray2,
					head_line1);
		}
	}
}


template <typename Index>
inline void Level3ParallelEchelon::saveBackAndReduce(const Modular<uint16>& R,
			SparseMultilineMatrix<uint16, Index>& A,
			uint32 row_idx,
			uint64 *tmpDenseArray1,
			uint64 *tmpDenseArray2,
			bool reduce)
{
	uint32 coldim = A.coldim ();

	if(reduce)
	{
		typedef Modular<uint16> Ring;
		long head_line1=-1, head_line2=-1;
		typename Ring::Element h_a2 = R.one ();

		head_line1 = Level1Ops::normalizeDenseArray(R, tmpDenseArray1, coldim);
		head_line2 = Level1Ops::normalizeDenseArray(R, tmpDenseArray2, coldim);

		//assert(h_a1 == 1 || h_a1 == 0);
		//assert(h_a2 == 1 || h_a2 == 0);

		//reduce by same MultilineVector
		if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
		{
			//assert(head_line1 >= 0 && head_line1 < coldim);
			if(tmpDenseArray2[head_line1] % R._modulus != 0)
			{
				register uint32 h = tmpDenseArray2[head_line1] % R._modulus;
				//R.negin(h);
				h = R._modulus - h;

				register uint32 v__;
				for(uint32 x = head_line1; x<coldim; ++x)
				{
					v__ = tmpDenseArray1[x] & 0x000000000000ffff;
					tmpDenseArray2[x] += h * v__;
				}
			}
		}

		//head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
		head_line2 = Level1Ops::headDenseArray(R, tmpDenseArray2, coldim, h_a2);		//only this can change

		if(head_line1 == -1 ||
			((head_line2 != -1) && (head_line1 > head_line2)))			//saves the line with the smallest column entry first
			Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray2, tmpDenseArray1, coldim, A[row_idx]);
		else
			Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray1, tmpDenseArray2, coldim, A[row_idx]);

		Level1Ops::normalizeMultiLineVector(R, A[row_idx]); //should be already normalized at this step if not reduce
	}
	
	if(!reduce)	//do not switch order unless the row is fully reduced
		Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray1, tmpDenseArray2, coldim, A[row_idx]);
	
}



template <typename Index>
uint32 Level3ParallelEchelon::echelonizeRowUpTo_Sequential(const Modular<uint16>& R,
		SparseMultilineMatrix<uint16, Index>& A,
		uint32 from_row,
		uint32 to_row)
{
	typedef Modular<uint16> Ring;
	typedef uint16 Element;

	uint32 coldim = A.coldim ();
	//uint32 npiv = 0;
	uint32 npiv_real = 0;
	uint32 N = A.rowdim()/NB_ROWS_PER_MULTILINE + A.rowdim()%NB_ROWS_PER_MULTILINE;

	//uint32 *piv;
	//posix_memalign((void**)&piv, 16, N * sizeof(uint32));
	//Level1Ops::memsetToZero(&piv, 1, N);

	MultiLineVector<Element, Index> *rowA;

	uint64 *tmpDenseArray1;
	posix_memalign((void**)&tmpDenseArray1, 16, coldim * sizeof(uint64));

	uint64 *tmpDenseArray2;
	posix_memalign((void**)&tmpDenseArray2, 16, coldim * sizeof(uint64));


	uint32 i=0;

	long head_line1=-1, head_line2=-1;
	uint32 head_line1_idx=0, head_line2_idx=0;

	Level1Ops::normalizeMultiLineVector(R, A[from_row]);

	if(A.rowdim() == 0)
		return 0;

	for(i=from_row; i <= min(to_row, N-1) ; ++i)
	{
        Level1Ops::memsetToZero(&tmpDenseArray1, 1, coldim);
        Level1Ops::memsetToZero(&tmpDenseArray2, 1, coldim);

        Level1Ops::copyMultiLineVectorToDenseArray(A[i], tmpDenseArray1, tmpDenseArray2, coldim);

        typename Ring::Element h_a1 = R.one (), h_a2 = R.one ();

		uint16 v1col1=0, v2col1=0, v1col2=0, v2col2=0;
		uint32 tmp=0;

		for(uint32 j=from_row; j<i; ++j)
		{
			rowA = &(A[j]);
			if(rowA->empty ())
			{
				continue;
			}

			head_line1 = Level1Ops::headMultiLineVectorHybrid(*rowA, 0, h_a1, head_line1_idx, coldim);
			head_line2 = Level1Ops::headMultiLineVectorHybrid(*rowA, 1, h_a2, head_line2_idx, coldim);

			if(head_line1 != -1)
			{
				//R.copy(v1col1, tmpDenseArray1[head_line1] % R._modulus);
				//R.copy(v2col1, tmpDenseArray2[head_line1] % R._modulus);
				v1col1 = tmpDenseArray1[head_line1] % R._modulus;
				v2col1 = tmpDenseArray2[head_line1] % R._modulus;

				if(v1col1 != 0)
					v1col1 = R._modulus - v1col1;
				if(v2col1 != 0)
					v2col1 = R._modulus - v2col1;
			}
			else
			{
				v1col1 = 0;
				v2col1 = 0;
			}

			if(head_line2 != -1)
			{
				//R.copy(v1col2, tmpDenseArray1[head_line2] % R._modulus);
				//R.copy(v2col2, tmpDenseArray2[head_line2] % R._modulus);
				v1col2 = tmpDenseArray1[head_line2] % R._modulus;
				v2col2 = tmpDenseArray2[head_line2] % R._modulus;

				uint16 val = rowA->at_unchecked(0, head_line2_idx);

				tmp = v1col2 + (uint32)v1col1*val;
				R.init(v1col2, tmp);
				tmp = v2col2 + (uint32)v2col1*val;
				R.init(v2col2, tmp);

				if (v1col2 != 0)
					v1col2 = R._modulus - v1col2;
				if (v2col2 != 0)
					v2col2 = R._modulus - v2col2;
			}
			else
			{
				v1col2 = 0;
				v2col2 = 0;
			}

			if(rowA->is_sparse (coldim))
			{
				Level2Ops::SparseScalMulSub__two_rows__vect_array(
						v1col1,
						v2col1,
						v1col2,
						v2col2,
						*rowA,
						tmpDenseArray1,
						tmpDenseArray2);
			}
			else
			{
				Level2Ops::DenseScalMulSub__two_rows__vect_array__variable_size(
						v1col1,
						v2col1,
						v1col2,
						v2col2,
						*rowA,
						tmpDenseArray1,
						tmpDenseArray2,
						head_line1);
			}
		}

		head_line1 = Level1Ops::normalizeDenseArray(R, tmpDenseArray1, coldim);
		head_line2 = Level1Ops::normalizeDenseArray(R, tmpDenseArray2, coldim);

        //assert(h_a1 == 1 || h_a1 == 0);
        //assert(h_a2 == 1 || h_a2 == 0);

        //reduce by same MultilineVector
        if(head_line1 >= head_line2 && head_line1 != -1 && head_line2 != -1)
        {
        	//assert(head_line1 >= 0 && head_line1 < coldim);
        	if(tmpDenseArray2[head_line1] % R._modulus != 0)
			{
				register uint32 h = tmpDenseArray2[head_line1] % R._modulus;
				//R.negin(h);
				h = R._modulus - h;

				register uint32 v__;
				for(uint32 x = head_line1; x<coldim; ++x)
				{
					v__ = tmpDenseArray1[x] & 0x000000000000ffff;
					tmpDenseArray2[x] += h * v__;
				}
			}
        }

        //head_line1 = head(R, tmpDenseArray1, coldim, h_a1);
		head_line2 = Level1Ops::headDenseArray(R, tmpDenseArray2, coldim, h_a2);		//only this can change

		if(head_line1 == -1 ||
        		((head_line2 != -1) && (head_line1 > head_line2)))			//saves the line with the smallest column entry first
			Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray2, tmpDenseArray1, coldim, A[i]);
		else
			Level1Ops::copyDenseArraysToMultilineVectorHybrid(R, tmpDenseArray1, tmpDenseArray2, coldim, A[i]);

		Level1Ops::normalizeMultiLineVector(R, A[i]);

		//npiv++;
		/*if(!A[i].empty ())
		{
			piv[npiv] = i;
			npiv++;
		}*/

		if(head_line1 != -1)
			npiv_real++;
		if(head_line2 != -1)
			npiv_real++;

	}

	free(tmpDenseArray1);
	free(tmpDenseArray2);

	return npiv_real;
}


template<typename Index>
void Level3ParallelEchelon::copyBlocMatrixToMultilineMatrix(
	SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
	SparseMultilineMatrix<uint16>& outMatrix,
	bool destruct_in_matrix, int NB_THREADS)
{
	typedef uint16 Element;

	//TODO: can be performed in PARALLEL

	omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
	{
	uint32 curr_row_base = 0;

#pragma omp for schedule(dynamic) nowait
	for(uint32 i=0; i<inMatrix.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
		curr_row_base = i * inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;

		for (uint32 j = 0; j < inMatrix[i].size(); ++j) //for each bloc in the row
		{
			if(inMatrix[i][j].empty ())
				continue;

			//check_equal_or_raise_exception(inMatrix.FirstBlocsColumIndexes[i], 0);
			uint32 bloc_idx = inMatrix.FirstBlocsColumIndexes[i] + (inMatrix.bloc_width () * j);

			uint32 idx;
			Element val1, val2;
			for (uint16 k = 0;
					k < inMatrix[i][j].bloc_height(); ++k) //for each row in the bloc
			{
				if(inMatrix[i][j][k].empty ())
					continue;

				if(inMatrix[i][j][k].is_sparse ())
				{
					for (uint32 p = 0; p < inMatrix[i][j][k].size (); ++p)
					{
						idx = inMatrix[i][j][k].IndexData[p];
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + idx);
						outMatrix[curr_row_base+k].ValuesData.push_back(val1);
						outMatrix[curr_row_base+k].ValuesData.push_back(val2);
					}
				}
				else
				{
					for (uint32 p = 0; p < inMatrix.bloc_width (); ++p)
					{
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						if(val1 != 0 || val2 != 0)
						{
							outMatrix[curr_row_base+k].IndexData.push_back(bloc_idx + p);
							outMatrix[curr_row_base+k].ValuesData.push_back(val1);
							outMatrix[curr_row_base+k].ValuesData.push_back(val2);
						}
					}
				}

				if(destruct_in_matrix)
					inMatrix[i][j][k].free ();
			}

			if(destruct_in_matrix)
				inMatrix[i][j].free ();
		}

	}
	}

	if(destruct_in_matrix)
		inMatrix.free ();
}




#endif /* ECHELON_C_ */
