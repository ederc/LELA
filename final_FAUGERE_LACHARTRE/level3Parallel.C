/*
 * level3Parallel.C
 *
 *  Created on: 4 août 2012
 *      Author: martani
 */

#include "level3Parallel.h"
#include <pthread.h>

#include "thrpool-0.8/src/TThreadPool.hh"
//#include "ThreadPool.h"

#ifndef LEVEL3PARALLEL_C_
#define LEVEL3PARALLEL_C_




static uint32 __reducePivotsByPivots_next_column_to_reduce;

#ifdef USE_MUTEX
	static pthread_mutex_t __reducePivotsByPivots_mutex_lock;
#else
	static pthread_spinlock_t __reducePivotsByPivots_spinlock_lock;
#endif


template<typename Index>
void Level3ParallelOps::reducePivotsByPivots__Parallel(const Modular<uint16>& R, const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B, int NB_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3ParallelOps::reducePivotsByPivots__Parallel] NB THREADS " << NB_THREADS << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;
	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());


	__reducePivotsByPivots_next_column_to_reduce = 0;
	ReducePivotsByPivots_Params_t params;
	params.A = &A;
	params.B = &B;
	params.R = &R;

	pthread_t threads[NB_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&__reducePivotsByPivots_mutex_lock, NULL);
#else
	pthread_spin_init(&__reducePivotsByPivots_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
#endif

	for(t=0; t<NB_THREADS; t++){
      //report << "Creating thread " << t << "\n";
      rc = pthread_create(&threads[t], NULL, reducePivotsByPivots__Parallel_in, &params);

      if (rc){
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         throw std::runtime_error ("Cannot create thread");
      }
   }

	/* Free attribute and wait for the other threads */
	for(t=0; t<NB_THREADS; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);

		}

		//report << "COMPLETED: thread " << t << " - handled columns: " << (long)status << std::endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&__reducePivotsByPivots_mutex_lock);
#else
	pthread_spin_destroy(&__reducePivotsByPivots_spinlock_lock);
#endif
}


void* Level3ParallelOps::reducePivotsByPivots__Parallel_in(void* p_params)
{
	struct ReducePivotsByPivots_Params_t params = *(struct ReducePivotsByPivots_Params_t *)p_params;

#define CHACHE_LINE_SIZE	64 //bytes

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
	{
		uint32 size = params.B->bloc_width() * sizeof(uint64);
		//size += CHACHE_LINE_SIZE;

		posix_memalign((void**)&dense_bloc[i], 16, size);
	}


	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) params.B->coldim() / params.B->bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) params.A->rowdim() / params.A->bloc_height());


	uint32 local_columns_idx;
	long nb_columns_handled=0;

#ifdef SHOW_PROGRE___SS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	uint32 C_rowdim_multiline = params.C->multiline_rowdim ();
#endif

	while(true)
	{
		LOCK(__reducePivotsByPivots);
			if(__reducePivotsByPivots_next_column_to_reduce < nb_column_blocs_B)
			{
				local_columns_idx = __reducePivotsByPivots_next_column_to_reduce;
				++__reducePivotsByPivots_next_column_to_reduce;
				++nb_columns_handled;
			}
			else
			{
				UNLOCK(__reducePivotsByPivots);
				break;
			}
		UNLOCK(__reducePivotsByPivots);

		//report << "Column B " << i << std::endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 j = 0;
				j < nb_row_blocs_A;
				++j)
		{
			const uint32 first_bloc_idx = params.A->FirstBlocsColumIndexes[j] / params.A->bloc_width();
			const uint32 last_bloc_idx = MIN((*params.A)[j].size () - 1, j);

			//report << "\tRow A " << j << "\tfirst bloc " << first_bloc_idx << " last " << last_bloc_idx << std::endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			Level1Ops::memsetToZero(dense_bloc);

			Level1Ops::copySparseBlocToDenseBlocArray(*params.R, (*params.B)[j][local_columns_idx], dense_bloc);

#ifdef SHOW_PROGRESS
		//report << "                                                                                    \r";
		//report << "\tcolumn\t" << i << "/" << nb_column_blocs_B << "\trow\t" << j << std::ends;
#endif
			//for all the blocs in the current row of A
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in A and B " << k + first_bloc_idx << std::endl;
				Level2Ops::reduceBlocByRectangularBloc(*params.R, (*params.A)[j][k], (*params.B)[k + first_bloc_idx][local_columns_idx], dense_bloc);
			}
			Level2Ops::reduceBlocByTriangularBloc(*params.R, (*params.A)[j][last_bloc_idx], dense_bloc);

			Level1Ops::copyDenseBlocArrayToSparseBloc(*params.R, dense_bloc, (*params.B)[j][local_columns_idx], false);


		}
	}
#ifdef SHOW_PROGRESS
	//report << "\r                                                                                    \n";
#endif

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);

	return (void*) nb_columns_handled;
}






static uint32 __reduceNonPivotsByPivots_next_column_to_reduce;
#ifdef USE_MUTEX
	static pthread_mutex_t __reduceNonPivotsByPivots_mutex_lock;
#else
	static pthread_spinlock_t __reduceNonPivotsByPivots_spinlock_lock;
#endif


template<typename Index>
void Level3ParallelOps::reduceNonPivotsByPivots__Parallel(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
			bool invert_scalars,
			int NB_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3ParallelOps::reduceNonPivotsByPivots__Parallel] NB THREADS " << NB_THREADS << std::endl;


	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(D.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(C.rowdim(), D.rowdim());
	check_equal_or_raise_exception(C.coldim(), B.rowdim());
	check_equal_or_raise_exception(B.coldim(), D.coldim());


	__reduceNonPivotsByPivots_next_column_to_reduce = 0;
	ReduceNonPivotsByPivots_Params_t params;
	params.C = &C;
	params.B = &B;
	params.D = &D;
	params.R = &R;
	params.invert_scalars = invert_scalars;

	pthread_t threads[NB_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&__reduceNonPivotsByPivots_mutex_lock, NULL);
#else
	pthread_spin_init(&__reduceNonPivotsByPivots_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
#endif

	for(t=0; t<NB_THREADS; t++){
      //report << "Creating thread " << t << "\n";
      rc = pthread_create(&threads[t], NULL, reduceNonPivotsByPivots__Parallel_in, &params);

      if (rc){
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         throw std::runtime_error ("Cannot create thread");
      }
   }

	/* Free attribute and wait for the other threads */
	for(t=0; t<NB_THREADS; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);

		}

		//report << "COMPLETED: thread " << t << " - handled columns: " << (long)status << std::endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&__reduceNonPivotsByPivots_mutex_lock);
#else
	pthread_spin_destroy(&__reduceNonPivotsByPivots_spinlock_lock);
#endif
}



void* Level3ParallelOps::reduceNonPivotsByPivots__Parallel_in(void* p_params)
{
	struct ReduceNonPivotsByPivots_Params_t params = *(struct ReduceNonPivotsByPivots_Params_t *)p_params;

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, params.D->bloc_width() * sizeof(uint64));

	const uint32 nb_column_blocs_D = (uint32) std::ceil((double) params.D->coldim() / params.D->bloc_width());
	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)params.C->rowdim() / params.C->bloc_height());


	uint32 local_columns_idx;
	long nb_columns_handled=0;

#ifdef SHOW_PROGRE__SS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif

	while(true)
	{
		LOCK(__reduceNonPivotsByPivots);
			if(__reduceNonPivotsByPivots_next_column_to_reduce < nb_column_blocs_D)
			{
				local_columns_idx = __reduceNonPivotsByPivots_next_column_to_reduce;
				++__reduceNonPivotsByPivots_next_column_to_reduce;
				++nb_columns_handled;
			}
			else
			{
				UNLOCK(__reduceNonPivotsByPivots);
				break;
			}
		UNLOCK(__reduceNonPivotsByPivots);

		//report << "Colum D|B\t" << i << endl;
		//for all rows of blocs in C
		for (uint32 j = 0; j < nb_row_blocs_C; ++j)
		{
			const uint32 first_bloc_idx = params.C->FirstBlocsColumIndexes[j] / params.C->bloc_width();
			const uint32 last_bloc_idx = (*params.C)[j].size ();

			//report << "\tRow C\t" << j << "\tfirst bloc: " << first_bloc_idx << " - last: " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			Level1Ops::memsetToZero(dense_bloc);
			Level1Ops::copySparseBlocToDenseBlocArray(*params.R, (*params.D)[j][local_columns_idx], dense_bloc);

#ifdef SHOW_PROGRE__SS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << local_columns_idx << "/" << nb_column_blocs_D << "\trow\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B\t" << k + first_bloc_idx << endl;
				Level2Ops::reduceBlocByRectangularBloc(*params.R, (*params.C)[j][k], (*params.B)[k + first_bloc_idx][local_columns_idx], dense_bloc, params.invert_scalars);
			}

			Level1Ops::copyDenseBlocArrayToSparseBloc(*params.R, dense_bloc, (*params.D)[j][local_columns_idx]);
		}
	}

#ifdef SHOW_PROGRE__SS
	report << "\r                                                                                    \n";
#endif

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);

	return (void*) nb_columns_handled;
}



static uint32 __reduceNonPivotsByPivots_horizontal_next_row_to_reduce;
#ifdef USE_MUTEX
	static pthread_mutex_t __reduceNonPivotsByPivots_horizontal_mutex_lock;
#else
	static pthread_spinlock_t __reduceNonPivotsByPivots_horizontal_spinlock_lock;
#endif

template<typename Index>
void Level3ParallelOps::reduceNonPivotsByPivots__Parallel_horizontal(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
			bool invert_scalars ,
			int NB_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3ParallelOps::reduceNonPivotsByPivots__Parallel_horizontal] NB THREADS " << NB_THREADS << std::endl;


	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(C.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(D.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(C.rowdim(), D.rowdim());
	check_equal_or_raise_exception(C.coldim(), B.rowdim());
	check_equal_or_raise_exception(B.coldim(), D.coldim());


	__reduceNonPivotsByPivots_horizontal_next_row_to_reduce = 0;
	ReduceNonPivotsByPivots_Params_t params;
	params.C = &C;
	params.B = &B;
	params.D = &D;
	params.R = &R;
	params.invert_scalars = invert_scalars;

	pthread_t threads[NB_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&__reduceNonPivotsByPivots_horizontal_mutex_lock, NULL);
#else
	pthread_spin_init(&__reduceNonPivotsByPivots_horizontal_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
#endif

	for(t=0; t<NB_THREADS; t++){
		//report << "Creating thread " << t << "\n";
		rc = pthread_create(&threads[t], NULL, reduceNonPivotsByPivots__Parallel_horizontal_in, &params);

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

		}

		//report << "COMPLETED: thread " << t << " - handled columns: " << (long)status << std::endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&__reduceNonPivotsByPivots_horizontal_mutex_lock);
#else
	pthread_spin_destroy(&__reduceNonPivotsByPivots_horizontal_spinlock_lock);
#endif
}



void* Level3ParallelOps::reduceNonPivotsByPivots__Parallel_horizontal_in(void* p_params)
{
	struct ReduceNonPivotsByPivots_Params_t params = *(struct ReduceNonPivotsByPivots_Params_t *)p_params;

	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**)&dense_bloc[i], 16, params.D->bloc_width() * sizeof(uint64));

	const uint32 nb_column_blocs_D = (uint32) std::ceil((double) params.D->coldim() / params.D->bloc_width());
	const uint32 nb_row_blocs_C = (uint32)std::ceil((double)params.C->rowdim() / params.C->bloc_height());


	uint32 local_row_idx;
	long nb_rows_handled=0;

	while(true)
	{
		LOCK(__reduceNonPivotsByPivots_horizontal);
			if(__reduceNonPivotsByPivots_horizontal_next_row_to_reduce < nb_row_blocs_C)
			{
				local_row_idx = __reduceNonPivotsByPivots_horizontal_next_row_to_reduce;
				++__reduceNonPivotsByPivots_horizontal_next_row_to_reduce;
				++nb_rows_handled;
			}
			else
			{
				UNLOCK(__reduceNonPivotsByPivots_horizontal);
				break;
			}
		UNLOCK(__reduceNonPivotsByPivots_horizontal);

//		for(local_row_idx = 0; local_row_idx < nb_row_blocs_C; ++local_row_idx)
//		{

		//const uint32 first_bloc_idx = params.C->FirstBlocsColumIndexes[local_row_idx] / params.C->bloc_width();
		const uint32 last_bloc_idx = (*params.C)[local_row_idx].size ();

	
		//report << "Colum D|B\t" << i << endl;
		//for all rows of blocs in C
		for (uint32 j = 0; j < nb_column_blocs_D; ++j)
		{
			//report << "\tRow C\t" << j << "\tfirst bloc: " << first_bloc_idx << " - last: " << last_bloc_idx << endl;
			//1. RazBloc
			//2. copy sparse bloc to Bloc_i_j

			Level1Ops::memsetToZero(dense_bloc);
			Level1Ops::copySparseBlocToDenseBlocArray(*params.R, (*params.D)[local_row_idx][j], dense_bloc);

#ifdef SHOW_PROGRE__SS
		report << "                                                                                    \r";
		report << "\tcolumn\t" << i << "/" << nb_column_blocs_D << "\trow\t" << j << std::ends;
#endif

			//for all the blocs in the current row of C (column of B)
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//report << "\t\tBloc in C and B\t" << k + first_bloc_idx << endl;
				Level2Ops::reduceBlocByRectangularBloc(*params.R, 
								       (*params.C)[local_row_idx][k], 
								       (*params.B)[k][j],
								       dense_bloc, 
								       params.invert_scalars);
			}

			Level1Ops::copyDenseBlocArrayToSparseBloc(*params.R, dense_bloc, (*params.D)[local_row_idx][j]);
		}
	}

#ifdef SHOW_PROGRE__SS
	report << "\r                                                                                    \n";
#endif

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);

	return (void*) nb_rows_handled;
}



















static uint32 __reduceC_next_row_to_reduce;

#ifdef USE_MUTEX
	static pthread_mutex_t __reduceC_mutex_lock;
#else
	static pthread_spinlock_t __reduceC_spinlock_lock;
#endif

template<typename Ring>
void Level3ParallelOps::reduceC__Parallel(const Ring& R, const SparseMultilineMatrix<uint16>& A,
				SparseMultilineMatrix<uint16>& C, int NB_THREADS)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "[Level3ParallelOps::reduceC__Parallel] NB THREADS " << NB_THREADS << std::endl;

	__reduceC_next_row_to_reduce = 0;
	ReduceC_Params_t params;
	params.A = &A;
	params.C = &C;
	params.R = &R;

	pthread_t threads[NB_THREADS];
	void *status;
	int rc;
	int t;

#ifdef USE_MUTEX
	pthread_mutex_init(&__reduceC_mutex_lock, NULL);
#else
	pthread_spin_init(&__reduceC_spinlock_lock, PTHREAD_PROCESS_PRIVATE);
#endif

	for(t=0; t<NB_THREADS; t++){
      //report << "Creating thread " << t << "\n";
      rc = pthread_create(&threads[t], NULL, reduceC__Parallel_in, &params);

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
		}

		//report << "COMPLETED: thread " << t << " - handled rows: " << (long)status << std::endl;
	}

#ifdef USE_MUTEX
	pthread_mutex_destroy(&__reduceC_mutex_lock);
#else
	pthread_spin_destroy(&__reduceC_spinlock_lock);
#endif
}

void* Level3ParallelOps::reduceC__Parallel_in(void* p_params)
{
	ReduceC_Params_t params = *(ReduceC_Params_t *)p_params;

	uint32 C_coldim = params.C->coldim ();

	uint64 *tmpDenseArray1C;
		posix_memalign((void**)&tmpDenseArray1C, 16, C_coldim * sizeof(uint64));

	uint64 *tmpDenseArray2C;
		posix_memalign((void**)&tmpDenseArray2C, 16, C_coldim * sizeof(uint64));

	uint32 local_row_idx;
	long nb_rows_handled=0;

#ifdef SHOW_PROGRE___SS
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	uint32 C_rowdim_multiline = params.C->multiline_rowdim ();
#endif

	while(true)
	{

		LOCK(__reduceC);
			if(__reduceC_next_row_to_reduce < params.C->multiline_rowdim ())
			{
				local_row_idx = __reduceC_next_row_to_reduce;
				++__reduceC_next_row_to_reduce;
				++nb_rows_handled;
			}
			else
			{
				UNLOCK(__reduceC);
				break;
			}
		UNLOCK(__reduceC);
	
	
#ifdef SHOW_PROGRE___SS
		report << "                                                                    \r";
		report << "\trow\t" << next_row_to_reduce << " / " << C_rowdim_multiline << std::ends;
#endif

		Level1Ops::memsetToZero(&tmpDenseArray1C, 1, C_coldim);
		Level1Ops::memsetToZero(&tmpDenseArray2C, 1, C_coldim);

		Level1Ops::copyMultiLineVectorToDenseArray((*params.C)[local_row_idx], tmpDenseArray1C, tmpDenseArray2C, C_coldim);

		uint32 start_idx;
		if((*params.C)[local_row_idx].is_sparse (C_coldim))
			start_idx = (*params.C)[local_row_idx].IndexData[0];
		else
			start_idx = 0;		//should not happen!

		ReduceC_Params_t::Ring::Element Cv1_col1=0, Cv2_col1=0;
		uint32 Cp1=0;
		ReduceC_Params_t::Ring::Element Cv1_col2=0, Cv2_col2=0;
		uint32 Cp2=0;

		uint32 tmp=0;
		uint32 row_in_A;

		for(uint32 j=start_idx; j<C_coldim; ++j)
		{
			Cp1 = j;
			Cv1_col1 = tmpDenseArray1C[j] % params.R->_modulus;
			Cv2_col1 = tmpDenseArray2C[j] % params.R->_modulus;

			if(Cv1_col1 == 0 && Cv2_col1 == 0)
				continue;

			if (Cv1_col1 != 0)
				Cv1_col1 = params.R->_modulus - Cv1_col1;
			if (Cv2_col1 != 0)
				Cv2_col1 = params.R->_modulus - Cv2_col1;

			row_in_A = (params.A->rowdim() - 1 - Cp1) / NB_ROWS_PER_MULTILINE;

			if ((params.A->rowdim() - 1 - Cp1) % 2 == 1)
			{
				Cp2 = j+1;
				Cv1_col2 = tmpDenseArray1C[Cp2] % params.R->_modulus;
				Cv2_col2 = tmpDenseArray2C[Cp2] % params.R->_modulus;

				uint16 v__ = (*params.A)[row_in_A].at_unchecked(1, 1);

				tmp = Cv1_col2 + (uint32)Cv1_col1 * v__;
				params.R->init(Cv1_col2, tmp);
				tmp = Cv2_col2 + (uint32)Cv2_col1 * v__;
				params.R->init(Cv2_col2, tmp);

				if (Cv1_col2 != 0)
					Cv1_col2 = params.R->_modulus - Cv1_col2;
				if (Cv2_col2 != 0)
					Cv2_col2 = params.R->_modulus - Cv2_col2;

				Level2Ops::SparseScalMulSub__two_rows__vect_array(
						Cv1_col2,
						Cv2_col2,
						Cv1_col1,
						Cv2_col1,
						(*params.A)[row_in_A],
						tmpDenseArray1C,
						tmpDenseArray2C);

				tmpDenseArray1C[Cp1] = Cv1_col1;
				tmpDenseArray2C[Cp1] = Cv2_col1;

				tmpDenseArray1C[Cp2] = Cv1_col2;
				tmpDenseArray2C[Cp2] = Cv2_col2;

				++j;
			}
			else
			{
					Level2Ops::SparseScalMulSub__one_row__vect_array(
							Cv1_col1,
							Cv2_col1,
							(*params.A)[row_in_A],
							(params.A->rowdim() - 1 - Cp1) % NB_ROWS_PER_MULTILINE,
							tmpDenseArray1C,
							tmpDenseArray2C);

					tmpDenseArray1C[Cp1] = Cv1_col1;
					tmpDenseArray2C[Cp1] = Cv2_col1;
			}
		}

		Level1Ops::copyDenseArraysToMultilineVector(*params.R, tmpDenseArray1C, tmpDenseArray2C, C_coldim, (*params.C)[local_row_idx]);

	}

#ifdef SHOW_PROGRE___SS
	report << "\r                                                                    \n";
#endif

		free(tmpDenseArray1C);
		free(tmpDenseArray2C);

		return (void*) nb_rows_handled;
}


static void report_num_threads(int level)
{
    #pragma omp single
    {
         printf("\nLevel %d: number of threads in the team - %d\n\n", level, omp_get_num_threads());
    }
}


/************************************************************************************************/
template<typename Index>
void Level3ParallelOps::reducePivotsByPivots_2_Level_Parallel(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "reducePivotsByPivots_2_Level_Parallel" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());

	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) B.coldim() / B.bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) A.rowdim() / A.bloc_height());


#ifdef L1	//LEVEL 1 PARALLELISM
omp_set_nested(1);
omp_set_dynamic(0);
#pragma omp parallel num_threads(NUM_THREADS_OMP_MASTER)
{
#endif
	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT] __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**) &dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	report_num_threads(1);

	//for all columns of blocs of B
#ifdef L1	//LEVEL 1 PARALLELISM
	#pragma omp for schedule(dynamic) nowait
#endif
	for (uint32 ii = 0; ii < nb_column_blocs_B; ++ii)
	{
		//report << "Column B " << i << std::endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 jj = 0; jj < nb_row_blocs_A; ++jj)
		{
			const uint32 first_bloc_idx = A.FirstBlocsColumIndexes[jj] / A.bloc_width();
			const uint32 last_bloc_idx = MIN(A[jj].size () - 1, jj);

			Level1Ops::memsetToZero(dense_bloc);
			Level1Ops::copySparseBlocToDenseBlocArray(R, B[jj][ii], dense_bloc);

			//for all the blocs in the current row of A

			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//Level2Ops::reduceBlocByRectangularBloc(R, A[j][k], B[k + first_bloc_idx][i], dense_bloc);

				if (A[jj][k].empty() || B[k + first_bloc_idx][ii].empty())
					continue;

#ifdef L2 //LEVEL 2!
				#pragma omp parallel num_threads(NUM_THREADS_OMP_SLAVES_PER_MASTER)
				{
					//report_num_threads(2);

				#pragma omp for schedule(dynamic) nowait
#endif
				for (int i = 0; i < DEFAULT_BLOC_HEIGHT / 2; ++i)
				{
					uint8 is_sparse = 0;

					if (A[jj][k][i].is_sparse())
						is_sparse = 1;
					else
						is_sparse = 0;

					const Index N = is_sparse == 1 ? A[jj][k][i].size() : DEFAULT_BLOC_WIDTH;

					for (uint32 j = 0; j < N; ++j)
					{
						const Index Ap1 = (is_sparse == 1 ? A[jj][k][i].IndexData[j] : j);

						//R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
						register uint32 Av1_col1 = A[jj][k][i].at_unchecked(0, j);
						//R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));
						register uint32 Av2_col1 = A[jj][k][i].at_unchecked(1, j);

						if (Av1_col1 != 0)
							Av1_col1 = (uint32) R._modulus - Av1_col1;
						if (Av2_col1 != 0)
							Av2_col1 = (uint32) R._modulus - Av2_col1;

						if (((Ap1 % 2) == 0) && (j < (uint32) (N - 1)))
						{
							const Index Ap2 = (is_sparse == 1 ? A[jj][k][i].IndexData[j + 1] : j + 1);
							if (Ap2 == Ap1 + 1) //axpy 2 ROWS
							{
								// 					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
								// 					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
								register uint32 Av1_col2 = A[jj][k][i].at_unchecked(0, j + 1);
								register uint32 Av2_col2 = A[jj][k][i].at_unchecked(1, j + 1);

								if (Av1_col2 != 0)
									Av1_col2 = (uint32) R._modulus - Av1_col2;
								if (Av2_col2 != 0)
									Av2_col2 = (uint32) R._modulus - Av2_col2;

								++j;

								if (B[k + first_bloc_idx][ii][Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
								{
									Level2Ops::SparseScalMulSub__two_rows__vect_array(
											Av1_col1, Av2_col1, Av1_col2,
											Av2_col2,
											B[k + first_bloc_idx][ii][Ap1
													/ NB_ROWS_PER_MULTILINE],
											dense_bloc[i * 2],
											dense_bloc[i * 2 + 1]);
								}
								else
								{
									Level2Ops::DenseScalMulSub__two_rows__vect_array(
											Av1_col1, Av2_col1, Av1_col2,
											Av2_col2,
											B[k + first_bloc_idx][ii][Ap1
													/ NB_ROWS_PER_MULTILINE],
											dense_bloc[i * 2],
											dense_bloc[i * 2 + 1]);
								}
							}
							else //axpy ONE ROW
							{
								if (B[k + first_bloc_idx][ii][Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
								{
									Level2Ops::SparseScalMulSub__one_row__vect_array(
											Av1_col1, Av2_col1,
											B[k + first_bloc_idx][ii][Ap1
													/ NB_ROWS_PER_MULTILINE],
											Ap1 % NB_ROWS_PER_MULTILINE,
											dense_bloc[i * 2],
											dense_bloc[i * 2 + 1]);
								}
								else
								{
									Level2Ops::DenseScalMulSub__one_row__vect_array(
											Av1_col1, Av2_col1,
											B[k + first_bloc_idx][ii][Ap1
													/ NB_ROWS_PER_MULTILINE],
											Ap1 % NB_ROWS_PER_MULTILINE,
											dense_bloc[i * 2],
											dense_bloc[i * 2 + 1]);
								}
							}
						}
						else //axpy ONE ROW
						{
							if (B[k + first_bloc_idx][ii][Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
							{
								Level2Ops::SparseScalMulSub__one_row__vect_array(
										Av1_col1, Av2_col1,
										B[k + first_bloc_idx][ii][Ap1
												/ NB_ROWS_PER_MULTILINE],
										Ap1 % NB_ROWS_PER_MULTILINE,
										dense_bloc[i * 2],
										dense_bloc[i * 2 + 1]);
							}
							else
							{
								Level2Ops::DenseScalMulSub__one_row__vect_array(
										Av1_col1, Av2_col1,
										B[k + first_bloc_idx][ii][Ap1
												/ NB_ROWS_PER_MULTILINE],
										Ap1 % NB_ROWS_PER_MULTILINE,
										dense_bloc[i * 2],
										dense_bloc[i * 2 + 1]);
							}

						}
					}
				}

			}
#ifdef L2 //LEVEL 2
				} //pragma openmp parallel
#endif


			Level2Ops::reduceBlocByTriangularBloc(R, A[jj][last_bloc_idx], dense_bloc);
			Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[jj][ii], false);

		}
	}

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);

#ifdef L1 //LEVEL 1
}
#endif

}











class reduceBlocByRectangularBlocJob: public ThreadPool::TPool::TJob
{
public:
	reduceBlocByRectangularBlocJob ( int p ) : ThreadPool::TPool::TJob(p) {}

	void run(void * arg)
	{
		Level3ParallelOps::reduceBlocByRectangularBloc_thread_pool_Params_t params =
				*(Level3ParallelOps::reduceBlocByRectangularBloc_thread_pool_Params_t *) arg;

		typedef Modular<uint16> Ring;

		if (params.bloc_A->empty() || params.bloc_B->empty())
			return;

		uint64 **Bloc_acc = (uint64**) params.Bloc_acc;

		for (int i = params.from; i < params.to; ++i)
		{
			uint8 is_sparse = 0;

			const MultiLineVector<uint16, IndexType> *rowA =
					&(*params.bloc_A)[i];

			if (rowA->is_sparse())
				is_sparse = 1;
			else
				is_sparse = 0;

			const IndexType N =
					is_sparse == 1 ? rowA->size() : DEFAULT_BLOC_WIDTH;

			for (uint32 j = 0; j < N; ++j)
			{
				const IndexType Ap1 = (is_sparse == 1 ? rowA->IndexData[j] : j);

				//R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
				register uint32 Av1_col1 = rowA->at_unchecked(0, j);
				//R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));
				register uint32 Av2_col1 = rowA->at_unchecked(1, j);

				if (params.invert_scalars)
				{
					if (Av1_col1 != 0)
						Av1_col1 = (uint32) params.R->_modulus - Av1_col1;
					if (Av2_col1 != 0)
						Av2_col1 = (uint32) params.R->_modulus - Av2_col1;
				}

				if (((Ap1 % 2) == 0) && (j < (uint32) (N - 1)))
				{
					const IndexType Ap2 = (
							is_sparse == 1 ? rowA->IndexData[j + 1] : j + 1);
					if (Ap2 == Ap1 + 1) //axpy 2 ROWS
					{
// 					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
// 					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
						register uint32 Av1_col2 = rowA->at_unchecked(0, j + 1);
						register uint32 Av2_col2 = rowA->at_unchecked(1, j + 1);

						if (params.invert_scalars)
						{
							if (Av1_col2 != 0)
								Av1_col2 = (uint32) params.R->_modulus
										- Av1_col2;
							if (Av2_col2 != 0)
								Av2_col2 = (uint32) params.R->_modulus
										- Av2_col2;
						}

						++j;

						if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
						{
							Level2Ops::SparseScalMulSub__two_rows__vect_array(
									Av1_col1, Av2_col1, Av1_col2, Av2_col2,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
						else
						{
							Level2Ops::DenseScalMulSub__two_rows__vect_array(
									Av1_col1, Av2_col1, Av1_col2, Av2_col2,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
					}
					else //axpy ONE ROW
					{
						if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
						{
							Level2Ops::SparseScalMulSub__one_row__vect_array(
									Av1_col1, Av2_col1,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
						else
						{
							Level2Ops::DenseScalMulSub__one_row__vect_array(
									Av1_col1, Av2_col1,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
					}
				}
				else //axpy ONE ROW
				{
					if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
					{
						Level2Ops::SparseScalMulSub__one_row__vect_array(
								Av1_col1, Av2_col1,
								(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE, Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
					else
					{
						Level2Ops::DenseScalMulSub__one_row__vect_array(
								Av1_col1, Av2_col1,
								(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE, Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}

				}
				//__builtin_prefetch (&(bloc_A[i].IndexData[j+1]), __PREFETCH_READ, __PREFETCH_LOCALITY_LOW);
				//__builtin_prefetch (&(bloc_A[i].ValuesData[(j+1)*NB_ROWS_PER_MULTILINE]), __PREFETCH_READ, __PREFETCH_LOCALITY_MODERATE);
			}
		}

	}
};



void* reduceBlocByRectangularBlocJob_func(void * arg)
{
		Level3ParallelOps::reduceBlocByRectangularBloc_thread_pool_Params_t params =
				*(Level3ParallelOps::reduceBlocByRectangularBloc_thread_pool_Params_t *) arg;

		typedef Modular<uint16> Ring;

		if (params.bloc_A->empty() || params.bloc_B->empty())
			return 0;

		uint64 **Bloc_acc = (uint64**) params.Bloc_acc;

		for (int i = params.from; i < params.to; ++i)
		{
			uint8 is_sparse = 0;

			const MultiLineVector<uint16, IndexType> *rowA =
					&(*params.bloc_A)[i];

			if (rowA->is_sparse())
				is_sparse = 1;
			else
				is_sparse = 0;

			const IndexType N =
					is_sparse == 1 ? rowA->size() : DEFAULT_BLOC_WIDTH;

			for (uint32 j = 0; j < N; ++j)
			{
				const IndexType Ap1 = (is_sparse == 1 ? rowA->IndexData[j] : j);

				//R.copy(Av1_col1, bloc_A[i].at_unchecked(0, j));
				register uint32 Av1_col1 = rowA->at_unchecked(0, j);
				//R.copy(Av2_col1, bloc_A[i].at_unchecked(1, j));
				register uint32 Av2_col1 = rowA->at_unchecked(1, j);

				if (params.invert_scalars)
				{
					if (Av1_col1 != 0)
						Av1_col1 = (uint32) params.R->_modulus - Av1_col1;
					if (Av2_col1 != 0)
						Av2_col1 = (uint32) params.R->_modulus - Av2_col1;
				}

				if (((Ap1 % 2) == 0) && (j < (uint32) (N - 1)))
				{
					const IndexType Ap2 = (
							is_sparse == 1 ? rowA->IndexData[j + 1] : j + 1);
					if (Ap2 == Ap1 + 1) //axpy 2 ROWS
					{
// 					R.copy(Av1_col2, bloc_A[i].at_unchecked(0, j+1));
// 					R.copy(Av2_col2, bloc_A[i].at_unchecked(1, j+1));
						register uint32 Av1_col2 = rowA->at_unchecked(0, j + 1);
						register uint32 Av2_col2 = rowA->at_unchecked(1, j + 1);

						if (params.invert_scalars)
						{
							if (Av1_col2 != 0)
								Av1_col2 = (uint32) params.R->_modulus
										- Av1_col2;
							if (Av2_col2 != 0)
								Av2_col2 = (uint32) params.R->_modulus
										- Av2_col2;
						}

						++j;

						if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
						{
							Level2Ops::SparseScalMulSub__two_rows__vect_array(
									Av1_col1, Av2_col1, Av1_col2, Av2_col2,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
						else
						{
							Level2Ops::DenseScalMulSub__two_rows__vect_array(
									Av1_col1, Av2_col1, Av1_col2, Av2_col2,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
					}
					else //axpy ONE ROW
					{
						if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
						{
							Level2Ops::SparseScalMulSub__one_row__vect_array(
									Av1_col1, Av2_col1,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
						else
						{
							Level2Ops::DenseScalMulSub__one_row__vect_array(
									Av1_col1, Av2_col1,
									(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
									Ap1 % NB_ROWS_PER_MULTILINE,
									Bloc_acc[i * 2], Bloc_acc[i * 2 + 1]);
						}
					}
				}
				else //axpy ONE ROW
				{
					if ((*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE].is_sparse())
					{
						Level2Ops::SparseScalMulSub__one_row__vect_array(
								Av1_col1, Av2_col1,
								(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE, Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}
					else
					{
						Level2Ops::DenseScalMulSub__one_row__vect_array(
								Av1_col1, Av2_col1,
								(*params.bloc_B)[Ap1 / NB_ROWS_PER_MULTILINE],
								Ap1 % NB_ROWS_PER_MULTILINE, Bloc_acc[i * 2],
								Bloc_acc[i * 2 + 1]);
					}

				}
				//__builtin_prefetch (&(bloc_A[i].IndexData[j+1]), __PREFETCH_READ, __PREFETCH_LOCALITY_LOW);
				//__builtin_prefetch (&(bloc_A[i].ValuesData[(j+1)*NB_ROWS_PER_MULTILINE]), __PREFETCH_READ, __PREFETCH_LOCALITY_MODERATE);
			}
		}

		return 0;
	}







template<typename Index>
void Level3ParallelOps::reducePivotsByPivots__Parallel_thread_pool(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	report << "reducePivotsByPivots__Parallel_thread_pool" << std::endl;

	typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

	check_equal_or_raise_exception(A.blocArrangement, Matrix::ArrangementDownTop_RightLeft);
	check_equal_or_raise_exception(B.blocArrangement, Matrix::ArrangementDownTop_LeftRight);
	check_equal_or_raise_exception(A.rowdim(), B.rowdim());
	check_equal_or_raise_exception(A.rowdim(), A.coldim());

	check_equal_or_raise_exception(A.bloc_height(), B.bloc_height());
	check_equal_or_raise_exception(A.bloc_height(), A.bloc_width());

// 	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT]  __attribute__((aligned(0x1000)));
// 	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
// 		posix_memalign((void**)&dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	TIMER_DECLARE_(memsetBlocToZero);
	TIMER_DECLARE_(copySparseBlocToDenseBlocArray);
	TIMER_DECLARE_(reduceBlocByRectangularBloc);
	TIMER_DECLARE_(reduceBlocByTriangularBloc);
	TIMER_DECLARE_(copyDenseBlocArrayToSparseBloc);

	const uint32 nb_column_blocs_B = (uint32) std::ceil((double) B.coldim() / B.bloc_width());
	const uint32 nb_row_blocs_A = (uint32) std::ceil((double) A.rowdim() / A.bloc_height());

omp_set_dynamic(0);
#pragma omp parallel num_threads(NUM_THREADS_OMP_MASTER)
{
	uint64 *dense_bloc[DEFAULT_BLOC_HEIGHT] __attribute__((aligned(0x1000)));
	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		posix_memalign((void**) &dense_bloc[i], 16, B.bloc_width() * sizeof(uint64));

	report_num_threads(1);

	ThreadPool::TPool *pool = new ThreadPool::TPool(NUM_THREADS_OMP_SLAVES_PER_MASTER);
	reduceBlocByRectangularBlocJob *runner1 = new reduceBlocByRectangularBlocJob (0);
	reduceBlocByRectangularBlocJob *runner2 = new reduceBlocByRectangularBlocJob (1);
//	ThreadPool pool;
//	pool.initThreadPool(NUM_THREADS_OMP_SLAVES_PER_MASTER);


	reduceBlocByRectangularBloc_thread_pool_Params_t params_struct[2];
	params_struct[0].R = &R;
	params_struct[0].Bloc_acc = dense_bloc;
	params_struct[0].invert_scalars = true;
	params_struct[0].from = 0;
	params_struct[0].to = DEFAULT_BLOC_WIDTH / 4;

	params_struct[1].R = &R;
	params_struct[1].Bloc_acc = dense_bloc;
	params_struct[1].invert_scalars = true;
	params_struct[1].from = DEFAULT_BLOC_WIDTH / 4;
	params_struct[1].to = DEFAULT_BLOC_WIDTH / 2;


	#pragma omp for schedule(dynamic)  nowait
	for (uint32 i = 0; i < nb_column_blocs_B; ++i)
	{
		//report << "Column B " << i << std::endl;
		//for all rows of blocs in A (starting from the end)
		for (uint32 j = 0;
				j < nb_row_blocs_A;
				++j)
		{
			const uint32 first_bloc_idx = A.FirstBlocsColumIndexes[j] / A.bloc_width();
			const uint32 last_bloc_idx = min(A[j].size () - 1, j);

			Level1Ops::memsetToZero(dense_bloc);
			Level1Ops::copySparseBlocToDenseBlocArray(R, B[j][i], dense_bloc);

			//for all the blocs in the current row of A
			for (uint32 k = 0; k < last_bloc_idx; ++k)
			{
				//Level2Ops::reduceBlocByRectangularBloc(R, A[j][k], B[k + first_bloc_idx][i], dense_bloc);
				params_struct[0].bloc_A = &(A[j][k]);
				params_struct[0].bloc_B = &(B[k + first_bloc_idx][i]);

				params_struct[1].bloc_A = &(A[j][k]);
				params_struct[1].bloc_B = &(B[k + first_bloc_idx][i]);

				pool->run(runner1, &params_struct[0], false);
				pool->run(runner2, &params_struct[1], false);

				pool->sync_all();
//				pool.queueTask(reduceBlocByRectangularBlocJob_func, &params_struct[0]);
//				pool.queueTask(reduceBlocByRectangularBlocJob_func, &params_struct[1]);
//
//				pool.waitAllTasks();

			}

			Level2Ops::reduceBlocByTriangularBloc(R, A[j][last_bloc_idx], dense_bloc);
			Level1Ops::copyDenseBlocArrayToSparseBloc(R, dense_bloc, B[j][i], false);

		}
	}
#ifdef SHOW_PROGRESS
	report << "\r                                                                                    \n";
#endif

	for (uint32 i = 0; i < DEFAULT_BLOC_HEIGHT; ++i)
		free(dense_bloc[i]);

}
	std::cout << "HERE" << std::endl;
}



template<typename Element, typename Index>
void Level3ParallelOps::copyMultilineMatrixToBlocMatrixRTL__Parallel(SparseMultilineMatrix<Element>& A,
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B, bool destruct_riginal, int NUM_THREADS)
{
	typedef SparseBlocMatrix<SparseMultilineBloc<Element, Index> > Matrix;
	B =  Matrix (A.rowdim(), A.coldim(), Matrix::ArrangementDownTop_RightLeft);


	const uint32 nb_blocs_per_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());

	//INIT BLOCS
	for(uint32 i=0; i<=A.multiline_rowdim()/(B.bloc_height() / NB_ROWS_PER_MULTILINE); ++i)
	{
		B[i].resize (nb_blocs_per_dim);

		for(uint32 k=0; k<nb_blocs_per_dim; ++k)
			B[i][k].init(B.bloc_height (), B.bloc_width ());

		B.FirstBlocsColumIndexes[i] = 0;
	}

#pragma omp parallel num_threads(NUM_THREADS)
	{
		uint32 row_bloc_idx = 0;
		uint32 idx;
		Element e1, e2;
		Index elt_idx_in_line;
		uint32 bloc_idx_in_row, line_idx_in_bloc;

#pragma omp for schedule(dynamic) nowait ordered
		for(uint32 i = 0; i<A.multiline_rowdim(); ++i)
		{
			row_bloc_idx = i / (B.bloc_height() / NB_ROWS_PER_MULTILINE);
			line_idx_in_bloc = i % (B.bloc_height() / NB_ROWS_PER_MULTILINE);
			//cout << "i " << i << " row bloc " << row_bloc_idx << endl;

//			if(line_idx_in_bloc == 0)
//			{
//				B[row_bloc_idx].resize (nb_blocs_per_dim);
//
//				for(uint32 k=0; k<nb_blocs_per_dim; ++k)
//					B[row_bloc_idx][k].init(B.bloc_height (), B.bloc_width ());
//
//				B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
//			}

			for(int j=A[i].size ()-1; j>=0; --j)	//sparse
			{
				idx = A[i].IndexData[j];
				e1 = A[i].at_unchecked(0, j);
				e2 = A[i].at_unchecked(1, j);

				bloc_idx_in_row = (A.coldim () - 1 - idx) / B.bloc_width();
				elt_idx_in_line = (A.coldim () - 1 - idx) % B.bloc_width();

				B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back(elt_idx_in_line);
				B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back(e1);
				B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back(e2);
			}

			if(destruct_riginal)
				A[i].free ();

		}


	}

	if(destruct_riginal)
		A.free ();

}















#endif /* LEVEL3PARALLEL_C_ */
