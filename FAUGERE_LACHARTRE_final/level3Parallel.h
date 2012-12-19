/*
 * level3Parallel.h
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 *  Created on: 4 august 2012
 *      Author: martani (UPMC University Paris 06 / INRIA)
 */

#ifndef LEVEL3PARALLEL_H_
#define LEVEL3PARALLEL_H_


#include "consts-macros.h"

#include "lela/matrix/sparse.h"
#include "level2-ops.h"

using namespace LELA;

class Level3ParallelOps
{

public:
		template<typename Index>
		static void reducePivotsByPivots__Parallel(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			int NB_THREADS);

		template<typename Index>
		static void reduceNonPivotsByPivots__Parallel(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
			bool invert_scalars ,
			int NB_THREADS);
		
		template<typename Index>
		static void reduceNonPivotsByPivots__Parallel_horizontal(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
			bool invert_scalars ,
			int NB_THREADS);

		template<typename Ring>
		static void reduceC__Parallel(const Ring& R,
				const SparseMultilineMatrix<uint16>& A,
				SparseMultilineMatrix<uint16>& C, int NB_THREADS);


		template<typename Index>
		static void reducePivotsByPivots_2_Level_Parallel(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B);

		template<typename Index>
		static void reducePivotsByPivots__Parallel_thread_pool(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B);


		template<typename Element, typename Index>
		static void copyMultilineMatrixToBlocMatrixRTL__Parallel(SparseMultilineMatrix<Element>& A,
				SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B, bool destruct_riginal, int NUM_THREADS);

		typedef struct reduceBlocByRectangularBloc_thread_pool_Params_t {
			const Modular<uint16>* R;
			const SparseMultilineBloc<uint16, IndexType>* bloc_A;
			const SparseMultilineBloc<uint16, IndexType>* bloc_B;
			uint64 **Bloc_acc;
			bool invert_scalars;
			int from;
			int to;
		} reduceBlocByRectangularBloc_thread_pool_Params_t;

private:
		typedef struct ReduceC_Params_t {
			typedef Modular<uint16> Ring;
			const Ring* R;
			const SparseMultilineMatrix<uint16>* A;
			SparseMultilineMatrix<uint16>* C;
		} ReduceC_Params_t;


		struct ReducePivotsByPivots_Params_t {
			typedef Modular<uint16> Ring;
			typedef SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> > BlocMatrix;
			const Ring* R;
			const SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> >* A;
			SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> >* B;
		};

		struct ReduceNonPivotsByPivots_Params_t {
			typedef Modular<uint16> Ring;
			typedef SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> > BlocMatrix;
			const Ring* R;
			const SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> >* C;
			const SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> >* B;
			SparseBlocMatrix<SparseMultilineBloc<uint16, IndexType> >* D;
			bool invert_scalars;
		};

		static void* reduceC__Parallel_in(void* p_params);

		static void* reducePivotsByPivots__Parallel_in(void* p_params);

		static void* reduceNonPivotsByPivots__Parallel_in(void* p_params);
		
		static void* reduceNonPivotsByPivots__Parallel_horizontal_in(void* p_params);

		static void* reduceBlocByRectangularBloc_thread_pool(void *p_params);
};



#include "level3Parallel.C"

#endif /* LEVEL3PARALLEL_H_ */
