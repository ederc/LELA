/*
 * echelon.h
 *
 *  Created on: 4 août 2012
 *      Author: martani
 */

#ifndef ECHELON_H_
#define ECHELON_H_



#include "lela/matrix/sparse.h"
#include "level2-ops.h"

using namespace LELA;

class Level3ParallelEchelon
{

public:
	template<typename Index>
	static uint32 echelonize__Parallel(const Modular<uint16>& R,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
			SparseMultilineMatrix<uint16>& outMatrix, bool destruct_in_matrix,
			int NB_THREADS);

	typedef struct waiting_row_t {
			uint32 row_idx;
			uint32 last_pivot_reduced_by;
	} waiting_row_t;

private:
		typedef struct echelonize_Params_t {
			typedef Modular<uint16> Ring;
			const Ring* R;
			SparseMultilineMatrix<uint16>* A;
		} echelonize_Params_t;
		
		static void* echelonize__Parallel_in(void* p_params);

		template<typename Index>
		static void copyBlocMatrixToMultilineMatrix(SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
			SparseMultilineMatrix<uint16>& outMatrix, bool destruct_in_matrix, int NB_THREADS);

		static bool getSmallestWaitingRow(waiting_row_t* elt);

		static void pushRowToWaitingList(uint32 row_idx, uint32 last_pivot_reduced_by);

		template <typename Index>
		static inline void echelonize_one_row(const Modular<uint16>& R,
					SparseMultilineMatrix<uint16, Index>& A,
					uint32 idx_row_to_echelonize,
					uint32 first_pivot,
					uint32 last_pivot,
					uint64 *tmpDenseArray1,
					uint64 *tmpDenseArray2);

		template <typename Index>
		static inline void saveBackAndReduce(const Modular<uint16>& R,
					SparseMultilineMatrix<uint16, Index>& A,
					uint32 row_idx,
					uint64 *tmpDenseArray1,
					uint64 *tmpDenseArray2,
					bool reduce);

		template <typename Index>
		static uint32 echelonizeRowUpTo_Sequential(const Modular<uint16>& R,
				SparseMultilineMatrix<uint16, Index>& A,
				uint32 from_row,
				uint32 to_row);
};



#include "level3Parallel_echelon.C"

#endif /* ECHELON_H_ */
