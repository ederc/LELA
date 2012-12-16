/*
 * level3-ops.h
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */

#ifndef LEVEL3_OPS_H_
#define LEVEL3_OPS_H_

#include "lela/matrix/sparse.h"
#include "level2-ops.h"

using namespace LELA;

class Level3Ops
{
public:

	template<typename Index>
	static void reducePivotsByPivots(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B);


	template<typename Index>
	static void reduceNonPivotsByPivots(const Modular<uint16>& R,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
		const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
		SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
		bool invert_scalars = true);

	template<typename Index>
	static void reduceC(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C);

	template<typename Index>
	static void reduceC(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseMultilineMatrix<uint16>& C);

	template<typename Ring>
	static void reduceC(const Ring& R,
			const SparseMultilineMatrix<uint16>& A,
			SparseMultilineMatrix<uint16>& C);

	template<typename Index>
	static uint32 echelonize(const Modular<uint16>& R,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& inMatrix,
			SparseMultilineMatrix<uint16>& outMatrix, bool destruct_in_matrix = true, bool use_hybrid_method = true);

	template<typename Element, typename Index>
	static void copyMultilineMatrixToBlocMatrixRTL(SparseMultilineMatrix<Element>& A,
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B, bool destruct_riginal = true);


	/***************************************************************************************************/
	template<typename Index>
	static void reducePivotsByPivots_horizontal(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B);
	
	template<typename Index>
	static void reduceNonPivotsByPivots_horizontal(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D,
			bool invert_scalars = true);
	
private:
	Level3Ops() {}
	Level3Ops(const Level3Ops& other) {}

	template <typename Index>
	static uint32 echelonize_in_sparse(const Modular<uint16>& R,
			SparseMultilineMatrix<uint16, Index>& A);

	template <typename Index>
	static uint32 echelonize_in_hybrid(const Modular<uint16>& R,
			SparseMultilineMatrix<uint16, Index>& A);

};
















#include "level3-ops.C"

#endif /* LEVEL3_OPS_H_ */
