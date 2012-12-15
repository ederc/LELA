/*
 * matrix-ops.h
 *
 *  Created on: 6 juil. 2012
 *      Author: martani
 */

#ifndef MATRIX_OPS_H_
#define MATRIX_OPS_H_

#include "lela/matrix/sparse.h"

using namespace LELA;

class MatrixOps {

public:

	/**
	 * Trsm solve: B = A^-1 B
	 */
	/*template<typename Element>
	static void reducePivotsByPivots(const Modular<Element>& R,
			const SparseBlocMatrix<Element>& A, SparseBlocMatrix<Element>& B);*/

	/**
	 * Trsm solve: B = A^-1 B
	 */
	template<typename Index>
	static void reducePivotsByPivots(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B);

	/*template<typename Element>
	static void reduceNonPivotsByPivots(const Modular<Element>& R,
		const SparseBlocMatrix<Element>& C, const SparseBlocMatrix<Element>& B,
		SparseBlocMatrix<Element>& D);*/

	/*template <typename Index>
	static void reduceNonPivotsByPivots(const Modular<uint16>& R,
			const SparseBlocMatrix<uint16, Index>& C, const SparseBlocMatrix<uint16, Index>& B,
			SparseBlocMatrix<uint16, Index>& D);*/

	template<typename Index>
	static void reduceNonPivotsByPivots(const Modular<uint16>& R,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& C,
			const SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<uint16, Index> >& D);
};

#include "matrix-ops.C"

#endif /* MATRIX_OPS_H_ */
