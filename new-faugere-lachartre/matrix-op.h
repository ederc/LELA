/*
 * matrix-op.h
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 */

#ifndef MATRIX_OP_H_
#define MATRIX_OP_H_


#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace LELA;

class MatrixOp {
public:

	template <typename Ring>
	static inline void razArray(const Ring& R, typename Ring::Element arr[], uint32 arrSize);

	template <typename Ring, typename Iterator>
	static inline void copySparseVectorToDenseArray(Ring R, const Iterator& v_start, const Iterator& v_end, typename Ring::Element array[]);

	template <typename Ring, typename Vector>
	static inline void copyDenseArrayToSparseVector(Ring& R, typename Ring::Element array[], uint32 arraySize, Vector& v);

	//B <- A^-1 x B <=> TRSM
	//template <typename Matrix, typename Ring>
	//static void reducePivotsByPivots(Ring& R, const Matrix& A, Matrix& B);
	
	template <typename Matrix>
	static void reducePivotsByPivots(Modular<uint16>& R, const Matrix& A, Matrix& B);

	// D <- D - CB <=> GEMM
	//template <typename Matrix, typename Ring>
	//static void reduceNonPivotsByPivots(Ring& R, const Matrix& C, const Matrix& B, Matrix& D);

	template <typename Matrix>
	static void reduceNonPivotsByPivots(Modular<uint16>& R, const Matrix& C, const Matrix& B, Matrix& D);
	
	//template <typename Ring, typename Matrix>
	//static size_t  echelonize(const Ring& R, Matrix& A);
	
	template <typename Matrix>
	static size_t  echelonize(const Modular<uint16>& R, Matrix& A);

	template <typename Matrix, typename Ring>
	static void normalizeRows(Ring& R, Matrix& A);


private:
	//reduce row so that the first entry is equal to unity
	template <typename Ring, typename Vector>
	static void normalizeRow(const Ring& R, Vector& x);

};

#include "matrix-op.C"

#endif /* MATRIX_OP_H_ */
