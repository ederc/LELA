/*
 * structured-gauss-lib.h
 * Copyright 2012 Martani Fayssal (LIP6 / UPMC University Paris06)
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 * 
 *  An implementation of structured Gaussian elimination
 */

#ifndef STRUCTURED_GAUSS_LIB_H_
#define STRUCTURED_GAUSS_LIB_H_

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace LELA;

class StructuredGauss {

public:

	template <typename Ring, typename Matrix>
	static size_t  echelonize_reduced(const Ring& R, Matrix& A);

private:

	template <typename Vector>
	static inline bool isRowNull(Vector& v);

	template <typename Vector>
	static inline void swapRows(Vector&v, Vector& w);

	template <typename Ring, typename Vector>
	static void normalizeRow(const Ring& R, Vector& x);

	template <typename Array, typename Matrix>
	static void sortRows(Matrix& A, Array pivots);

	template <typename Ring>
	static inline void razArray_in(const Ring& R, typename Ring::Element arr[], uint32 arrSize);

	template <typename Ring, typename Vector>
	static inline void copySparseVectorToDenseArray_in(Ring R, const Vector& v, typename Ring::Element array[]);

	template <typename Ring, typename Iterator>
	static inline void copySparseVectorToDenseArray_in(Ring R, const Iterator& v_start, const Iterator& v_end, typename Ring::Element array[]);

	template <typename Ring, typename Vector>
	static inline void copyDenseArrayToSparseVector_in(Ring& R, typename Ring::Element array[], uint32 arraySize, Vector& v);

	template<typename Context, typename Element, typename Vector>
	static inline typename Vector::value_type::first_type head_generic(Context& ctx, Element& a, const Vector& v);

};


#include "structured-gauss-lib.C"

#endif /* STRUCTURED_GAUSS_LIB_H_ */
