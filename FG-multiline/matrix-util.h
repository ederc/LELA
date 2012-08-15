/*
 * matrix-util.h
 *
 *  Created on: 19 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_UTIL_H_
#define MATRIX_UTIL_H_

using namespace LELA;
using namespace std;

template<typename Element>
static void dumpMatrixAsPbmImage(const SparseMultilineMatrix<Element>& A, const char *outputFileName);

template<typename Matrix, typename Ring>
static bool equal(const Ring& R, const SparseMultilineMatrix<typename Ring::Element>& A, const Matrix& B);

template <typename Matrix, typename Ring>
static bool equal_reverse(const Ring& R, const SparseMultilineMatrix<typename Ring::Element>& A, const Matrix& B);


#include "matrix-util.C"

#endif /* MATRIX_UTIL_H_ */
