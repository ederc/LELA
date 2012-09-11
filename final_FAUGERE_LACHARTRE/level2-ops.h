/*
 * level2-ops.h
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */


#ifndef LEVEL2_OPS_H_
#define LEVEL2_OPS_H_

#include "consts-macros.h"
#include "types.h"

class Level2Ops
{
public:

	static void DenseScalMulSub__one_row__array_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const uint64 *arr_source,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));

	static void DenseScalMulSub__two_rows__array_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const uint32 av1_col2,
			const uint32 av2_col2,
			const uint64 *arr_source1,
			const uint64 *arr_source2,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));

	template <typename Index>
	static void SparseScalMulSub__one_row__vect_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const MultiLineVector<uint16, Index>& v,
			const uint32 line,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));

	template <typename Index>
	static void DenseScalMulSub__one_row__vect_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const MultiLineVector<uint16, Index>& v,
			const uint32 line,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));


	template <typename Index>
	static void DenseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const uint32 av1_col2,
			const uint32 av2_col2,
			const MultiLineVector<uint16, Index>& v,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));

	template <typename Index>
	static void SparseScalMulSub__two_rows__vect_array(const uint32 av1_col1,
			const uint32 av2_col1,
			const uint32 av1_col2,
			const uint32 av2_col2,
			const MultiLineVector<uint16, Index>& v,
			uint64 *arr1,
			uint64 *arr2) __attribute__((noinline));


	template <typename Index>
	static void reduceBlocByRectangularBloc(const Modular<uint16>& R,
			const SparseMultilineBloc<uint16, Index>& bloc_A,
			const SparseMultilineBloc<uint16, Index>& bloc_B,
			uint64 **Bloc_acc, bool invert_scalars = true) __attribute__((noinline));


	/**
	 * Reduce the rows inside the bloc by themselves
	 */
	template<typename Index>
	static void reduceBlocByTriangularBloc(const Modular<uint16>& R,
			const SparseMultilineBloc<uint16, Index>& bloc_A,
			uint64** Bloc_acc, bool invert_scalars = true) __attribute__((noinline));

	template<typename Index>
	static void reduceBlocByRectangularBloc_C(const Modular<uint16>& R,
			const SparseMultilineBloc<uint16, Index>& bloc_C,
			const SparseMultilineBloc<uint16, Index>& bloc_A, uint64 **bloc_dense);


	template <typename Index>
	static void reduceBlocByTriangularBloc_C(const Modular<uint16>& R,
			const SparseMultilineBloc<uint16, Index>& bloc_A,
			uint64 **bloc_dense);


	template <typename Index>
	static void DenseScalMulSub__one_row__vect_array__variable_size(const uint32 av1_col1,
			const uint32 av2_col1,
			const MultiLineVector<uint16, Index>& v,
			const uint32 line,
			uint64 *arr1,
			uint64 *arr2,
			const int start_offset = 0);


	template <typename Index>
	static void DenseScalMulSub__two_rows__vect_array__variable_size(const uint32 av1_col1,
			const uint32 av2_col1,
			const uint32 av1_col2,
			const uint32 av2_col2,
			const MultiLineVector<uint16, Index>& v,
			uint64 *arr1,
			uint64 *arr2,
			const int start_offset = 0);

private:
	Level2Ops() {}
	Level2Ops(const Level2Ops& other) {}
};



#include "level2-ops.C"


#endif /* LEVEL2_OPS_H_ */

