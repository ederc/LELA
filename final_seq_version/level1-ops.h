/*
 * level1-ops.h
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */

#ifndef LEVEL1_OPS_H_
#define LEVEL1_OPS_H_

#include "consts-macros.h"
#include "types.h"

class Level1Ops
{
public:
	/**
	 * Given a bloc of dense rows (an array of arrays), Zeros its memory
	 */
	template <typename Element>
	static inline void memsetToZero(Element** arr);

	template <typename Element>
	static inline void memsetToZero(Element** arr, const uint32 nb_lines, const uint32 line_size);


	template<typename Ring, typename Index, typename DoubleFlatElement>
	static void copySparseBlocToDenseBlocArray(const Ring& R,
		const SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		DoubleFlatElement** arr);


	template<typename Ring, typename DoubleFlatElement, typename Index>
	static void copyDenseBlocArrayToSparseBloc(const Ring& R,
		DoubleFlatElement** arr,
		SparseMultilineBloc<typename Ring::Element, Index>& bloc,
		bool reduce_in_Ring = true);


	template <typename Ring, typename DoubleFlatElement>
	static inline void reduceDenseArrayModulo(const Ring& R, DoubleFlatElement* arr);


	template <typename Element, typename Index>
	static inline long headMultiLineVector(const MultiLineVector<Element, Index>& v,
			const uint16 line_idx, uint16& a, uint32& head_idx);

	template <typename Element, typename Index>
	static long headMultiLineVectorHybrid(const MultiLineVector<Element, Index>& v,
			const uint16 line_idx, uint16& a, uint32& head_idx, const size_t SIZE_DENSE_VECTOR);

	template<typename Ring, typename Index>
	static inline void normalizeMultiLineVector(const Ring& R,
			MultiLineVector<typename Ring::Element, Index>& v);


	template <typename Element, typename Index, typename Element2>
	static inline void copyMultiLineVectorToDenseArray(const MultiLineVector<Element, Index>& v,
		Element2 *arr1, Element2 *arr2, size_t arr_size);

	template<typename Ring, typename DoubleFlatElement>
	static long headDenseArray(const Ring& R, DoubleFlatElement arr[],
		const size_t size, typename Ring::Element& a);


	template <typename Ring, typename DoubleFlatElement>
	static long normalizeDenseArray(const Ring& R, DoubleFlatElement arr[], const size_t size);


	template<typename Ring, typename Index, typename DoubleFlatElement>
	static void copyDenseArraysToMultilineVector(const Ring& R,
			const DoubleFlatElement *arr1,
			const DoubleFlatElement *arr2,
			const uint32 size,
			MultiLineVector<typename Ring::Element, Index>& v);

	template<typename Ring, typename Index, typename DoubleFlatElement>
	static void copyDenseArraysToMultilineVectorHybrid(const Ring& R,
		const DoubleFlatElement *arr1,
		const DoubleFlatElement *arr2,
		const uint32 size,
		MultiLineVector<typename Ring::Element, Index>& v);

private:
	Level1Ops() {}
	Level1Ops(const Level1Ops& other) {}
};


#include "level1-ops.C"

#endif /* LEVEL1_OPS_H_ */
