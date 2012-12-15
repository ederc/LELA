/*
 * matrix-util.h
 *
 *  Created on: 4 juil. 2012
 *      Author: martani
 */

#ifndef MATRIX_UTILS_H_
#define MATRIX_UTILS_H_

#include <unistd.h>

#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace LELA;


#ifdef DETAILED_PROFILE_TIMERS
#  define TIMER_DECLARE_(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_RESET_(part) part##_time = 0.0;
#  define TIMER_START_(part) part##_timer.start ()
#  define TIMER_STOP_(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT_(part) \
		commentator.report (Commentator::LEVEL_ALWAYS, TIMING_MEASURE)		\
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
//#  define TIMER_REPORT_(part)
#else
#  define TIMER_DECLARE_(part)
#  define TIMER_RESET_(part)
#  define TIMER_START_(part)
#  define TIMER_STOP_(part)
#  define TIMER_REPORT_(part)
#endif //DETAILED_PROFILE_TIMERS

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

#define check_equal_or_raise_exception(a, b)	\
	if (!((a) == (b)))	{							\
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)	\
		 << "at " << __FILE__ << ":" << __LINE__ << " : __func__ " << __func__ 	\
		 << std::endl << #a " must be equal to " #b << std::endl;			\
		 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	\
		throw std::logic_error (#a " must be equal to " #b);				\
	}


class MatrixUtils {

public:
	template <typename Element, typename Index>
	static void copy(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B);

	template <typename Element, typename Index>
	static void copy(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B);

	///dumps the matrix content as a PBM image. null elements are represented with a white pixel
	///other elements are represented with a black pixel
	template <typename Element, typename Index>
	static void dumpMatrixAsPbmImage(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, const char *outputFileName);

	///dumps the matrix content as a PBM image. null elements are represented with a white pixel
	///other elements are represented with a black pixel
	template <typename Element>
	static void dumpMatrixAsPbmImage(const SparseMatrix<Element>& A, const char *outputFileName);

	template <typename Element, typename Index>
	static bool equal(const Modular<Element>& R, const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, const SparseMatrix<Element>& B);

	template <typename Element, typename Index>
	static bool equal(const Modular<Element>& R, const SparseMatrix<Element>& A, const SparseBlocMatrix<ContiguousBloc<Element, Index> >& B);

	static void show_mem_usage(std::string msg);

	static uint32 loadF4Modulus(const char *fileName);

	template <class Ring>
	static void loadF4Matrix(const Ring &R, SparseMatrix<typename Ring::Element>& A, const char *fileName);

	template<class Ring>
	static void loadF4Matrix__low_memory(const Ring &R, SparseMatrix<typename Ring::Element>& A, const char *fileName);

private:
	/*template <typename Element, typename Index>
	static void sparse_to_bloc_copyDownTopRightLeft(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B);*/
	template <typename Element, typename SparseBlocMatrix_>
	static void sparse_to_bloc_copyDownTopRightLeft(const SparseMatrix<Element>& A, SparseBlocMatrix_& B);

	template <typename Element, typename Index>
	static void sparse_to_bloc_copyTopDowLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B);

	template <typename Element, typename Index>
	static void sparse_to_bloc_copyDownTopLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B);



	template <typename Element, typename Index>
	static void bloc_to_sparse_copyDownTopRightLeft(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B);

	template <typename Element, typename Index>
	static void bloc_to_sparse_copyTopDowLeftRight(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B);

	template <typename Element, typename Index>
	static void bloc_to_sparse_copyDownTopLeftRight(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B);


	static void process_mem_usage(double& vm_usage, double& resident_set);
};

#include "matrix-utils.C"

#endif /* MATRIX_UTIL_HS_ */
