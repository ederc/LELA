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
#include "types.h"

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
//#define MIN(a, b) (((a) < (b)) ? (a) : (b))
//#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define check_equal_or_raise_exception(a, b)	\
	if (!((a) == (b)))	{							\
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)	\
		 << "at " << __FILE__ << ":" << __LINE__ << " : __func__ " << __func__ 	\
		 << std::endl << #a " must be equal to " #b << std::endl;			\
		 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	\
		throw std::logic_error (#a " must be equal to " #b);				\
	}

/*#ifdef DETAILED_PROFILE_TIMERS
 #  define TIMER_DECLARE__(timer)	\
	timeval timer##_start__, timer##_end__; double timer##_diff__ = 0.0;

 #  define TIMER_START__(timer)	\
	 gettimeofday(&timer##_start__, NULL);

 #  define TIMER_STOP__(timer)	\
	 gettimeofday(&timer##_end__, NULL);		\
	 timer##_diff__ += (timer##_end__.tv_sec - timer##_start__.tv_sec) * 1000.0;	\
	 timer##_diff__ += (timer##_end__.tv_usec - timer##_start__.tv_usec) / 1000.0;

 #  define TIMER_REPORT__(timer) \
		commentator.report (Commentator::LEVEL_ALWAYS, TIMING_MEASURE)		\
		<< "Total " #timer " time: " << timer##_diff__ << "ms." << std::endl;

 #  define TIMER_RESET__(timer) 	timer##_diff__ = 0.0;
 #else
 #  define TIMER_DECLARE__(part)
 #  define TIMER_RESET__(part)
 #  define TIMER_START__(part)
 #  define TIMER_STOP__(part)
 #  define TIMER_REPORT__(part)
 #endif //DETAILED_PROFILE_TIMERS*/

class MatrixUtils
{

public:
	template<typename Element, typename Index>
	static void copy(const SparseMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B);

	//WRANING: these functions have side effects on entry matrices, generaly the data is changed
	//from haybrid to sparse a vice versa, or destructed when specified
	template<typename Element, typename Index>
	static void copy(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseMatrix<Element>& B, bool destruct_original = false);

	template<typename Element, typename Index>
	static void transformSparseToHybridRows(
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A);

	template<typename Element, typename Index>
	static void transformHybridToSparseRows(
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A);

	///dumps the matrix content as a PBM image. null elements are represented with a white pixel
	///other elements are represented with a black pixel
	template<typename Element, typename Index>
	static void dumpMatrixAsPbmImage(
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			const char *outputFileName);

	///dumps the matrix content as a PBM image. null elements are represented with a white pixel
	///other elements are represented with a black pixel
	template<typename Element>
	static void dumpMatrixAsPbmImage(const SparseMatrix<Element>& A,
			const char *outputFileName);

	template<typename Element, typename Index>
	static bool equal(const Modular<Element>& R,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			const SparseMatrix<Element>& B);

	template<typename Element, typename Index>
	static bool equal(const Modular<Element>& R, const SparseMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B);

	static void show_mem_usage(std::string msg);

	static uint32 loadF4Modulus(const char *fileName);

	template<class Ring>
	static void loadF4Matrix(const Ring &R,
			SparseMatrix<typename Ring::Element>& A, const char *fileName);

	template<class Ring>
	static void loadF4Matrix__low_memory(const Ring &R,
			SparseMatrix<typename Ring::Element>& A, const char *fileName);

	template<class Ring>
	static void loadF4Matrix__low_memory_syscall_no_checks(const Ring &R,
			SparseMatrix<typename Ring::Element>& A, const char *fileName);

	template<typename Element>
	static std::pair<uint64, double> getMatrixSizeAndDensity(
			const SparseMatrix<Element>& A, bool exact);

	template<typename Element, typename Index>
	static std::pair<uint64, double> getMatrixSizeAndDensity(
			const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			bool exact = false);

	template <typename Element>
	static std::pair<uint64, double> getMatrixSizeAndDensity(
			const SparseMultilineMatrix<Element>& A, bool exact);

	template<typename Element, typename Index>
	static void dumpMatrixAsPbmImage(
			const SparseMultilineMatrix<Element, Index>& A,
			const char *outputFileName);

	template<typename Element, typename Index>
	static void transformHybridMatrixToSparseMatrix(
			SparseMultilineMatrix<Element, Index>& A);

	template<typename Element, typename Index>
	static bool equal(const SparseMultilineMatrix<Element, Index>& A,
			const SparseMultilineMatrix<Element, Index>& B);

	template<typename Element, typename Index>
	static void copyBlocToMultilineMatrix_DownTop_RightLeft(
			const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& inMatrix,
			SparseMultilineMatrix<Element>& outMatrix);

	template <typename Matrix>
	static void invertMatrixRows(Matrix& A);

private:
	/*template <typename Element, typename Index>
	 static void sparse_to_bloc_copyDownTopRightLeft(const SparseMatrix<Element>& A, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B);*/
	template<typename Element, typename SparseBlocMatrix_>
	static void sparse_to_bloc_copyDownTopRightLeft(
			const SparseMatrix<Element>& A, SparseBlocMatrix_& B);

	template<typename Element, typename Index>
	static void sparse_to_bloc_copyTopDowLeftRight(
			const SparseMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B);

	template<typename Element, typename Index>
	static void sparse_to_bloc_copyDownTopLeftRight(
			const SparseMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B);

	template<typename Element, typename Index>
	static void bloc_to_sparse_copyDownTopRightLeft(
			const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseMatrix<Element>& B);

	template<typename Element, typename Index>
	static void bloc_to_sparse_copyTopDowLeftRight(
			const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseMatrix<Element>& B);

	template<typename Element, typename Index>
	static void bloc_to_sparse_copyDownTopLeftRight(
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseMatrix<Element>& B, bool destruct_original = false);

	static void process_mem_usage(double& vm_usage, double& resident_set);
};

#include "matrix-utils.C"

#endif /* MATRIX_UTIL_HS_ */
