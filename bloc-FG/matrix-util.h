/*
 * matrix-util.h
 *
 *  Created on: 12 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_UTIL_H_
#define MATRIX_UTIL_H_

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

#include "bloc-types.h"

using namespace LELA;


#ifndef DETAILED_PROFILE_TIMERS
#define DETAILED_PROFILE_TIMERS
#  define TIMER_DECLARE_(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START_(part) part##_timer.start ()
#  define TIMER_STOP_(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT_(part) \
		commentator.report (Commentator::LEVEL_ALWAYS, TIMING_MEASURE)		\
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
//#  define TIMER_REPORT_(part)
#else
#  define TIMER_DECLARE(part)
#  define TIMER_START(part)
#  define TIMER_STOP(part)
#  define TIMER_REPORT(part)
#endif //DETAILED_PROFILE_TIMERS

class MatrixUtil {

public:
	static void show_mem_usage(std::string msg);

	static uint32 loadF4Modulus(const char *fileName);

	template <class Ring>
	static SparseMatrix<typename Ring::Element> loadF4Matrix(const Ring &R, const char *fileName);

	template <typename Matrix>
	static void dumpMatrixAsPbmImage(const Matrix& A, const char *outputFileName);

	template <typename Element, typename Index>
	static void MatrixUtil::dumpMatrixAsPbmImage(const FG_Types::SparseBlocMatrix<Element, Index>& A, const char *outputFileName);

private:

	static void process_mem_usage(double& vm_usage, double& resident_set);
};

#include "matrix-util.C"

#endif /* MATRIX_UTIL_H_ */
