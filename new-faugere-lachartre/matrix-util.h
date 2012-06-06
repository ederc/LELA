#ifndef MATRIX_UTIL_H_
#define MATRIX_UTIL_H_

/*
 * matrix-util.h
 *
 *  Created on: 30 mai 2012
 *      Author: martani
 */
#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"


using namespace LELA;



#ifndef DETAILED_PROFILE_TIMERS
#define DETAILED_PROFILE_TIMERS
#  define TIMER_DECLARE_(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START_(part) part##_timer.start ()
#  define TIMER_STOP_(part) part##_timer.stop (); part##_time += part##_timer.time ()
/*#  define TIMER_REPORT_(part) \
		commentator.report (Commentator::LEVEL_ALWAYS, TIMING_MEASURE)		\
		<< "Total " #part " time: " << part##_time << "s" << std::endl;*/
#  define TIMER_REPORT_(part)
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

	template <class Ring, typename Matrix>
	static void writeF4MatrixToFile(const Ring &R, const char *fileName, const Matrix& A);

	template <typename Matrix>
	static void dumpMatrixAsPbmImage(const Matrix& A, const char *outputFileName);

	static std::string getOutputFileNameWithExtension(const char *inputfile, const char* prefix, const char *suffix);

	template <typename Matrix, typename Ring>
	static bool verifyMatrixRowsAreUnitary(Ring& R, const Matrix& A, typename Ring::Element& det);

	template <typename Matrix, typename Ring>
	static bool verifyEveryRowEntryIsNoGreaterThanThePreceedingRows(Ring& R, const Matrix& A);

	//verifies the equality of matrices by the hash of their rows (if rows are not in the same order, the matrices are still equal)
	template <typename Matrix, typename Ring>
	static bool matrixEqualUsingHash(Ring& R, const Matrix& A, const Matrix& B);

	template <typename Matrix>
	static void invertMatrixRows(Matrix& A);

	// For each r in A: r <- entry(r)^-1 * r
	template <typename Ring, typename Matrix>
	static void makeRowsUnitary(const Ring& R, Matrix& A);

	template <typename Ring>
	static SparseMatrix<typename Ring::Element> generateIDMatrix(const Ring R, size_t size);

	template <typename Matrix>
	static void freeMatrixMemory(Matrix& A);

	template <typename Matrix>
	static std::pair<uint64, double> getMatrixSizeAndDensity(const Matrix& A);

private:

	static void process_mem_usage(double& vm_usage, double& resident_set);

	template <typename T>
	static std::size_t hasharray(const T arr[], int N);

	template <typename Matrix>
	static void sortMatrixPivotsByIncrementingRowEntry(Matrix& A);

	template <typename Array, typename Matrix>
	static void sortRows_in(Matrix& A, Array pivots);

};

#include "matrix-util.C"

#endif /* MATRIX_UTIL_H_ */
