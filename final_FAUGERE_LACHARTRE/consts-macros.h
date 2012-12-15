/*
 * consts-macros.h
 *
 *  Created on: 30 juil. 2012
 *      Author: martani
 */

#ifndef CONSTS_MACROS_H_
#define CONSTS_MACROS_H_

#include <string>
#include "lela/ring/modular.h"

using namespace LELA;

#ifndef DEFAULT_BLOC_SIZE
#define DEFAULT_BLOC_SIZE 0
#error "must define DEFAULT_BLOC_SIZE"
#endif


//#ifndef DEFAULT_BLOC_HEIGHT
#define DEFAULT_BLOC_HEIGHT DEFAULT_BLOC_SIZE
//#error "must define DEFAULT_BLOC_HEIGHT"
//#endif

//#ifndef DEFAULT_BLOC_WIDTH
#define DEFAULT_BLOC_WIDTH DEFAULT_BLOC_SIZE
//#error "must define DEFAULT_BLOC_WIDTH"
//#endif

//#ifndef HYBRID_REPRESENTATION_THRESHOLD
#define HYBRID_REPRESENTATION_THRESHOLD 0.5f
//#error "must define HYBRID_REPRESENTATION_THRESHOLD"
//#endif

#ifndef NB_ROWS_PER_MULTILINE
#define NB_ROWS_PER_MULTILINE 2
#endif

#if NB_ROWS_PER_MULTILINE != 2
#error "NB_ROWS_PER_MULTILINE must be equal to 2. Other sizes are not supported in this version!"
#endif

//#ifndef NUM_THREADS
//#define NUM_THREADS		2
//#endif

#ifndef NUM_THREADS_OMP
#define NUM_THREADS_OMP		4
#endif

#ifndef NUM_THREADS_OMP_MASTER
#define NUM_THREADS_OMP_MASTER		4
#endif

#ifndef NUM_THREADS_OMP_SLAVES_PER_MASTER
#define NUM_THREADS_OMP_SLAVES_PER_MASTER		2
#endif


#define UNROLL_STEP__64		16
#define UNROLL_STEP__16		16

#if DEFAULT_BLOC_WIDTH % UNROLL_STEP__64 != 0
#error "DEFAULT_BLOC_WIDTH must devide UNROLL_STEP__64 = 16 "
#endif


#ifdef DEFAULT_BLOC_HEIGHT
#if DEFAULT_BLOC_HEIGHT>=256
	typedef uint16 IndexType;
	std::string __IndexType__ = "uint16";
#else
	typedef uint8 IndexType;
	std::string __IndexType__ = "uint8";
#endif
#else
	typedef uint16 IndexType;
	std::string __IndexType__ = "uint16";
#endif


#define __PREFETCH_WRITE	1
#define __PREFETCH_READ		0
#define __PREFETCH_LOCALITY_NO_LOCALITY	0
#define __PREFETCH_LOCALITY_LOW	1
#define __PREFETCH_LOCALITY_MODERATE	2
#define __PREFETCH_LOCALITY_HIGH	3

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

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

#ifdef USE_MUTEX					
#define LOCK(lock)					\
	pthread_mutex_lock(&lock##_mutex_lock);		
#else							
#define LOCK(lock)					\
	pthread_spin_lock(&lock##_spinlock_lock);	
#endif							

#ifdef USE_MUTEX					
#define UNLOCK(lock)					\
	pthread_mutex_unlock(&lock##_mutex_lock);		
#else							
#define UNLOCK(lock)					\
	pthread_spin_unlock(&lock##_spinlock_lock);	
#endif
		
		
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define check_equal_or_raise_exception(a, b)	\
	if (!((a) == (b)))	{							\
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)	\
		 << "at " << __FILE__ << ":" << __LINE__ << " : __func__ " << __func__ 	\
		 << std::endl << #a " must be equal to " #b << std::endl;			\
		 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	\
		throw std::logic_error (#a " must be equal to " #b);				\
	}

#define SHOW_MATRIX_INFO_BLOC(M)			\
		{									\
	std::pair<uint64, double> ___info___bloc_ = MatrixUtils::getMatrixSizeAndDensity(M, true);						\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << " (" << right << setw(5) << M.rowdim() << "x " << left << setw(5) << M.coldim() << ")\t"	\
			<< "- nb elements " << ___info___bloc_.first << " - density: " << ___info___bloc_.second << " % \t bloc dim "	\
			<< M.bloc_height() << " x " << M.bloc_width() << "  ****" << endl;						\
		}

#define SHOW_MATRIX_INFO_SPARSE(M)			\
		{									\
	std::pair<uint64, double> ___info___sparse_ = MatrixUtils::getMatrixSizeAndDensity(M, true);						\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << " (" << right << setw(5) << M.rowdim() << "x " << left << setw(5) << M.coldim() << ")\t"	\
			<< "- nb elements " << ___info___sparse_.first << " - density: " << ___info___sparse_.second << " % \t bloc dim "	\
			<< "  ****" << endl;			\
		}

#define SHOW_MATRIX_INFO_MULTILINE(M)			\
		{									\
	std::pair<uint64, double> ___info___sparse_ = MatrixUtils::getMatrixSizeAndDensity(M, true);						\
	commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)		\
			<< "****\t" << #M << " (" << right << setw(5) << M.rowdim() << "x " << left << setw(5) << M.coldim() << ")\t"	\
			<< "- nb elements " << ___info___sparse_.first << " - density: " << ___info___sparse_.second << " % \t bloc dim "	\
			<< "  ****" << endl;			\
		}


#endif /* CONSTS_MACROS_H_ */
