/*
 * test-Thread-Pool.C
 *
 *  Created on: 9 août 2012
 *      Author: martani
 */

#include <iostream>
#include <stdlib.h>

#include "consts-macros.h"

#include "ThreadPool.h"


//#include "thrpool-0.8/src/TThreadPool.hh"
//
//#ifndef NULL
//#define NULL   ((void *) 0)
//#endif
//
//
//class TMyJob : public ThreadPool::TPool::TJob
//{
//public:
//  TMyJob ( int p ) : ThreadPool::TPool::TJob(p) {}
//
//  void run ( void * arg )
//  {
//		int i = *(int*) arg;
//		fprintf(stdout, "  Thread \t\t\t[%lu]\n", pthread_self());
//		fprintf(stdout, "  handling %d\n", i);
//		sleep(rand() % 3);
//		fprintf(stdout, "  Thread \t\t\t[%lu] Done!\n", pthread_self());
//  }
//};
//
//
//class TMyJob2 : public ThreadPool::TPool::TJob
//{
//public:
//  TMyJob2 ( int p ) : ThreadPool::TPool::TJob(p) {}
//
//  void run ( void * arg )
//  {
//		std::cout << ("Task 2") << std::endl;
//		sleep(1);
//  }
//};


using namespace std;

void* first_task(void *arg)
{
	int i = *(int*) arg;
	fprintf(stdout, "  Thread \t\t\t[%lu]\n", pthread_self());
	fprintf(stdout, "  handling %d\n", i);
	//sleep(rand() % 3);
	fprintf(stdout, "  Thread \t\t\t[%lu] Done!\n", pthread_self());

	return NULL;
}

void* second_task(void *arg)
{
	cout << ("Task 2") << endl;
	//sleep(1);

	return NULL;
}




int main()
{

omp_set_dynamic(0);
#pragma omp parallel num_threads(NUM_THREADS_OMP_MASTER)
{
	cout << "PARALLEL" << endl;
	ThreadPool pool;

	pool.initThreadPool(3);

	for(int i=0; i<1000; ++i)
		pool.queueTask(first_task, &i);

	//pool.waitAllTasks();

	for(int i=0; i<1000; ++i)
		pool.queueTask(second_task, &i);
}
	//ThreadPool::init(3);

//	ThreadPool::TPool *pool = new ThreadPool::TPool (3);
//	TMyJob *job1 = new TMyJob(1);
//	TMyJob *job2 = new TMyJob(2);
//
//	for(int i=0; i<10; ++i)
//	{
//		//TMyJob* job = new TMyJob(i);
//		pool->run(job1, &i, false);
//		pool->run(job2, &i, false);
//	}
//
//	std::cout << "SECOND STEP" << endl;
//
//	for(int i=0; i<10; ++i)
//	{
//		TMyJob* job = new TMyJob(i);
//		ThreadPool::run(job, &i, false);
//	}
//
//	for(int i=0; i<10; ++i)
//	{
//		TMyJob2* job = new TMyJob2(i);
//		ThreadPool::run(job, &i, false);
//	}

	cout << "OK" << endl;

	return 0;
}
