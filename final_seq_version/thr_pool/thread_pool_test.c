/**
 * threadpool_test.c, copyright 2001 Steve Gribble
 * Just a regression test for the threadpool code.
 */

//gcc threadpool.h threadpool.c thread_pool_test.c -lpthread

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>
#include "threadpool.h"
#include "external.h"

//extern int errno;
//extern int add_memory_access(memaddr_t addr, line_type_e type, pthread_t thread);

void dispatch_to_me(void *arg)
{
	int seconds = (int) arg;
	fprintf(stdout, "  in dispatch %d\n", seconds);
	fprintf(stdout, "  done dispatch %d\n", seconds);
}

void handle(void *arg)
{
	int i = (int) arg;
	fprintf(stdout, "  Thread \t\t\t[%lu]\n", pthread_self()%10);
	fprintf(stdout, "  handling %d\n", i);
	sleep(rand() % 3);
	fprintf(stdout, "  Thread \t\t\t[%lu] Done!\n", pthread_self()%10);
}


int main(int argc, char **argv)
{
	threadpool tp;
	int i = 0;
	//construct_locksmith(122);

	/*
	 locksmith_add_memory_access(pthread_self(),2,3,4,5);
	 locksmith_add_lock_to_thread(1);
	 */

	tp = create_threadpool(3);
	for (; i < 5; i++)
	{
		dispatch(tp, handle, (void *) i);
	}
	fprintf(stdout, "**main** done first, waiting tasks\n");

	wait_all_tasks(tp);

	fprintf(stdout, "**main** after wait\n");

	wait_all_tasks(tp);

	sleep(1);

	fprintf(stdout, "Starting second\n\n");
	for (i = 0; i < 10; i++)
	{
		dispatch(tp, dispatch_to_me, (void *) i);
	}

	fprintf(stdout, "**main done second\n");
	//destroy_locksmith();
	destroy_threadpool(tp);
	sleep(5);
	exit(-1);
}
