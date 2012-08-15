/*
 * test.c
 *
 *  Created on: 9 août 2012
 *      Author: martani
 */

#include "stdlib.h"
#include "stdio.h"
#include "unistd.h"

#include "thr_pool.h"


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
	fprintf(stdout, "Task 2 [%lu] Done!\n", pthread_self());
	//sleep(1);

	return NULL;
}

int main()
{
	thr_pool_t *pool = thr_pool_create(2, 3, 3, NULL);

	for(int i=0; i<10; ++i)
		thr_pool_queue(pool, first_task, &i);

	thr_pool_wait(pool);

	for(int i=0; i<10; ++i)
		thr_pool_queue(pool, second_task, &i);

	thr_pool_destroy(pool);
}
