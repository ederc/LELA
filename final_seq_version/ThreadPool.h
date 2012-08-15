/*
 * ThreadPool.h
 *
 *  Created on: 9 août 2012
 *      Author: martani
 */

#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#define MAXT_IN_POOL 200

#include <pthread.h>
#include <malloc.h>

#ifndef NULL
#define NULL   ((void *) 0)
#endif

class ThreadPool {

public:

	void initThreadPool (int num_threads_in_pool);

	void destroyThreadPool ();

	void queueTask (void *(*func)(void *), void* args);

	void waitAllTasks ();

	void* spin ();

	~ThreadPool ()
	{
		destroyThreadPool ();
	}

private:

	typedef struct job
	{
		void *(*func)(void *);
		void * arg;
		struct job* next;
	} job;

	int num_threads;			//number of active threads
	int qsize; 					//number in the queue
	pthread_t *threads; 		//pointer to threads
	job* qhead; 				//queue head pointer
	job* qtail; 				//queue tail pointer
	pthread_mutex_t qlock; 		//lock on the queue list
	pthread_cond_t q_not_empty; //non empty and empty condidtion vairiables
	pthread_cond_t q_empty;
	int shutdown;
	int dont_accept;

	pthread_cond_t task_finished;
	int nb_tasks_running;

};

// _run_thread(void*) is supplied which does nothing but call the member function run():
extern "C" void *
_run_spin ( void * arg )
{  if (arg != NULL)
	 ((ThreadPool*) arg)->spin();

   return NULL;
}

void* ThreadPool::spin()
{
	job* cur; //The q element


	while (1)
	{
		//pool->qsize = pool->qsize;
		pthread_mutex_lock(&(this->qlock)); //get the q lock.

		while (this->qsize == 0)
		{ //if the size is 0 then wait.
			if (this->shutdown)
			{
				pthread_mutex_unlock(&(this->qlock));
				pthread_cond_signal(&(this->task_finished));
				pthread_exit(NULL);
			}
			//wait until the condition says its no emtpy and give up the lock.
			pthread_mutex_unlock(&(this->qlock)); //get the qlock.
			pthread_cond_wait(&(this->q_not_empty), &(this->qlock));

			//check to see if in shutdown mode.
			if (this->shutdown)
			{
				pthread_mutex_unlock(&(this->qlock));
				pthread_cond_signal(&(this->task_finished));
				pthread_exit(NULL);
			}
		}

		cur = this->qhead; //set the cur variable.
		this->nb_tasks_running++;

		this->qsize--; //decriment the size.

		if (this->qsize == 0)
		{
			this->qhead = NULL;
			this->qtail = NULL;
		}
		else
		{
			this->qhead = cur->next;
		}

		if (this->qsize == 0 && !this->shutdown)
		{
			//the q is empty again, now signal that its empty.
			pthread_cond_signal(&(this->q_empty));
		}
		pthread_mutex_unlock(&(this->qlock));

		(cur->func)(cur->arg); //actually do work.
		free(cur); //free the work storage.

		pthread_mutex_lock(&(this->qlock));
		this->nb_tasks_running--;
		pthread_mutex_unlock(&(this->qlock));

		pthread_cond_signal(&(this->task_finished));
	}

	return NULL;
}


void ThreadPool::initThreadPool(int num_threads_in_pool)
{
	int i;

	if ((num_threads_in_pool <= 0) || (num_threads_in_pool > MAXT_IN_POOL))
		return;

//	pool = (_threadpool *) malloc(sizeof(_threadpool));
//	if (pool == NULL)
//	{
//		fprintf(stderr, "Out of memory creating a new threadpool!\n");
//		return NULL;
//	}

	this->threads = (pthread_t*) malloc(sizeof(pthread_t) * num_threads_in_pool);

	//TODO: throw exceptions instead of error messages

	if (!this->threads)
	{
		fprintf(stderr, "Out of memory creating a new threadpool!\n");
		return;
	}

	this->num_threads = num_threads_in_pool; //set up structure members
	this->nb_tasks_running = 0;
	this->qsize = 0;
	this->qhead = NULL;
	this->qtail = NULL;
	this->shutdown = 0;
	this->dont_accept = 0;

	//initialize mutex and condition variables.
	if (pthread_mutex_init(&(this->qlock), NULL))
	{
		fprintf(stderr, "Mutex initiation error!\n");
		return;
	}
	if (pthread_cond_init(&(this->q_empty), NULL))
	{
		fprintf(stderr, "CV initiation error!\n");
		return;
	}
	if (pthread_cond_init(&(this->q_not_empty), NULL))
	{
		fprintf(stderr, "CV initiation error!\n");
		return;
	}

	//XXX: modif
	if (pthread_cond_init(&(this->task_finished), NULL))
	{
		fprintf(stderr, "CV initiation error!\n");
		return;
	}

	//make threads

	for (i = 0; i < num_threads_in_pool; i++)
	{
		if (pthread_create(&(this->threads[i]), NULL, _run_spin, this))
		{
			fprintf(stderr, "Thread initiation error!\n");
			return ;
		}
	}
}



void ThreadPool::queueTask(void *(*func)(void *), void* args)
{
	job * cur;

	//make a work queue element.
	cur = (job*) malloc(sizeof(job));
	if (cur == NULL)
	{
		fprintf(stderr, "Out of memory creating a work struct!\n");
		return;
	}

	cur->func = func;
	cur->arg = args;
	cur->next = NULL;

	pthread_mutex_lock(&(this->qlock));

	if (this->dont_accept)
	{ //Just incase someone is trying to queue more
		free(cur); //work structs.
		return;
	}
	if (this->qsize == 0)
	{
		this->qhead = cur; //set to only one
		this->qtail = cur;
		pthread_cond_signal(&(this->q_not_empty)); //I am not empty.
	}
	else
	{
		this->qtail->next = cur; //add to end;
		this->qtail = cur;
	}
	this->qsize++;
	//this->nb_tasks_running++;

	pthread_mutex_unlock(&(this->qlock)); //unlock the queue.
}


void ThreadPool::destroyThreadPool()
{
	int i = 0;
	void *status;

	pthread_mutex_lock(&(this->qlock));
	this->dont_accept = 1;
	while (this->qsize != 0)
	{
		pthread_cond_wait(&(this->q_empty), &(this->qlock)); //wait until the q is empty.
	}

	this->shutdown = 1; //allow shutdown
	pthread_cond_broadcast(&(this->q_not_empty)); //allow code to return NULL;
	pthread_mutex_unlock(&(this->qlock));

	//kill everything.
	for (; i < this->num_threads; i++)
	{
		//		pthread_cond_broadcast(&(pool->q_not_empty));
		//allowcode to return NULL;/
		pthread_join(this->threads[i], &status);
	}


	//this->waitAllTasks();

	free(this->threads);

	pthread_mutex_destroy(&(this->qlock));
	pthread_cond_destroy(&(this->q_empty));
	pthread_cond_destroy(&(this->q_not_empty));

	//XXX
	pthread_cond_destroy(&(this->task_finished));

}

void ThreadPool::waitAllTasks()
{
	fprintf(stderr, "[thread pool %lu] Waiting for tasks to finish - tasks running %d\n", pthread_self(), this->nb_tasks_running);

	pthread_mutex_lock(&(this->qlock)); //get the q lock.

		while (this->nb_tasks_running > 0)
		{
			//wait until the condition says its emtpy and give up the lock.
			pthread_mutex_unlock(&(this->qlock)); //get the qlock.

			fprintf(stderr, "[thread pool %lu] Waiting on condition - tasks running %d\n", pthread_self(), this->nb_tasks_running);
			pthread_cond_wait(&(this->task_finished), &(this->qlock));
		}

	//this->nb_tasks_running = 0;
	pthread_mutex_unlock(&(this->qlock));

	fprintf(stderr, "[thread pool %lu] All tasks done\n", pthread_self());
}

#endif /* THREADPOOL_H_ */
