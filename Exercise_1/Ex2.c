// gcc pc.c -o pc -pthread -lm

#include <sys/time.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#define QUEUESIZE 10
#define LOOP 5000        
#define NUM_PRODUCERS 3   
#define NUM_CONSUMERS 4 

struct workFunction {
  void * (*work)(void *);
  void * arg;
};

typedef struct {
  struct workFunction buf[QUEUESIZE];
  long head, tail;
  int full, empty;
  pthread_mutex_t *mut;
  pthread_cond_t *notFull, *notEmpty;
} queue;

queue *queueInit (void);
void queueDelete (queue *q);
void queueAdd (queue *q, struct workFunction in);
void queueDel (queue *q, struct workFunction *out);

void *producer (void *args);
void *consumer (void *args);

void *calculate_sine(void *arg) {
  double *start_angle = (double *)arg;
  double angle = *start_angle;
  double result;
  
  for (int i = 0; i < 10; i++) {
    result = sin(angle + (i * 0.1));
  }
  
  free(arg); 
  return NULL;
}

int main ()
{
  queue *fifo;
  pthread_t pro[NUM_PRODUCERS];
  pthread_t con[NUM_CONSUMERS];

  fifo = queueInit ();
  if (fifo ==  NULL) {
    fprintf (stderr, "main: Queue Init failed.\n");
    exit (1);
  }

  for (int i = 0; i < NUM_PRODUCERS; i++) {
    pthread_create (&pro[i], NULL, producer, fifo);
  }

  for (int i = 0; i < NUM_CONSUMERS; i++) {
    pthread_create (&con[i], NULL, consumer, fifo);
  }

  for (int i = 0; i < NUM_PRODUCERS; i++) {
    pthread_join (pro[i], NULL);
  }

  printf("Όλοι οι producers ολοκλήρωσαν την εργασία τους.\n");
  

  usleep(100000); 

  queueDelete (fifo);
  return 0;
}

void *producer (void *q)
{
  queue *fifo;
  int i;

  fifo = (queue *)q;

  for (i = 0; i < LOOP; i++) {
    struct workFunction task;
    task.work = calculate_sine;
    
    double *arg = (double *)malloc(sizeof(double));
    *arg = (double)i * 0.05; 
    task.arg = arg;

    pthread_mutex_lock (fifo->mut);
    while (fifo->full) {
      pthread_cond_wait (fifo->notFull, fifo->mut);
    }
    
    queueAdd (fifo, task);
    
    pthread_mutex_unlock (fifo->mut);
    pthread_cond_signal (fifo->notEmpty);
  }
  
  return (NULL);
}

void *consumer (void *q)
{
  queue *fifo;
  struct workFunction task;

  fifo = (queue *)q;

  while (1) {
    pthread_mutex_lock (fifo->mut);
    while (fifo->empty) {
      pthread_cond_wait (fifo->notEmpty, fifo->mut);
    }
    
    queueDel (fifo, &task);
    
    pthread_mutex_unlock (fifo->mut);
    pthread_cond_signal (fifo->notFull);

    task.work(task.arg);
  }
  
  return (NULL);
}


queue *queueInit (void)
{
  queue *q;

  q = (queue *)malloc (sizeof (queue));
  if (q == NULL) return (NULL);

  q->empty = 1;
  q->full = 0;
  q->head = 0;
  q->tail = 0;
  q->mut = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t));
  pthread_mutex_init (q->mut, NULL);
  q->notFull = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
  pthread_cond_init (q->notFull, NULL);
  q->notEmpty = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
  pthread_cond_init (q->notEmpty, NULL);
    
  return (q);
}

void queueDelete (queue *q)
{
  pthread_mutex_destroy (q->mut);
  free (q->mut);    
  pthread_cond_destroy (q->notFull);
  free (q->notFull);
  pthread_cond_destroy (q->notEmpty);
  free (q->notEmpty);
  free (q);
}

void queueAdd (queue *q, struct workFunction in)
{
  q->buf[q->tail] = in;
  q->tail++;

  if (q->tail == QUEUESIZE)
    q->tail = 0;

  if (q->tail == q->head)
    q->full = 1;

  q->empty = 0;

  return;
}

void queueDel (queue *q, struct workFunction *out)
{
  *out = q->buf[q->head];

  q->head++;

  if (q->head == QUEUESIZE)
    q->head = 0;

  if (q->head == q->tail)
    q->empty = 1;

  q->full = 0;

  return;
}