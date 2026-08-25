/*     SPAX: a real mutex for the U6 patch map.
 *
 *     ccx parallelises with pthreads and compiles the Fortran WITHOUT
 *     -fopenmp (see FFLAGS), so every !$omp directive in a .f file is an
 *     inert comment.  u6patch.f guarded its one-time node->element map
 *     build with !$omp critical, which therefore synchronised nothing:
 *     all assembly threads allocated and filled nstart/nlist at once.
 *     Measured on the sphere cell -- equilibrium_gap 8.7e-08 with one
 *     assembly thread, 0.5..1.4 and non-reproducible with eight.
 *
 *     A pthread mutex is what this threading model actually needs, and it
 *     carries the memory barrier that makes the finished map visible to
 *     every other thread.
 */

#include <pthread.h>

static pthread_mutex_t u6mutex = PTHREAD_MUTEX_INITIALIZER;

void u6lock_(void)   { pthread_mutex_lock(&u6mutex); }
void u6unlock_(void) { pthread_mutex_unlock(&u6mutex); }
