/* $begin mountainmain */
#include <stdlib.h>
#include <stdio.h>
#include "Rdtsc.h"              /* startTimer, stopTimer */
#include <sys/time.h>           /* gettimeofday */
#include "x86intrin.h"              /* startTimer, stopTimer */

#define L1 (1<<15)    /* Working set size for L1 cache 32KB */
#define L2 (1<<18)    /* Working set size for L2 cache 256KB */
#define L3 (1<<20)*2.5    /* Working set size for L3 cache 2.5MB */
#define LLC (1<<20)*55    /* Working set size for LLC cache 55MB */
//#define MAXELEMS 60000000 
#define MAXELEMS 6000 
#define random(x) (rand()%x)
void init_data(int *data, int n, int cardinality);
void test(int *data, int *vector, int n,int vec_len);
long run(int *data, int *vector, int n,int vec_len);
/* $begin mountainmain */
int main()
{
    int * data=(int *)malloc(sizeof(int) * MAXELEMS);      /* foreign key column*/
    printf("SIMD nest loop join(cycles)\n");
    printf("\t");
    int step;
    for (step = 10; step <= 100; step+=10)
	printf("%d%\t", step);
    printf("\n");
	/*SIMD nest loop join for L1 cache sets*/
    printf("L1\t");
    for (step = 10; step <= 100; step+=10){
	int i,vector_len=L1/4*step/100;
	init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
	int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
	for (i=0;i<vector_len;i++)
	    vector[i]=1;
	    printf("%ld\t", run(data,vector,MAXELEMS, vector_len));
	}
	printf("\n");
	/*SIMD nest loop join for L2 cache sets*/
	printf("L2\t");
        for (step = 10; step <= 100; step+=10){
    	    int i,vector_len=L2/4*step/100;
	    init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
    	    int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
	    for (i=0;i<vector_len;i++)
	        vector[i]=1;
	    printf("%ld\t", run(data,vector,MAXELEMS, vector_len));
	}
	printf("\n");
	/*SIMD nest loop join for L3 cache sets*/
	printf("L3\t");
        for (step = 10; step <= 100; step+=10){
	    int i,vector_len=L3/4*step/100;
	    init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
	    int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
	    for (i=0;i<vector_len;i++)
		vector[i]=1;
	    printf("%ld\t", run(data,vector,MAXELEMS,vector_len));
	}
	printf("\n");
	/*SIMD nest loop join for LLC cache sets*/
	printf("LLC\t");
        for (step = 10; step <= 100; step+=10){
	    int i,vector_len=LLC/4*step/100;
	    init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
	    int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
	    for (i=0;i<vector_len;i++)
	        vector[i]=1;
	    printf("%ld\t", run(data,vector,MAXELEMS,vector_len));
	}
	printf("\n");
    exit(0);
}
/* init_data - initializes the array */
void init_data(int *data, int n, int cardinality)
{
    int i;
    for (i = 0; i < n; i++)
	data[i] = random(cardinality);
}
void test(int *data, int *vector, int n, int vec_len) /* The test function */
{
    int i, j, k, s;
    int result = 0; 
    volatile int sink=0; 

    int partBlock = L1 / (2 * sizeof(int));
    int blockNum = vec_len / partBlock;
    int remNum = vec_len % partBlock;

    for (k = 0; k < blockNum; k++) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < partBlock ; i++) {
                if (data[j] == (i + k * partBlock))
                    result++;
            }
        }
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < remNum; i++) {
            if (data[j] == (i + k * partBlock)) 
                result++;
        }
    }
    sink = result; /* So compiler doesn't optimize away the loop */
//  printf("[%ld]",sink);
}

long run(int *data, int *vector, int n, int vec_len)
{   
    long cycles;
    uint64_t timer;
    test(data, vector,n,vec_len);                     /* warm up the cache */       //line:mem:warmup
    startTimer(&timer);
    test(data, vector,n,vec_len);                     /* test for cache locality */    
    stopTimer(&timer); 
    cycles = timer;  
    return cycles; 
}
/* $end mountainfuns */
