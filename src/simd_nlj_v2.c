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
    size_t nBlockWidth = 32 / (sizeof(float));
    size_t cntBlock = n / (nBlockWidth * 4);
    size_t cntRem = n % (nBlockWidth * 4);

    //printf("n=%d, vec_len=%d, cntBlock=%d, cntRem=%d\n", n, vec_len, cntBlock, cntRem);

    int *p1 = (int*)malloc(nBlockWidth * sizeof(int));	
    int *p2 = data;
	
    int *q = (int*)malloc(nBlockWidth * sizeof(int));    //将AVX变量上的多个数值合并时所用指针.
    __m256 outVecs, outVecs1, outVecs2, outVecs3; 
    __m256 innerVecs, innerVecs1, innerVecs2, innerVecs3;
    // set rst to zero
    __m256 cmpRst;	
    __m256 rstSum = _mm256_setzero_ps();	

    for (i = 0; i < vec_len; i++) {
        for (k = 0; k < nBlockWidth; k++) {
	    p1[k] = i;
        }	
        outVecs = _mm256_loadu_ps((float*)p1);	
        p2 = data;
        for (j = 0; j < cntBlock; j++) {
            innerVecs = _mm256_loadu_ps((float*)p2);
	    innerVecs1 = _mm256_loadu_ps((float*)p2 + nBlockWidth);
	    innerVecs2 = _mm256_loadu_ps((float*)p2 + nBlockWidth * 2);
	    innerVecs3 = _mm256_loadu_ps((float*)p2 + nBlockWidth * 3);

	    cmpRst = _mm256_cmp_ps(outVecs, innerVecs, _CMP_EQ_OQ);
	    _mm256_storeu_ps((float*)q, cmpRst);
            result += -(int)q[0] + -(int)q[1] + -(int)q[2] + -(int)q[3] + -(int)q[4] + -(int)q[5] + -(int)q[6] + -(int)q[7];
		
	    cmpRst = _mm256_cmp_ps(outVecs, innerVecs1, _CMP_EQ_OQ);
	    _mm256_storeu_ps((float*)q, cmpRst);
            result += -(int)q[0] + -(int)q[1] + -(int)q[2] + -(int)q[3] + -(int)q[4] + -(int)q[5] + -(int)q[6] + -(int)q[7];

	    cmpRst = _mm256_cmp_ps(outVecs, innerVecs2, _CMP_EQ_OQ);
	    _mm256_storeu_ps((float*)q, cmpRst);
            result += -(int)q[0] + -(int)q[1] + -(int)q[2] + -(int)q[3] + -(int)q[4] + -(int)q[5] + -(int)q[6] + -(int)q[7];

	    cmpRst = _mm256_cmp_ps(outVecs, innerVecs3, _CMP_EQ_OQ);
	    _mm256_storeu_ps((float*)q, cmpRst);
            result += -(int)q[0] + -(int)q[1] + -(int)q[2] + -(int)q[3] + -(int)q[4] + -(int)q[5] + -(int)q[6] + -(int)q[7];
	    p2 += nBlockWidth * 4;
        }			
        for (s = 0; s < cntRem; s++) {
            if (i == *p2)
	        result++;
	    p2++;
        }
    }
	
    sink = result; /* So compiler doesn't optimize away the loop */
    //printf("[%ld]",sink);
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
