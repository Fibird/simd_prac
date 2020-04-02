/* $begin mountainmain */
#include <stdlib.h>
#include <stdio.h>
#include "Rdtsc.h"              /* startTimer, stopTimer */
#include <sys/time.h>           /* gettimeofday */
#include "x86intrin.h"
#include <immintrin.h>

#define L1 (1 << 15)    /* Working set size for L1 cache 32KB */
#define L2 (1 << 18)    /* Working set size for L2 cache 256KB */
#define L3 ((1 << 20) * 5 / 2)    /* Working set size for L3 cache 2.5MB */
#define LLC ((1 << 20) * 55)    /* Working set size for LLC cache 55MB */
#define MAXELEMS 6000
#define random(x) (rand()%(x))

void init_data(int *data, int n, int cardinality);
void test(int *data, int *vector, int n, int vec_len);
void test_SIMD(int *data, int *vector, int n, int vec_len);
void test_SIMD_partition(int *data, int *vector, int n, int vec_len);
double run(int *data, int *vector, int n, int vec_len);

/* $begin mountainmain */
int main()
{
	int * data=(int *)malloc(sizeof(int) * MAXELEMS);      /* foreign key column*/
    printf("nest loop join(cycles)\n");
    printf("\t");
	int step;
    for (step = 10; step <= 100; step+=10) {
		printf("%d\t", step);
	}
    printf("\n");
	/*vector join for L1 cache sets*/
	printf("L1\t");
    for (step = 10; step <= 100; step+=10){
		int i,vector_len=L1/4*step/100;
		init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
		int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
		for (i=0;i<vector_len;i++) {
			vector[i]=i;
		}
		printf("%lf\t", run(data,vector,MAXELEMS, vector_len) );
	}
	printf("\n");
	/*vector join for L2 cache sets*/
	printf("L2\t");
    for (step = 10; step <= 100; step+=10){
		int i,vector_len=L2/4*step/100;
		init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
		int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
		for (i=0;i<vector_len;i++) {
			vector[i]=i;
		}
		printf("%lf\t", run(data,vector,MAXELEMS, vector_len));
	}
	printf("\n");
	/*vector join for L3 cache sets*/
	printf("L3\t");
    for (step = 10; step <= 100; step+=10){
		int i,vector_len=L3/4*step/100;
		init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
		int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
		for (i=0;i<vector_len;i++) {
			vector[i]=i;
		}
		printf("%lf\t", run(data,vector,MAXELEMS, vector_len) );
	}
	printf("\n");
	/*vector join for LLC cache sets*/
	printf("LLC\t");
    for (step = 10; step <= 100; step+=10){
		int i,vector_len=LLC/4*step/100;
		init_data(data, MAXELEMS,vector_len); /* Initialize foreign key elements in data */
		int * vector=(int *)malloc(sizeof(int) * vector_len);   /* Initialize primary key vector elements in vector */
		for (i=0;i<vector_len;i++) {
			vector[i]=i;
		}
		printf("%lf\t", run(data,vector,MAXELEMS, vector_len));
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
    int i, j;
    int result = 0; 
    volatile int sink=0; 
    for (i = 0; i < n; i++) {
        for (j = 0; j < vec_len; ++j) {
            if (data[i] == vector[j]) {
                ++result;
            }
        }
    }
    sink = result; /* So compiler doesn't optimize away the loop */
	// printf("[%d]",sink);
}

void test_SIMD(int *data, int *vector, int n, int vec_len) {
	__m256i yidOuter;
	__m256i yidInner;
	__m256i yidResult = _mm256_set1_epi32(0);
	__m256i yidTmp;

	volatile int sink = 0;
	int i, j;
	int ret = 0;
	int q[8];
	int bound = vec_len / 8 * 8;
	void *p;
	for (i = 0; i < n; ++i) {
		yidOuter = _mm256_set1_epi32(data[i]);
		p = vector;
		for (j = 0; j < bound; j+= 8) {
			yidInner = _mm256_loadu_si256(p);
			yidTmp = _mm256_cmpeq_epi32(yidOuter, yidInner);
			yidResult = _mm256_sub_epi32(yidResult, yidTmp);
			p += 32;
		}
		for (; j < vec_len; ++j) {
			if (data[i] == vector[j]) {
                ++ret;
            }
		}
	}
	_mm256_storeu_si256((void *)q, yidResult);
	ret += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
	sink = ret;
	// printf("[%d]", sink);
}

void test_SIMD_partition(int *data, int *vector, int n, int vec_len) {
	__m256i yidOuter;
	__m256i yidInner;
	__m256i yidResult = _mm256_set1_epi32(0);
	__m256i yidTmp;

	volatile int sink = 0;
	int i, j, k;
	int ret = 0;
	int q[8];
	int partition = L1/ 2 / 4;
	void *p;
	for (k = 0; k < vec_len; k += partition) {
		int end = k + partition < vec_len ? k + partition : vec_len;
		int bound = (end - k) / 8 * 8 + k;
		for (i = 0; i < n; ++i) {
			yidOuter = _mm256_set1_epi32(data[i]);
			p = vector + k;
			for (j = k; j < bound; j+= 8) {
				yidInner = _mm256_loadu_si256(p);
				yidTmp = _mm256_cmpeq_epi32(yidOuter, yidInner);
				yidResult = _mm256_sub_epi32(yidResult, yidTmp);
				p += 32;
			}
			for (; j < end; ++j) {
				if (data[i] == vector[j]) {
					++ret;
				}
			}
		}
	}
	_mm256_storeu_si256((void *)q, yidResult);
	ret += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
	sink = ret;
	// printf("[%d]", sink);
}

double run(int *data, int *vector, int n, int vec_len)
{   
    double cycles;
	uint64_t timer;
    test_SIMD(data, vector, n, vec_len);                     /* warm up the cache */       //line:mem:warmup
	startTimer(&timer);
    test_SIMD(data, vector, n, vec_len);                     /* test for cache locality */    
	stopTimer(&timer); 
    cycles = (double)timer;
    return cycles; 
}
/* $end mountainfuns */


