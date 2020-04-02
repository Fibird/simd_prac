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
void test_SIMD_partition(void *args); 
int SPNLJ(int *pk, int pk_len, int *fk, int fk_len, int tnum); 

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

typedef struct {
    int *pk;
    int pk_len;
    int *fk;
    int fk_len;
    int result;
} Args;

int SPNLJ(int *pk, int pk_len, int *fk, int fk_len, int tnum) {
    int part_size = (pk_len) / tnum ;
    int remain_size = pk_len % tnum;
    int i = 0; 
    pthread_t tid[tnum];
    Args args[tnum];
    for (i = 0; i < tnum; i++) {
        args[i].pk = i * part_size + pk;
        args[i].pk_len = part_size;
        args[i].fk = fk;  args[i].fk_len = fk_len;
        args[i].result = 0;
        if (i == (tnum - 1) && remain_size != 0) args->pk_len = remain_size;
        pthread_create(&tid[i], NULL, test_SIMD_partition, (void*)&args[i]);
    }

    int result = 0;
    for (i = 0; i < tnum; i++) {
        pthread_join(tid[i], NULL);
        result += args[i].result; 
    }
    printf("[%d]", result);
}

void test_SIMD_partition(void *args) {
        int *pk = ((Args*)args)->pk;
        int pk_len = ((Args*)args)->pk_len;
        int *fk = ((Args*)args)->fk;
        int fk_len = ((Args*)args)->fk_len;
        ((Args*)args)->result = 0;

        printf("pk len:%d\n", pk_len);
        printf("fk len:%d\n", fk_len);
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
	for (k = 0; k < pk_len; k += partition) {
		int end = k + partition < pk_len ? k + partition : pk_len;
		int bound = (end - k) / 8 * 8 + k;
		for (i = 0; i < fk_len; ++i) {
			yidOuter = _mm256_set1_epi32(fk[i]);
			p = pk + k;
			for (j = k; j < bound; j+= 8) {
				yidInner = _mm256_loadu_si256(p);
				yidTmp = _mm256_cmpeq_epi32(yidOuter, yidInner);
				yidResult = _mm256_sub_epi32(yidResult, yidTmp);
				p += 32;
			}
			for (; j < end; ++j) {
				if (fk[i] == pk[j]) {
					++ret;
				}
			}
		}
	}
	_mm256_storeu_si256((void *)q, yidResult);
	ret += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
	sink = ret;
        ((Args*)args)->result = sink;
	//return sink;
}

long run(int *data, int *vector, int n, int vec_len)
{   
    long cycles;
    uint64_t timer;
    SPNLJ(vector,vec_len,data,n, 1);                     /* warm up the cache */       //line:mem:warmup
    startTimer(&timer);
    SPNLJ(vector,vec_len,data,n, 1);                     /* warm up the cache */       //line:mem:warmup
    stopTimer(&timer); 
    cycles = timer;  
    return cycles; 
}
/* $end mountainfuns */
