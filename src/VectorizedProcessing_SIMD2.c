/* Vectorized Processing Simulation */
#include <stdlib.h>
#include <stdio.h>
#include "Rdtsc.h"              /* startTimer, stopTimer */
#include <sys/time.h>           /* gettimeofday */
#include "x86intrin.h"              /* startTimer, stopTimer */

#define MAXELEMS 600000000 /* column length,base line 600000000 */
#define VECTORLEN (1<<20) /* maximal vector length */
#define MAXVECSTEP 20
#define TAX 0.2 /* define initiate column value, constang value for simplicity and verification */
#define EPRICE 100
#define QUANTITY 10
//typedef int data_t;     /* define column data type */
//typedef long data_t;
typedef float data_t;
//typedef double data_t;

void init_cols(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength);
double test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len);
double test_O(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len);
double test_O4(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len);
double SIMD_test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len);
double SIMD2_test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len);

/* $begin mountainmain */
int main()
{
    int columnlength=MAXELEMS; /*column length*/
    data_t *l_tax=(data_t *)malloc(sizeof(data_t) * MAXELEMS);      /*malloc array for three columns*/
    data_t *l_extendedprice=(data_t *)malloc(sizeof(data_t) * MAXELEMS);
    data_t *l_quantity=(data_t *)malloc(sizeof(data_t) * MAXELEMS);
    
    printf("initiate columns of l_tax,l_extendedprice,l_quantity...\n");
    init_cols(l_tax,l_extendedprice,l_quantity,MAXELEMS);
    
    double finalresult=0;
    finalresult=test(l_tax,l_extendedprice,l_quantity,MAXELEMS,1);  /*warm up*/
    
    double cycles;
    uint64_t timer;
    int i,vec_len=1;
    /**/
    printf("column-wise processing performance(cycle-per-tuple):\n");
    startTimer(&timer);
    finalresult=test(l_tax,l_extendedprice,l_quantity,MAXELEMS,MAXELEMS);  //warm up//
    stopTimer(&timer);
    cycles = timer/MAXELEMS;
    printf("vector length:%d result:%lf cycles:%lf\n", MAXELEMS,finalresult,cycles);
    
    printf("naive vectorized processing performance(cycle-per-tuple):\n");
    for(i=0;i<=MAXVECSTEP;i++){
        startTimer(&timer);
        finalresult=test(l_tax,l_extendedprice,l_quantity,MAXELEMS,vec_len);
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", vec_len,finalresult,cycles);
        vec_len=vec_len*2;
    }
    printf("optimized vectorized processing performance(cycle-per-tuple):\n");
    vec_len=1;
    for(i=0;i<=MAXVECSTEP;i++){
        startTimer(&timer);
        finalresult=test_O(l_tax,l_extendedprice,l_quantity,MAXELEMS,vec_len);
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", vec_len,finalresult,cycles);
        vec_len=vec_len*2;
    }

    printf("optimized*4 vectorized processing performance(cycle-per-tuple):\n");
    vec_len=1;
    for(i=0;i<=MAXVECSTEP;i++){
        startTimer(&timer);
        finalresult=test_O4(l_tax,l_extendedprice,l_quantity,MAXELEMS,vec_len);
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", vec_len,finalresult,cycles);
        vec_len=vec_len*2;
    }
    
    printf("SIMD vectorized processing performance(cycle-per-tuple):\n");
    vec_len=1;finalresult=0;
    for(i=0;i<=MAXVECSTEP;i++){
        startTimer(&timer);
        finalresult=SIMD_test(l_tax,l_extendedprice,l_quantity,MAXELEMS,vec_len);  /*warm up*/
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", vec_len,finalresult,cycles);
        vec_len=vec_len*2;
    }

    printf("SIMD processing performance(cycle-per-tuple):\n");
    vec_len=1;finalresult=0;
        startTimer(&timer);
        finalresult=SIMD_test(l_tax,l_extendedprice,l_quantity,MAXELEMS,MAXELEMS);  /*warm up*/
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", MAXELEMS,finalresult,cycles);

    printf("SIMD2 vectorized processing performance(cycle-per-tuple):\n");
    vec_len=1;finalresult=0;
    for(i=0;i<=MAXVECSTEP;i++){
        startTimer(&timer);
        finalresult=SIMD2_test(l_tax,l_extendedprice,l_quantity,MAXELEMS,vec_len);  /*warm up*/
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", vec_len,finalresult,cycles);
        vec_len=vec_len*2;
    }

    printf("SIMD2 processing performance(cycle-per-tuple):\n");
    vec_len=1;finalresult=0;
        startTimer(&timer);
        finalresult=SIMD2_test(l_tax,l_extendedprice,l_quantity,MAXELEMS,MAXELEMS);  /*warm up*/
        stopTimer(&timer);
        cycles = (double)timer/MAXELEMS;
        printf("vector length:%d result:%lf cycles:%4.2f\n", MAXELEMS,finalresult,cycles);
    exit(0);
}
/* init_data - initializes the array */
void init_cols(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength)
{
    int i;
    for (i = 0; i < columnlength; i++){
        l_tax[i] = TAX;
        l_extendedprice[i] = EPRICE;
        l_quantity[i] = QUANTITY;
    }
}
double test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len)
{
    data_t *netto_value=(data_t *)malloc(sizeof(data_t) * vector_len);  /*malloc 3 intermediator vector*/
    data_t *tax_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    data_t *total_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    //    memset(netto_value, 0, vector_len);     /*initiate array with 0*/
    //    memset(tax_value, 0, vector_len);
    //    memset(total_value, 0, vector_len);
    int i,j,k;
    double result = 0;
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*/
            netto_value[j]=l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len];
        }
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*(1+l_tax)*/
            tax_value[j]=1+l_tax[j+i*vector_len];
        }
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*l_tax*/
            total_value[j]=netto_value[j]*tax_value[j];
        }
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*l_tax*/
            result+=total_value[j];
        }
    }
    k=vector_len*i;
    for (; k<columnlength; k++){
        result+=l_extendedprice[k]*l_quantity[k]*(1+l_tax[k]);
    }
    return result; /* return final result */
}
double test_O(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len)
{
    data_t *netto_value=(data_t *)malloc(sizeof(data_t) * vector_len);  /*malloc 3 intermediator vector*/
    data_t *tax_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    data_t *total_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    //    memset(netto_value, 0, vector_len);     /*initiate array with 0*/
    //    memset(tax_value, 0, vector_len);
    //    memset(total_value, 0, vector_len);
    int i,j,k;
    double result = 0;
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*/
            netto_value[j]=l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len];
            tax_value[j]=1+l_tax[j+i*vector_len];
            total_value[j]=netto_value[j]*tax_value[j];
            result+=total_value[j];
        }
    }
    k=vector_len*i;
    for (; k<columnlength; k++){
        result+=l_extendedprice[k]*l_quantity[k]*(1+l_tax[k]);
    }
    return result; /* return final result */
}
double test_O4(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len)
{
    data_t *netto_value=(data_t *)malloc(sizeof(data_t) * vector_len);  /*malloc 3 intermediator vector*/
    data_t *tax_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    data_t *total_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    //    memset(netto_value, 0, vector_len);     /*initiate array with 0*/
    //    memset(tax_value, 0, vector_len);
    //    memset(total_value, 0, vector_len);
    int i,j,k;
    double result= 0;
	double res0=0,res1=0,res2=0,res3=0;
	if(vector_len<4){
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
        for(j=0; j<vector_len;j++){  /*computing with vectors for l_extendedprice*l_quantity*/
            netto_value[j]=l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len];
            tax_value[j]=1+l_tax[j+i*vector_len];
            total_value[j]=netto_value[j]*tax_value[j];
            result+=total_value[j];
        }  //end of loop j
    }  //end of loop i
	} //end of else
 	else{
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
        for(j=0; j<vector_len;j=j+4){  /*computing with vectors for l_extendedprice*l_quantity*/
            netto_value[j]=l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len];
            netto_value[j+1]=l_extendedprice[j+1+i*vector_len]*l_quantity[j+1+i*vector_len];
            netto_value[j+2]=l_extendedprice[j+2+i*vector_len]*l_quantity[j+2+i*vector_len];
            netto_value[j+3]=l_extendedprice[j+3+i*vector_len]*l_quantity[j+3+i*vector_len];
            tax_value[j]=1+l_tax[j+i*vector_len];
            tax_value[j+1]=1+l_tax[j+1+i*vector_len];
            tax_value[j+2]=1+l_tax[j+2+i*vector_len];
            tax_value[j+3]=1+l_tax[j+3+i*vector_len];
            total_value[j]=netto_value[j]*tax_value[j];
            total_value[j+1]=netto_value[j+1]*tax_value[j+1];
            total_value[j+2]=netto_value[j+2]*tax_value[j+2];
            total_value[j+3]=netto_value[j+3]*tax_value[j+3];
           res0=total_value[j];
           res1=total_value[j+1];
           res2=total_value[j+2];
           res3=total_value[j+3];
		   result+=res0+res1+res2+res3;
        }  //end of loop j
    }  //end of loop i
	} //end of else
    k=vector_len*i;
    for (; k<columnlength; k++){
        result+=l_extendedprice[k]*l_quantity[k]*(1+l_tax[k]);
    }
    return result; /* return final result */
}
double SIMD_test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len)
{
    
    int i,j,k,m,n;
    double result = 0;
    data_t *netto_value=(data_t *)malloc(sizeof(data_t) * vector_len);  /*malloc 3 intermediator vector*/
    data_t *tax_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    data_t *total_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    size_t nBlockWidth = 32.0/(sizeof(float));    //块宽. AVX寄存器能一次处理data_t类型数据的数量.
    size_t cntBlock = vector_len / nBlockWidth;    //向量中SIMD块数.
    size_t cntRem = vector_len % nBlockWidth;    //剩余数量.
    __m256 yfsLoad1;    //加载extendedprice向量.
    __m256 yfsLoad2;    //加载quantity向量.
    __m256 yfsLoad3;    //加载tax向量.
    __m256 yfsLoad4;    //加载1+tax向量.
    const data_t* p1 = l_extendedprice;    // AVX批量处理时所用的指针.
    const data_t* p2 = l_quantity;    // AVX批量处理时所用的指针.
    const data_t* p3 = l_tax;    // AVX批量处理时所用的指针.
    data_t* p4 = netto_value;    // AVX批量处理时所用的指针.
    data_t* p5 = tax_value;    // AVX批量处理时所用的指针.
    data_t* p6 = total_value;    // AVX批量处理时所用的指针.
    
    __m256 yfsSum1 = _mm256_setzero_ps();    //l_extendedprice*l_quantity临时结果变量。[AVX] 赋初值0
    __m256 yfsSum2 = _mm256_setzero_ps();    //1-l_tax临时结果变量。[AVX] 赋初值0
    __m256 yfsSum_plus_1 = _mm256_set1_ps(1);  // [AVX] 加载向量1
    __m256 yfsSum3 = _mm256_setzero_ps();    //聚集计算结果变量。[AVX] 赋初值0
    data_t* q = (data_t*)malloc(nBlockWidth * sizeof(data_t));    //将AVX变量上的多个数值合并时所用指针.
//    const data_t* q;
    
//    int cnt = 0;
//    double tmp = 0;
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
         /*computing with vectors for l_extendedprice*l_quantity*/
            if(vector_len<nBlockWidth) {
                for(j=0; j<vector_len;j++)
                {
                    result += l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len]*(1+l_tax[j+i*vector_len]);
                }
            } else{
                yfsSum3 = _mm256_setzero_ps();
                for(m=0; m<cntBlock; ++m){
                    //first, load vectors from columns//
                    yfsLoad1 = _mm256_loadu_ps(p1);    // [AVX] 加载l_extendedprice
                    yfsLoad2 = _mm256_loadu_ps(p2);    // [AVX] 加载l_quantity
                    yfsSum1 = _mm256_mul_ps(yfsLoad1, yfsLoad2);// [AVX] 单精浮点紧缩加法
                    //            _mm256_store_ps(p4,yfsSum1);
                    yfsLoad3 = _mm256_loadu_ps(p3);   // [AVX] 加载l_tax
                    yfsLoad4 = _mm256_add_ps(yfsSum_plus_1, yfsLoad3);// [AVX] 单精浮点紧缩加法执行1+l_tax
                    //            _mm256_store_ps(p5,yfsLoad4);
                    //             yfsSum2 = _mm256_mul_ps(_mm256_load_ps(p4), _mm256_load_ps(p5));// [AVX] 单精浮点紧缩加法
                    yfsSum2 = _mm256_mul_ps(yfsSum1, yfsLoad4);// [AVX] 单精浮点紧缩加法
                    //            _mm256_store_ps(p6,yfsSum2);
                    //            yfsSum3 = _mm256_add_ps(yfsSum3, _mm256_load_ps(p6));
                    yfsSum3 = _mm256_add_ps(yfsSum3, yfsSum2);
                    p1 += nBlockWidth;
                    p2 += nBlockWidth;
                    p3 += nBlockWidth;
                    p4 += nBlockWidth;
                    p5 += nBlockWidth;
                    p6 += nBlockWidth;
                }
                //aggregate vector results//
//                q = (const data_t*)&yfsSum3;
                _mm256_storeu_ps((float*)q, yfsSum3);
                result += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
                for(n=0; n<cntRem; ++n){
                    result +=(*p1)*(*p2)*(1+(*p3));
                    p1++;p2++;p3++;
                }
            }
//        if (i == 0) printf("result = %f\n", result);
//        if (i == columnlength/vector_len - 1) printf("tmp = %f\n", tmp);
//        if (i <= 10) printf("tmp = %f\n", tmp);
//        ++cnt;
    }
//    printf("cnt = %d\n", cnt);
        k=vector_len*i;
     for (; k<columnlength; k++){
     result+=l_extendedprice[k]*l_quantity[k]*(1+l_tax[k]);
     }
         return result; /* return final result */
}

double SIMD2_test(data_t *l_tax, data_t *l_extendedprice, data_t *l_quantity,int columnlength,int vector_len)
{
    
    int i,j,k,m,n;
    double result = 0;
    data_t *netto_value=(data_t *)malloc(sizeof(data_t) * vector_len);  /*malloc 3 intermediator vector*/
    data_t *tax_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    data_t *total_value=(data_t *)malloc(sizeof(data_t) * vector_len);
    size_t nBlockWidth = 32.0/(sizeof(float));    //块宽. AVX寄存器能一次处理data_t类型数据的数量.
    //size_t nBlockWidth = ;    //块宽. AVX寄存器能一次处理data_t类型数据的数量.
    size_t cntBlock = vector_len / (nBlockWidth * 4);    //向量中SIMD块数.
    size_t cntRem = vector_len % (nBlockWidth * 4);    //剩余数量.
    __m256 yfsLoad1;    //加载extendedprice向量.
    __m256 yfsLoad2;    //加载quantity向量.
    __m256 yfsLoad3;    //加载tax向量.
    __m256 yfsLoad4;    //加载1+tax向量.
    const data_t* p1 = l_extendedprice;    // AVX批量处理时所用的指针.
    const data_t* p2 = l_quantity;    // AVX批量处理时所用的指针.
    const data_t* p3 = l_tax;    // AVX批量处理时所用的指针.
    data_t* p4 = netto_value;    // AVX批量处理时所用的指针.
    data_t* p5 = tax_value;    // AVX批量处理时所用的指针.
    data_t* p6 = total_value;    // AVX批量处理时所用的指针.
    
    __m256 yfsSum1 = _mm256_setzero_ps();    //l_extendedprice*l_quantity临时结果变量。[AVX] 赋初值0
    __m256 yfsSum2 = _mm256_setzero_ps();    //1-l_tax临时结果变量。[AVX] 赋初值0
    __m256 yfsSum_plus_1 = _mm256_set1_ps(1);  // [AVX] 加载向量1
    __m256 yfsSum3 = _mm256_setzero_ps();    //聚集计算结果变量。[AVX] 赋初值0
    __m256 yfsSum3_1 = _mm256_setzero_ps();    //聚集计算结果变量。[AVX] 赋初值0
    __m256 yfsSum3_2 = _mm256_setzero_ps();    //聚集计算结果变量。[AVX] 赋初值0
    __m256 yfsSum3_3 = _mm256_setzero_ps();    //聚集计算结果变量。[AVX] 赋初值0
    data_t* q = (data_t*)malloc(nBlockWidth * sizeof(data_t));    //将AVX变量上的多个数值合并时所用指针.
//    const data_t* q;
    
//    int cnt = 0;
//    double tmp = 0;
    for (i = 0; i < columnlength/vector_len; i++) {  /*loop with vector granularity*/
         /*computing with vectors for l_extendedprice*l_quantity*/
            if(vector_len<nBlockWidth) {
                for(j=0; j<vector_len;j++)
                {
                    result += l_extendedprice[j+i*vector_len]*l_quantity[j+i*vector_len]*(1+l_tax[j+i*vector_len]);
                }
            } else{
                yfsSum3 = _mm256_setzero_ps();
                yfsSum3_1 = _mm256_setzero_ps();
                yfsSum3_2 = _mm256_setzero_ps();
                yfsSum3_3 = _mm256_setzero_ps();
                for(m=0; m<cntBlock; ++m){
                    //first, load vectors from columns//
                    yfsLoad1 = _mm256_loadu_ps(p1);    // [AVX] 加载l_extendedprice
                    yfsLoad2 = _mm256_loadu_ps(p2);    // [AVX] 加载l_quantity
                    yfsSum1 = _mm256_mul_ps(yfsLoad1, yfsLoad2);// [AVX] 单精浮点紧缩加法
                    yfsLoad3 = _mm256_loadu_ps(p3);   // [AVX] 加载l_tax
                    yfsLoad4 = _mm256_add_ps(yfsSum_plus_1, yfsLoad3);// [AVX] 单精浮点紧缩加法执行1+l_tax
                    yfsSum2 = _mm256_mul_ps(yfsSum1, yfsLoad4);// [AVX] 单精浮点紧缩加法
                    yfsSum3 = _mm256_add_ps(yfsSum3, yfsSum2);
		    
                    yfsLoad1 = _mm256_loadu_ps(p1+nBlockWidth);    // [AVX] 加载l_extendedprice
                    yfsLoad2 = _mm256_loadu_ps(p2+nBlockWidth);    // [AVX] 加载l_quantity
                    yfsSum1 = _mm256_mul_ps(yfsLoad1, yfsLoad2);// [AVX] 单精浮点紧缩加法
                    yfsLoad3 = _mm256_loadu_ps(p3+nBlockWidth);   // [AVX] 加载l_tax
                    yfsLoad4 = _mm256_add_ps(yfsSum_plus_1, yfsLoad3);// [AVX] 单精浮点紧缩加法执行1+l_tax
                    yfsSum2 = _mm256_mul_ps(yfsSum1, yfsLoad4);// [AVX] 单精浮点紧缩加法
                    yfsSum3_1 = _mm256_add_ps(yfsSum3_1, yfsSum2);

                    yfsLoad1 = _mm256_loadu_ps(p1+nBlockWidth*2);    // [AVX] 加载l_extendedprice
                    yfsLoad2 = _mm256_loadu_ps(p2+nBlockWidth*2);    // [AVX] 加载l_quantity
                    yfsSum1 = _mm256_mul_ps(yfsLoad1, yfsLoad2);// [AVX] 单精浮点紧缩加法
                    yfsLoad3 = _mm256_loadu_ps(p3+nBlockWidth*2);   // [AVX] 加载l_tax
                    yfsLoad4 = _mm256_add_ps(yfsSum_plus_1, yfsLoad3);// [AVX] 单精浮点紧缩加法执行1+l_tax
                    yfsSum2 = _mm256_mul_ps(yfsSum1, yfsLoad4);// [AVX] 单精浮点紧缩加法
                    yfsSum3_2 = _mm256_add_ps(yfsSum3_2, yfsSum2);

                    yfsLoad1 = _mm256_loadu_ps(p1+nBlockWidth*3);    // [AVX] 加载l_extendedprice
                    yfsLoad2 = _mm256_loadu_ps(p2+nBlockWidth*3);    // [AVX] 加载l_quantity
                    yfsSum1 = _mm256_mul_ps(yfsLoad1, yfsLoad2);// [AVX] 单精浮点紧缩加法
                    yfsLoad3 = _mm256_loadu_ps(p3+nBlockWidth*3);   // [AVX] 加载l_tax
                    yfsLoad4 = _mm256_add_ps(yfsSum_plus_1, yfsLoad3);// [AVX] 单精浮点紧缩加法执行1+l_tax
                    yfsSum2 = _mm256_mul_ps(yfsSum1, yfsLoad4);// [AVX] 单精浮点紧缩加法
                    yfsSum3_3 = _mm256_add_ps(yfsSum3_3, yfsSum2);

                    p1 += nBlockWidth * 4;
                    p2 += nBlockWidth * 4;
                    p3 += nBlockWidth * 4;
                }
                //aggregate vector results//
                //q = (const data_t*)&yfsSum3;
		//yfsSum_plus_1 = _mm256_add_ps(yfsSum_plus_1, yfsSum1);
		//yfsSum2 = _mm256_add_ps(yfsSum2, yfsSum3);
		//yfsSum_plus_1 = _mm256_add_ps(yfsSum_plus_1, yfsSum2);
		yfsSum3 = _mm256_add_ps(yfsSum3,yfsSum3_1);
		yfsSum3_2 = _mm256_add_ps(yfsSum3_2,yfsSum3_3);
		yfsSum3 = _mm256_add_ps(yfsSum3,yfsSum3_2);
                _mm256_storeu_ps((float*)q, yfsSum3);
		//q = (const float*)&yfsSum_plus_1;
                result += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
                for(n=0; n<cntRem; ++n){
                    result +=(*p1)*(*p2)*(1+(*p3));
                    p1++;p2++;p3++;
                }
            }
//        if (i == 0) printf("result = %f\n", result);
//        if (i == columnlength/vector_len - 1) printf("tmp = %f\n", tmp);
//        if (i <= 10) printf("tmp = %f\n", tmp);
//        ++cnt;
    }
//    printf("cnt = %d\n", cnt);
        k=vector_len*i;
     for (; k<columnlength; k++){
     result+=l_extendedprice[k]*l_quantity[k]*(1+l_tax[k]);
     }
         return result; /* return final result */
}
