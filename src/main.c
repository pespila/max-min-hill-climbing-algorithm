#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP MyMatrix(SEXP RinMatrix) {
    SEXP Rdim = getAttrib(RinMatrix,R_DimSymbol);
    int I = INTEGER(Rdim)[0];
    int J = INTEGER(Rdim)[1];

    RinMatrix = coerceVector(RinMatrix,REALSXP);

    SEXP Rval;
    PROTECT(Rval = allocMatrix(REALSXP,I,J));

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++)
        {
            REAL(Rval)[i+I*j] = REAL(RinMatrix)[i+I*j];
            // Rprintf("%f ",REAL(Rval)[i+I*j]);
            // if(i<=j) REAL(Rval)[i+I*j] = REAL(RinMatrix)[i+I*j];
            // else REAL(Rval)[i+I*j] = 0;
        }
        // Rprintf("\n");
    }

    UNPROTECT(1);
    return Rval;
}


// /* conditional mutual information, to be used for the asymptotic test. */
// SEXP cmi(SEXP x, SEXP y, SEXP z, SEXP gsquare) {

// int llx = NLEVELS(x), lly = NLEVELS(y), llz = NLEVELS(z);
// int num = length(x);
// int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
// double *res = NULL;
// SEXP result;

//    allocate and initialize result to zero. 
//   PROTECT(result = allocVector(REALSXP, 2));
//   res = REAL(result);
//   res[0] = c_cmi(xx, &llx, yy, &lly, zz, &llz, &num);
//   res[1] = (double)(llx - 1) * (double)(lly - 1) * (double)llz;

//   /* rescale to match the G^2 test. */
//   if (isTRUE(gsquare))
//     res[0] *= 2 * num;

//   UNPROTECT(1);

//   return result;

// }/*CMI*/




// #include "classes.h"

// int main(int argc, char const *argv[]) {
//     // const char *method;
//     int arg;
//     if(argc>1) arg=atoi(&*argv[1]);
//     else if(argc>1 && atoi(&*argv[1])>1000) arg=1000;
//     else arg=5;

//     // ScoreMatrix A(arg);
//     // A.Print();

//     int* inner;
//     inner=(int*)malloc(arg*sizeof(int));
//     int** matrix;
//     matrix=(int**)malloc(arg*sizeof(inner));

//     for (int i = 0; i < arg; i++)
//     {
//       matrix[i]=inner;
//     }
//     // inner=malloc(arg*sizeof(int));
//     // int** matrix=inner;
//     // matrix=malloc(arg*sizeof(inner));

//     for (int i = 0; i < arg; i++) {
//        	for (int j = 0; j < arg; j++)	{
//             matrix[i][j]=1;
//             printf("%d ", matrix[i][j]);
//         }
//         printf("\n");
//     }

//     // for (int i = 0; i < arg; i++)
//     // {
//     //     free(matrix[i][0]);
//     // }
//     // free(matrix);

//     free(inner);
//     free(matrix);

//     return EXIT_SUCCESS;


// 	   value = malloc(size*sizeof(int));
// 	   static jump_buf s_jumpBuffer;
// 	   if(setjmp(s_jumpBuffer)) {

// 	   }
//    if(NULL == value) {
//       printf("Fehler bei malloc...!!\n");
//       return EXIT_FAILURE;
//    }
//    while(i < size) {
//       printf("Wert für value[%d] eingeben : ",i);
//       scanf("%d",(value+i));
//       i++;
//    }
//    printf("Hier Ihre Werte\n");
//    for(i=0; i<size; i++)
//       printf("value[%d] = %d\n", i, *(value+i));
//    return EXIT_SUCCESS;
    
//     return EXIT_SUCCESS;
// }

// static jump_buf s_jumpBuffer;

// void Example() { 
//   if (setjmp(s_jumpBuffer)) {
//     // The longjmp was executed and returned control here
//     printf("Exception happened\n");
//   } else {
//     // Normal code execution starts here
//     Test();
//   }
// }