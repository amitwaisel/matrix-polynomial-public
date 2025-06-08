#pragma once

#define CALCULATE_REAL 1
#define CALCULATE_COMPLEX 0
#define USE_LAPACKE_FLAG 0

#define USE_CILK 1

#define MATRIX_ORDER LAPACK_COL_MAJOR

#define USE_LAPACKE  (USE_LAPACKE_FLAG || MATRIX_ORDER == LAPACK_ROW_MAJOR)

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
#define CBLAS_MATRIX_ORDER CblasRowMajor
#else
#define CBLAS_MATRIX_ORDER CblasColMajor
#endif

#if USE_CILK
#include <cilk/cilk.h>
#include <cilk/reducer_list.h>
#include <cilk/reducer_opadd.h>
#include <cilk/cilk_api.h>
#else
#define cilk_for for
#define cilk_spawn
#define cilk_sync do{} while(false)
#endif

#if CALCULATE_REAL && CALCULATE_COMPLEX
#error "Can't calculate both real and complex flavors, choose only one"
#endif
#if !(CALCULATE_REAL || CALCULATE_COMPLEX)
#error "You must choose either real or complex calculation flavors"
#endif