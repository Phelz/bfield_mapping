#pragma once

#include "constants.h"
#include <cmath>  // for sqrt()

inline void subtract(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result);
inline void cross_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result);
inline PRECISION_TYPE dot_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b);
inline PRECISION_TYPE norm(const PRECISION_TYPE* a);
inline void scale(const PRECISION_TYPE* a, PRECISION_TYPE scalar, PRECISION_TYPE* result);
inline void add(PRECISION_TYPE* a, const PRECISION_TYPE* b);



//vector<PRECISION_TYPE> cross_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
//PRECISION_TYPE dot_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
//vector<PRECISION_TYPE> subtract(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
//vector<PRECISION_TYPE> add(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
//vector<PRECISION_TYPE> scale(const vector<PRECISION_TYPE>& v, PRECISION_TYPE scalar);
//PRECISION_TYPE norm(const vector<PRECISION_TYPE>& v);



