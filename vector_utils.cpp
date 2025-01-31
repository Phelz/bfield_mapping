#include "constants.h"
#include "vector_utils.h"

// Subtract two vectors: result = a - b
inline void subtract(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

// Compute cross product: result = a x b
inline void cross_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

// Compute dot product: returns a . b
inline PRECISION_TYPE dot_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute norm (magnitude) of a vector
inline PRECISION_TYPE norm(const PRECISION_TYPE* a) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// Scale a vector by a scalar: result = scalar * a
inline void scale(const PRECISION_TYPE* a, PRECISION_TYPE scalar, PRECISION_TYPE* result) {
    result[0] = scalar * a[0];
    result[1] = scalar * a[1];
    result[2] = scalar * a[2];
}

// Add two vectors: a += b
inline void add(PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}


//PRECISION_TYPE norm(const vector<PRECISION_TYPE>& v) {
//    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
//}
//
//void cross_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
//    result[0] = a[1] * b[2] - a[2] * b[1];
//    result[1] = a[2] * b[0] - a[0] * b[2];
//    result[2] = a[0] * b[1] - a[1] * b[0];
//}
//
//PRECISION_TYPE dot_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b) {
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//void subtract(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
//    result[0] = a[0] - b[0];
//    result[1] = a[1] - b[1];
//    result[2] = a[2] - b[2];
//}
//
//void scale(const PRECISION_TYPE* a, PRECISION_TYPE factor, PRECISION_TYPE* result) {
//    result[0] = a[0] * factor;
//    result[1] = a[1] * factor;
//    result[2] = a[2] * factor;
//}
//
//void add(PRECISION_TYPE* a, const PRECISION_TYPE* b) {
//    a[0] += b[0];
//    a[1] += b[1];
//    a[2] += b[2];
//}
//
//PRECISION_TYPE norm(const PRECISION_TYPE* a) {
//    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
//}



//vector<PRECISION_TYPE> cross_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
//	return {
//		a[1] * b[2] - a[2] * b[1],
//		a[2] * b[0] - a[0] * b[2],
//		a[0] * b[1] - a[1] * b[0],
//	};
//}

//PRECISION_TYPE dot_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
//	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//vector<PRECISION_TYPE> subtract(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
//	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
//}
//vector<PRECISION_TYPE> add(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
//	return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
//}
//vector<PRECISION_TYPE> scale(const vector<PRECISION_TYPE>& v, PRECISION_TYPE scalar) {
//	return { v[0] * scalar, v[1] * scalar, v[2] * scalar };
//}
