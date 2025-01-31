#include "constants.h"
#include "vector_utils.h"

vector<PRECISION_TYPE> cross_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
	return {
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0],
	};
}

PRECISION_TYPE dot_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

vector<PRECISION_TYPE> subtract(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}
vector<PRECISION_TYPE> add(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b) {
	return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}
vector<PRECISION_TYPE> scale(const vector<PRECISION_TYPE>& v, PRECISION_TYPE scalar) {
	return { v[0] * scalar, v[1] * scalar, v[2] * scalar };
}

PRECISION_TYPE norm(const vector<PRECISION_TYPE>& v) {
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
