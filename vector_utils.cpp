#include "vector_utils.h"

vector<long double> cross_product(const vector<long double>& a, const vector<long double>& b) {
	return {
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0],
	};
}

long double dot_product(const vector<long double>& a, const vector<long double>& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

vector<long double> subtract(const vector<long double>& a, const vector<long double>& b) {
	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}
vector<long double> add(const vector<long double>& a, const vector<long double>& b) {
	return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}
vector<long double> scale(const vector<long double>& v, long double scalar) {
	return { v[0] * scalar, v[1] * scalar, v[2] * scalar };
}

long double norm(const vector<long double>& v) {
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
