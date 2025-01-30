#pragma once

#include <vector>
#include <cmath>
#include <cstdint>

using namespace std;

const long double PI = 3.141592653589793238462643383279502884197169399375105820974944L;
const long double MU0 = 1.2566370612720e-6L;
const long double TOO_SMALL = 1e-12L;

vector<long double> cross_product(const vector<long double>& a, const vector<long double>& b);
long double dot_product(const vector<long double>& a, const vector<long double>& b);
vector<long double> subtract(const vector<long double>& a, const vector<long double>& b);
vector<long double> add(const vector<long double>& a, const vector<long double>& b);
vector<long double> scale(const vector<long double>& v, long double scalar);
long double norm(const vector<long double>& v);



