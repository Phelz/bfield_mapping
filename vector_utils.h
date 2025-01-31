#pragma once

using namespace std;

#include <vector>
#include <cmath>
#include <cstdint>

#include "constants.h"

vector<PRECISION_TYPE> cross_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
PRECISION_TYPE dot_product(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
vector<PRECISION_TYPE> subtract(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
vector<PRECISION_TYPE> add(const vector<PRECISION_TYPE>& a, const vector<PRECISION_TYPE>& b);
vector<PRECISION_TYPE> scale(const vector<PRECISION_TYPE>& v, PRECISION_TYPE scalar);
PRECISION_TYPE norm(const vector<PRECISION_TYPE>& v);



