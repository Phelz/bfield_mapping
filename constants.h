#pragma once

const long double PI = 3.141592653589793238462643383279502884197169399375105820974944L;
const long double MU0 = 1.2566370612720e-6L;
const long double TOO_SMALL = 1e-12L;

//#define PRECISION_TYPE double
#define PRECISION_TYPE long double

constexpr auto OUTPUT_PRECISION = 20;
constexpr auto NUM_THREADS = 4;
constexpr auto TRACES_DIR = "D:\\OneDrive - University of Calgary\\PARA\\Projects\\Inventor\\Magnetic Quadruple\\comsol\\test_2\\particle_data";