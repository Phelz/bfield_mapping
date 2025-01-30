#pragma once

#include <vector>
#include <cstdint>
#include <functional>

using namespace std;

void biot_savart_chukman(const vector<long double>& s0,
    const vector<long double>& s1,
    const vector<long double>& r,
    vector<long double>& B);

void biot_savart_se(const vector<long double>& start,
    const vector<long double>& end,
    const vector<long double>& grid,
    vector<long double>& b1);

void calc_bfield_parallel(const uint32_t num_seg,                     // number of wire segments
    const vector<vector<long double>>& seg_start, // start of segments
    const vector<vector<long double>>& seg_end,   // end of segments
    const uint32_t num_grid,                    // number of grid points
    const vector<vector<long double>>& grid_xyz, // the grid
    vector<vector<long double>>& b1_xyz,
    function<void(const vector<long double>&,
        const vector<long double>&,
        const vector<long double>&,
        vector<long double>&)> biot_savart_method);        // Output B-field grid

