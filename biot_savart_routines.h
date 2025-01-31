#pragma once

using namespace std;

#include <vector>
#include <cstdint>
#include <functional>

#include "constants.h"

void biot_savart_chukman(	const vector<PRECISION_TYPE>& s0,
							const vector<PRECISION_TYPE>& s1,
							const vector<PRECISION_TYPE>& r,
							vector<PRECISION_TYPE>& B);

void biot_savart_se(const vector<PRECISION_TYPE>& start,
					const vector<PRECISION_TYPE>& end,
					const vector<PRECISION_TYPE>& grid,
					vector<PRECISION_TYPE>& b1);


// Bfield Parallel for Chunkman
void calc_bfield_parallel_chukman(const uint32_t num_seg,                     // number of wire segments
	const vector<vector<PRECISION_TYPE>>& seg_start, // start of segments
	const vector<vector<PRECISION_TYPE>>& seg_end,   // end of segments
	const uint32_t num_grid,                    // number of grid points
	const vector<vector<PRECISION_TYPE>>& grid_xyz, // the grid
	vector<vector<PRECISION_TYPE>>& b1_xyz);        // Output B-field grid

// Bfield Parallel for SE
void calc_bfield_parallel_se(const uint32_t num_seg,                     // number of wire segments
	const vector<vector<PRECISION_TYPE>>& seg_start, // start of segments
	const vector<vector<PRECISION_TYPE>>& seg_end,   // end of segments
	const uint32_t num_grid,                    // number of grid points
	const vector<vector<PRECISION_TYPE>>& grid_xyz, // the grid
	vector<vector<PRECISION_TYPE>>& b1_xyz);        // Output B-field grid
