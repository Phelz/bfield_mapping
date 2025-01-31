#pragma once

using namespace std;

#include <vector>
#include <cstdint>

#include "constants.h"


void biot_savart_chukman(const PRECISION_TYPE* s0,
    const PRECISION_TYPE* s1,
    const PRECISION_TYPE* r,
    PRECISION_TYPE* B);

void biot_savart_se(const PRECISION_TYPE* start,
    const PRECISION_TYPE* end,
    const PRECISION_TYPE* grid,
    PRECISION_TYPE* B);


// B-field Parallel for Chunkman
void calc_bfield_parallel_chukman(const uint32_t num_seg,                 // number of wire segments
    const PRECISION_TYPE* seg_start,      // start of segments
    const PRECISION_TYPE* seg_end,        // end of segments
    const uint32_t num_grid,              // number of grid points
    const PRECISION_TYPE* grid_xyz,      // the grid
    PRECISION_TYPE* bfield_xyz);         // Output B-field grid

// B-field Parallel for SE
void calc_bfield_parallel_se(const uint32_t num_seg,                     // number of wire segments
    const PRECISION_TYPE* seg_start,          // start of segments
    const PRECISION_TYPE* seg_end,            // end of segments
    const uint32_t num_grid,                  // number of grid points
    const PRECISION_TYPE* grid_xyz,          // the grid
    PRECISION_TYPE* bfield_xyz);