using namespace std;

#include <iostream>
#include <functional>
#include <omp.h>
#include <chrono>

#include "constants.h"
#include "vector_utils.h"  // Include the header where vector operations are defined FIRST
#include "biot_savart_routines.h"


void biot_savart_chukman(   const vector<PRECISION_TYPE>& s0,
                            const vector<PRECISION_TYPE>& s1,
                            const vector<PRECISION_TYPE>& r,
                            vector<PRECISION_TYPE>& B) {

    vector<PRECISION_TYPE> i, k0, k1, A;
    PRECISION_TYPE norm_i, norm_k0, norm_k1, norm_A_sqr, a0;

    i = subtract(s1, s0);
    norm_i = norm(i);
    if (norm_i < TOO_SMALL) return;
    norm_i = 1. / norm_i;
    i = scale(i, norm_i);

    k0 = subtract(s0, r);
    k1 = subtract(s1, r);
    norm_k0 = norm(k0);
    norm_k1 = norm(k1);
    if (norm_k0 < TOO_SMALL) return; norm_k0 = 1. / norm_k0;
    if (norm_k1 < TOO_SMALL) return; norm_k1 = 1. / norm_k1;

    A = cross_product(k0, i);
    norm_A_sqr = (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
    if (norm_A_sqr < TOO_SMALL) return;
    norm_A_sqr = 1. / norm_A_sqr;

    A = scale(A, norm_A_sqr);
    k0 = scale(k0, norm_k0);
    k1 = scale(k1, norm_k1);

    a0 = dot_product(subtract(k1, k0), i);
    //B = scale(A, a0);
    B = add(B, scale(A, a0));
}

void biot_savart_se(const vector<PRECISION_TYPE>& start,
                    const vector<PRECISION_TYPE>& end,
                    const vector<PRECISION_TYPE>& grid,
                    vector<PRECISION_TYPE>& b1) {
    // Implementation of the Biot-Savart SE method (as you provided)
    vector<PRECISION_TYPE> rf = subtract(end, grid); // r_(i+1)    
    vector<PRECISION_TYPE> ri = subtract(start, grid); // r_i
    vector<PRECISION_TYPE> dl = subtract(end, start); // dl

    PRECISION_TYPE norm_dl = norm(dl);
    PRECISION_TYPE norm_rf = norm(rf);
    PRECISION_TYPE norm_ri = norm(ri);

    vector<PRECISION_TYPE> dist = cross_product(rf, dl); // r cross dl
    dist = scale(dist, 1. / norm_dl);

    PRECISION_TYPE distance = norm(dist); // r

    PRECISION_TYPE cos_theta_i = dot_product(rf, dl) / (norm_rf * norm_dl);
    PRECISION_TYPE cos_theta_f = dot_product(ri, dl) / (norm_ri * norm_dl);

    PRECISION_TYPE absB1 = distance > TOO_SMALL ? (cos_theta_i - cos_theta_f) / distance : 0.;

    vector<PRECISION_TYPE> dir = cross_product(rf, dl);
    PRECISION_TYPE norm_dir = norm(dir);

    absB1 = norm_dir > TOO_SMALL ? absB1 / norm_dir : 0.;

    b1 = add(b1, scale(dir, absB1));
}


//void calc_bfield_parallel_se(   const uint32_t num_seg,
//                                const vector<vector<PRECISION_TYPE>>& seg_start,
//                                const vector<vector<PRECISION_TYPE>>& seg_end,
//                                const uint32_t num_grid,
//                                const vector<vector<PRECISION_TYPE>>& grid_xyz,
//                                vector<vector<PRECISION_TYPE>>& bfield_xyz) {
//    
//    // Reset magnetic field grid
//    for (auto& b : bfield_xyz) {
//        b = { 0.0, 0.0, 0.0 };
//    }
//
//    cout << "Starting B-field Calculation (SE)..." << endl;
//    auto start_time = chrono::high_resolution_clock::now();  // Start timer
//
//    //#pragma omp parallel for
//    for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
//        const auto& start_point = seg_start[seg_indx];
//        const auto& end_point = seg_end[seg_indx];
//
//        try {
//
//            // Parallel loop over grid points
//            #pragma omp parallel for
//            for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//                biot_savart_se(start_point, end_point, grid_xyz[grid_indx], bfield_xyz[grid_indx]);
//            }
//        }
//        catch (exception& ex) {
//            cout << "Simulation failed: " << ex.what() << endl;
//        }
//
//        if (seg_indx % 10 == 0) {
//            cout << "Progress (SE): " << seg_indx << " / " << num_seg << endl;
//        }
//    }
//
//    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
//    chrono::duration<double> elapsed = end_time - start_time;
//    cout << "Simulation finished (SE). Time elapsed: " << elapsed.count() << " seconds" << endl;
//}

void calc_bfield_parallel_se(const uint32_t num_seg,
                            const vector<vector<PRECISION_TYPE>>& seg_start,
                            const vector<vector<PRECISION_TYPE>>& seg_end,
                            const uint32_t num_grid,
                            const vector<vector<PRECISION_TYPE>>& grid_xyz,
                            vector<vector<PRECISION_TYPE>>& bfield_xyz) {

    // Reset magnetic field grid
    for (auto& b : bfield_xyz) {
        b = { 0.0, 0.0, 0.0 };
    }

    cout << "Starting B-field Calculation (SE)..." << endl;
    auto start_time = chrono::high_resolution_clock::now();  // Start timer

    //#pragma omp parallel for
    for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const auto& start_point = seg_start[seg_indx];
        const auto& end_point = seg_end[seg_indx];

        #pragma omp parallel for
        for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
            biot_savart_se(start_point, end_point, grid_xyz[grid_indx], bfield_xyz[grid_indx]);
        }

        if (seg_indx % 10 == 0) {
            cout << "Progress (SE): " << seg_indx << " / " << num_seg << endl;
        }
    }


    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Simulation finished (SE). Time elapsed: " << elapsed.count() << " seconds" << endl;
}




//void calc_bfield_parallel_se(const uint32_t num_seg,
//    const vector<vector<PRECISION_TYPE>>& seg_start,
//    const vector<vector<PRECISION_TYPE>>& seg_end,
//    const uint32_t num_grid,
//    const vector<vector<PRECISION_TYPE>>& grid_xyz,
//    vector<vector<PRECISION_TYPE>>& bfield_xyz) {
//
//    // Reset magnetic field grid
//    for (auto& b : bfield_xyz) {
//        b = { 0.0, 0.0, 0.0 };
//    }
//
//    cout << "Starting B-field Calculation (SE)..." << endl;
//    auto start_time = chrono::high_resolution_clock::now();  // Start timer
//
//    // Parallel processing for the segments
//#pragma omp parallel
//    {
//        // Create a local copy of bfield_xyz for each thread
//        vector<vector<PRECISION_TYPE>> local_bfield_xyz(num_grid, vector<PRECISION_TYPE>(3, 0.0));
//
//        // Loop over segments
//#pragma omp for
//        for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
//            const auto& start_point = seg_start[seg_indx];
//            const auto& end_point = seg_end[seg_indx];
//
//            // Loop over grid points
//            for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//                biot_savart_se(start_point, end_point, grid_xyz[grid_indx], local_bfield_xyz[grid_indx]);
//            }
//
//            if (seg_indx % 10 == 0) {
//                cout << "Progress (SE): " << seg_indx << " / " << num_seg << endl;
//            }
//        }
//
//        // Combine the local bfield_xyz results into the global bfield_xyz
//#pragma omp for
//        for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//            for (int i = 0; i < 3; ++i) {
//                // Sum the contributions from all threads
//#pragma omp atomic
//                bfield_xyz[grid_indx][i] += local_bfield_xyz[grid_indx][i];
//            }
//        }
//    }
//
//    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
//    chrono::duration<double> elapsed = end_time - start_time;
//    cout << "Simulation finished (SE). Time elapsed: " << elapsed.count() << " seconds" << endl;
//}



void calc_bfield_parallel_chukman(  const uint32_t num_seg,
                                    const vector<vector<PRECISION_TYPE>>& seg_start,
                                    const vector<vector<PRECISION_TYPE>>& seg_end,
                                    const uint32_t num_grid,
                                    const vector<vector<PRECISION_TYPE>>& grid_xyz,
                                    vector<vector<PRECISION_TYPE>>& bfield_xyz) {
    // Reset magnetic field grid
    for (auto& b : bfield_xyz) {
        b = { 0.0, 0.0, 0.0 };
    }

    cout << "Starting B-field Calculation (Chukman)..." << endl;
    auto start_time = chrono::high_resolution_clock::now();  // Start timer

    //#pragma omp parallel for
    for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const auto& start_point = seg_start[seg_indx];
        const auto& end_point = seg_end[seg_indx];

        try {
            // Parallel loop over grid points
            //#pragma omp parallel for
            for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++) {
                biot_savart_chukman(start_point, end_point, grid_xyz[grid_indx], bfield_xyz[grid_indx]);
            }
        }
        catch (exception& ex) {
            cout << "Simulation failed: " << ex.what() << endl;
        }

        if (seg_indx % 10 == 0) {
            cout << "Progress (Chukman): " << seg_indx << " / " << num_seg << endl;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Simulation finished (Chukman). Time elapsed: " << elapsed.count() << " seconds" << endl;
}
