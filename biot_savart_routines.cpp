using namespace std;

#include <iostream>
#include <functional>
#include <omp.h>
#include <chrono>

#include "constants.h"
//#include "vector_utils.h"  // Include the header where vector operations are defined FIRST
#include "biot_savart_routines.h"

inline void subtract(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

inline void cross_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

inline PRECISION_TYPE dot_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline PRECISION_TYPE norm(const PRECISION_TYPE* a) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline void scale(const PRECISION_TYPE* a, PRECISION_TYPE scalar, PRECISION_TYPE* result) {
    result[0] = scalar * a[0];
    result[1] = scalar * a[1];
    result[2] = scalar * a[2];
}

inline void add(PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}


void biot_savart_chukman(const PRECISION_TYPE* s0,
    const PRECISION_TYPE* s1,
    const PRECISION_TYPE* r,
    PRECISION_TYPE* B) {

    PRECISION_TYPE i[3], k0[3], k1[3], A[3];
    PRECISION_TYPE norm_i, norm_k0, norm_k1, norm_A_sqr, a0;

    // i = s1 - s0
    subtract(s1, s0, i);
    norm_i = norm(i);
    if (norm_i < TOO_SMALL) return;
    scale(i, 1.0 / norm_i, i); // Normalize i

    subtract(s0, r, k0); // k0 = s0 - r
    subtract(s1, r, k1); // k1 = s1 - r

    norm_k0 = norm(k0);
    norm_k1 = norm(k1);
    if (norm_k0 < TOO_SMALL || norm_k1 < TOO_SMALL) return;

    // A = cross(k0, i)
    cross_product(k0, i, A);
    norm_A_sqr = dot_product(A, A);
    if (norm_A_sqr < TOO_SMALL) return;

    scale(A, 1.0 / norm_A_sqr, A); 
    scale(k0, 1.0 / norm_k0, k0); 
    scale(k1, 1.0 / norm_k1, k1); 


    // a0 = dot(k1 - k0, i)
    PRECISION_TYPE k1_minus_k0[3];
    subtract(k1, k0, k1_minus_k0);
    a0 = dot_product(k1_minus_k0, i);

    // Accumulate the result: B = B + A * a0
    PRECISION_TYPE scaled_A[3];
    scale(A, a0, scaled_A); // Scale A by a0
    add(B, scaled_A);       // Add scaled A to B (accumulation)

}


void calc_bfield_parallel_chukman(  const uint32_t num_seg,
                                    const PRECISION_TYPE* seg_start,
                                    const PRECISION_TYPE* seg_end,
                                    const uint32_t num_grid,
                                    const PRECISION_TYPE* grid_xyz,
                                    PRECISION_TYPE* bfield_xyz) {

    // Reset magnetic field grid
    fill(bfield_xyz, bfield_xyz + 3 * num_grid, 0.0);

    cout << "Starting B-field Calculation (Chukman)..." << endl;
    auto start_time = chrono::high_resolution_clock::now();

    for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const PRECISION_TYPE* start_point = &seg_start[3 * seg_indx];
        const PRECISION_TYPE* end_point = &seg_end[3 * seg_indx];

        for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
            biot_savart_chukman(start_point, end_point, &grid_xyz[3 * grid_indx], &bfield_xyz[3 * grid_indx]);
        }

        if (seg_indx % 10 == 0) {
            cout << "Progress (Chukman): " << seg_indx << " / " << num_seg << endl;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Simulation finished (Chukman). Time elapsed: " << elapsed.count() << " seconds" << endl;
}

void biot_savart_se(const PRECISION_TYPE* start, const PRECISION_TYPE* end,
                    const PRECISION_TYPE* grid, PRECISION_TYPE* B) {
    
    // Compute vectors rf = end - grid, ri = start - grid, and dl = end - start
    PRECISION_TYPE rf[3], ri[3], dl[3];

    subtract(end, grid, rf);
    subtract(start, grid, ri);
    subtract(end, start, dl);

    // Compute norms
    PRECISION_TYPE norm_dl = norm(dl);
    PRECISION_TYPE norm_rf = norm(rf);
    PRECISION_TYPE norm_ri = norm(ri);

    if (norm_dl < TOO_SMALL) return; // Avoid division by zero

    // Compute cross product (rf × dl)
    PRECISION_TYPE cross_rf_dl[3];
    cross_product(rf, dl, cross_rf_dl);

    // Normalize cross product by norm_dl
    scale(cross_rf_dl, 1.0 / norm_dl, cross_rf_dl);
    PRECISION_TYPE distance = norm(cross_rf_dl);

    // Compute cosines
    PRECISION_TYPE cos_theta_i = dot_product(rf, dl) / (norm_rf * norm_dl);
    PRECISION_TYPE cos_theta_f = dot_product(ri, dl) / (norm_ri * norm_dl);

    // Compute Biot-Savart contribution
    PRECISION_TYPE absB1 = (distance > TOO_SMALL) ? (cos_theta_i - cos_theta_f) / distance : 0.0;

    // Compute direction (rf × dl) and normalize
    PRECISION_TYPE dir[3];
    cross_product(rf, dl, dir);
    PRECISION_TYPE norm_dir = norm(dir);

    absB1 = (norm_dir > TOO_SMALL) ? absB1 / norm_dir : 0.0;

    // Update B-field
    PRECISION_TYPE scaled_dir[3];
    scale(dir, absB1, scaled_dir);
    add(B, scaled_dir);
}

void calc_bfield_parallel_se(const uint32_t num_seg,
    const PRECISION_TYPE* seg_start,
    const PRECISION_TYPE* seg_end,
    const uint32_t num_grid,
    const PRECISION_TYPE* grid_xyz,
    PRECISION_TYPE* bfield_xyz) {
    // Reset magnetic field grid
    fill(bfield_xyz, bfield_xyz + 3 * num_grid, 0.0);

    cout << "Starting B-field Calculation (SE)..." << endl;
    auto start_time = chrono::high_resolution_clock::now();

    for (int seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const PRECISION_TYPE* start_point = &seg_start[3 * seg_indx];
        const PRECISION_TYPE* end_point = &seg_end[3 * seg_indx];

        for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
            biot_savart_se(start_point, end_point, &grid_xyz[3 * grid_indx], &bfield_xyz[3 * grid_indx]);
        }

        if (seg_indx % 10 == 0) {
            cout << "Progress (SE): " << seg_indx << " / " << num_seg << endl;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Simulation finished (SE). Time elapsed: " << elapsed.count() << " seconds" << endl;
}

// ---------------------------------------------------------












//void biot_savart_chukman(const vector<PRECISION_TYPE>& s0,
//                                const vector<PRECISION_TYPE>& s1,
//                                const vector<PRECISION_TYPE>& r,
//                                vector<PRECISION_TYPE>& B) {
//    
//        vector<PRECISION_TYPE> i, k0, k1, A;
//        PRECISION_TYPE norm_i, norm_k0, norm_k1, norm_A_sqr, a0;
//    
//        i = subtract(s1, s0);
//        norm_i = norm(i);
//        if (norm_i < TOO_SMALL) return;
//        norm_i = 1. / norm_i;
//        i = scale(i, norm_i);
//    
//        k0 = subtract(s0, r);
//        k1 = subtract(s1, r);
//        norm_k0 = norm(k0);
//        norm_k1 = norm(k1);
//        if (norm_k0 < TOO_SMALL) return; norm_k0 = 1. / norm_k0;
//        if (norm_k1 < TOO_SMALL) return; norm_k1 = 1. / norm_k1;
//    
//        A = cross_product(k0, i);
//        norm_A_sqr = (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
//        if (norm_A_sqr < TOO_SMALL) return;
//        norm_A_sqr = 1. / norm_A_sqr;
//    
//        A = scale(A, norm_A_sqr);
//        k0 = scale(k0, norm_k0);
//        k1 = scale(k1, norm_k1);
//    
//        a0 = dot_product(subtract(k1, k0), i);
//        //B = scale(A, a0);
//        B = add(B, scale(A, a0));
//    }

//void biot_savart_se(const vector<PRECISION_TYPE>& start,
//                    const vector<PRECISION_TYPE>& end,
//                    const vector<PRECISION_TYPE>& grid,
//                    vector<PRECISION_TYPE>& B) {
//
//    vector<PRECISION_TYPE> rf = subtract(end, grid); // r_(i+1)    
//    vector<PRECISION_TYPE> ri = subtract(start, grid); // r_i
//    vector<PRECISION_TYPE> dl = subtract(end, start); // dl
//
//    PRECISION_TYPE norm_dl = norm(dl);
//    PRECISION_TYPE norm_rf = norm(rf);
//    PRECISION_TYPE norm_ri = norm(ri);
//
//    vector<PRECISION_TYPE> dist = cross_product(rf, dl); // r cross dl
//    dist = scale(dist, 1. / norm_dl);
//
//    PRECISION_TYPE distance = norm(dist); // r
//
//    PRECISION_TYPE cos_theta_i = dot_product(rf, dl) / (norm_rf * norm_dl);
//    PRECISION_TYPE cos_theta_f = dot_product(ri, dl) / (norm_ri * norm_dl);
//
//    PRECISION_TYPE absB1 = distance > TOO_SMALL ? (cos_theta_i - cos_theta_f) / distance : 0.;
//
//    vector<PRECISION_TYPE> dir = cross_product(rf, dl);
//    PRECISION_TYPE norm_dir = norm(dir);
//
//    absB1 = norm_dir > TOO_SMALL ? absB1 / norm_dir : 0.;
//
//    B = add(B, scale(dir, absB1));
//}
//
//
//
//void calc_bfield_parallel_se(const uint32_t num_seg,
//                            const vector<vector<PRECISION_TYPE>>& seg_start,
//                            const vector<vector<PRECISION_TYPE>>& seg_end,
//                            const uint32_t num_grid,
//                            const vector<vector<PRECISION_TYPE>>& grid_xyz,
//                            vector<vector<PRECISION_TYPE>>& bfield_xyz) {
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
//        //#pragma omp parallel for
//        for (int grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//            biot_savart_se(start_point, end_point, grid_xyz[grid_indx], bfield_xyz[grid_indx]);
//        }
//
//        if (seg_indx % 10 == 0) {
//            cout << "Progress (SE): " << seg_indx << " / " << num_seg << endl;
//        }
//    }
//
//
//    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
//    chrono::duration<double> elapsed = end_time - start_time;
//    cout << "Simulation finished (SE). Time elapsed: " << elapsed.count() << " seconds" << endl;
//}







//void calc_bfield_parallel_chukman(  const uint32_t num_seg,
//                                    const vector<vector<PRECISION_TYPE>>& seg_start,
//                                    const vector<vector<PRECISION_TYPE>>& seg_end,
//                                    const uint32_t num_grid,
//                                    const vector<vector<PRECISION_TYPE>>& grid_xyz,
//                                    vector<vector<PRECISION_TYPE>>& bfield_xyz) {
//    // Reset magnetic field grid
//    for (auto& b : bfield_xyz) {
//        b = { 0.0, 0.0, 0.0 };
//    }
//
//    cout << "Starting B-field Calculation (Chukman)..." << endl;
//    auto start_time = chrono::high_resolution_clock::now();  // Start timer
//
//    //#pragma omp parallel for
//    for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
//        const auto& start_point = seg_start[seg_indx];
//        const auto& end_point = seg_end[seg_indx];
//
//        try {
//            // Parallel loop over grid points
//            //#pragma omp parallel for
//            for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//                biot_savart_chukman(start_point, end_point, grid_xyz[grid_indx], bfield_xyz[grid_indx]);
//            }
//        }
//        catch (exception& ex) {
//            cout << "Simulation failed: " << ex.what() << endl;
//        }
//
//        if (seg_indx % 10 == 0) {
//            cout << "Progress (Chukman): " << seg_indx << " / " << num_seg << endl;
//        }
//    }
//
//    auto end_time = chrono::high_resolution_clock::now();  // Stop timer
//    chrono::duration<double> elapsed = end_time - start_time;
//    cout << "Simulation finished (Chukman). Time elapsed: " << elapsed.count() << " seconds" << endl;
//}
