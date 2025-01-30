using namespace std;

#include <iostream>
#include <functional>

#include "vector_utils.h"  // Include the header where vector operations are defined FIRST
#include "biot_savart_routines.h"


void biot_savart_chukman(const vector<long double>& s0,
    const vector<long double>& s1,
    const vector<long double>& r,
    vector<long double>& B) {
    // Implementation of the Biot-Savart Chukman method (as you provided)
    vector<long double> i, k0, k1, A;
    long double norm_i, norm_k0, norm_k1, norm_A_sqr, a0;

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
    B = scale(A, a0);
}

void biot_savart_se(const vector<long double>& start,
    const vector<long double>& end,
    const vector<long double>& grid,
    vector<long double>& b1) {
    // Implementation of the Biot-Savart SE method (as you provided)
    vector<long double> rf = subtract(end, grid); // r_(i+1)    
    vector<long double> ri = subtract(start, grid); // r_i
    vector<long double> dl = subtract(end, start); // dl

    long double norm_dl = norm(dl);
    long double norm_rf = norm(rf);
    long double norm_ri = norm(ri);

    vector<long double> dist = cross_product(rf, dl); // r cross dl
    dist = scale(dist, 1. / norm_dl);

    long double distance = norm(dist); // r

    long double cos_theta_i = dot_product(rf, dl) / (norm_rf * norm_dl);
    long double cos_theta_f = dot_product(ri, dl) / (norm_ri * norm_dl);

    long double absB1 = distance > TOO_SMALL ? (cos_theta_i - cos_theta_f) / distance : 0.;

    vector<long double> dir = cross_product(rf, dl);
    long double norm_dir = norm(dir);

    absB1 = norm_dir > TOO_SMALL ? absB1 / norm_dir : 0.;

    b1 = add(b1, scale(dir, absB1));
}


//Unified wrapper for both Biot-Savart routines (SE and Chukman)
void calc_bfield_parallel(const uint32_t num_seg,
    const vector<vector<long double>>& seg_start,
    const vector<vector<long double>>& seg_end,
    const uint32_t num_grid,
    const vector<vector<long double>>& grid_xyz,
    vector<vector<long double>>& b1_xyz,
    function<void(const vector<long double>&,
        const vector<long double>&,
        const vector<long double>&,
        vector<long double>&)> biot_savart_method) {
    // Reset magnetic field grid
    for (auto& b : b1_xyz) {
        b = { 0.0, 0.0, 0.0 };
    }

    //#pragma omp parallel for
    for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
        const auto& start_point = seg_start[seg_indx];
        const auto& end_point = seg_end[seg_indx];

        try {
            // Parallel loop over grid points
            #pragma omp parallel for
            for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++) {
                biot_savart_method(start_point, end_point, grid_xyz[grid_indx], b1_xyz[grid_indx]);
            }
        }
        catch (exception& ex) {
            cout << "Simulation failed: " << ex.what() << endl;
        }

        if (seg_indx % 100 == 0) {
            cout << "Progress: " << seg_indx << " / " << num_seg << endl;
        }
    }

    cout << "Simulation finished." << endl;
}

//void calc_bfield_parallel(const uint32_t num_seg,
//    const vector<vector<long double>>& seg_start,
//    const vector<vector<long double>>& seg_end,
//    const uint32_t num_grid,
//    const vector<vector<long double>>& grid_xyz,
//    vector<vector<long double>>& b1_xyz,
//    function<void(const vector<long double>&,
//        const vector<long double>&,
//        const vector<long double>&,
//        vector<long double>&)> biot_savart_method) {
//    // Reset magnetic field grid
//    for (auto& b : b1_xyz) {
//        b = { 0.0, 0.0, 0.0 };
//    }
//
//    #pragma omp parallel for
//    for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
//        const auto& start_point = seg_start[seg_indx];
//        const auto& end_point = seg_end[seg_indx];
//
//        try {
//            // Parallel loop over grid points
//            #pragma omp parallel for
//            for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++) {
//                // Create a local temporary vector for each thread
//                vector<long double> B_temp(3, 0.0);
//
//                // Call the Biot-Savart method to compute the field for this segment and grid point
//                biot_savart_method(start_point, end_point, grid_xyz[grid_indx], B_temp);
//
//                // Accumulate the results in the global b1_xyz
//                #pragma omp atomic
//                b1_xyz[grid_indx][0] += B_temp[0];
//                #pragma omp atomic
//                b1_xyz[grid_indx][1] += B_temp[1];
//                #pragma omp atomic
//                b1_xyz[grid_indx][2] += B_temp[2];
//            }
//        }
//        catch (exception& ex) {
//            cout << "Simulation failed: " << ex.what() << endl;
//        }
//
//        if (seg_indx % 100 == 0) {
//            cout << "Progress: " << seg_indx << " / " << num_seg << endl;
//        }
//    }
//
//    cout << "Simulation finished." << endl;
//}


//extern "C" {
//    bool calc_bfield_se(const uint32_t num_seg,
//        const vector<vector<long double>>& seg_start,
//        const vector<vector<long double>>& seg_end,
//        const uint32_t num_grid,
//        const vector<vector<long double>>& grid_xyz,
//        vector<vector<long double>>& b1_xyz) {
//        // Implementation of the calculation of B-field over grid points
//        for (auto& b : b1_xyz) {
//            b = { 0.0, 0.0, 0.0 };
//        }
//
//        for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
//            const auto& start_point = seg_start[seg_indx];
//            const auto& end_point = seg_end[seg_indx];
//
//            try {
//                for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++)
//                    biot_savart_se(start_point, end_point, grid_xyz[grid_indx], b1_xyz[grid_indx]);
//            }
//            catch (exception& ex) {
//                cout << "Simulation failed: " << ex.what() << endl;
//                return false;
//            }
//
//            if (seg_indx % 100 == 0) {
//                cout << "Progress: " << seg_indx << " / " << num_seg << endl;
//            }
//        }
//
//        cout << "Simulation finished." << endl;
//        return true;
//    }
//}
