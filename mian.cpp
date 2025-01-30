using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

const long double PI = 3.141592653589793238462643383279502884197169399375105820974944L;
const long double MU0 = 1.2566370612720e-6L;
const long double TOO_SMALL = 1e-12L;


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



void generate_uniform_grid(vector<vector<long double>>& grid_xyz,
	long double x_min, long double x_max,
	long double y_min, long double y_max,
	long double z_min, long double z_max,
	uint32_t num_points) {

	grid_xyz.clear();
	grid_xyz.reserve(num_points * num_points * num_points);

	long double x_step = (x_max - x_min) / (num_points - 1);
	long double y_step = (y_max - y_min) / (num_points - 1);
	long double z_step = (z_max - z_min) / (num_points - 1);

	for (uint32_t i = 0; i < num_points; i++) {
		long double x = x_min + i * x_step;
		for (uint32_t j = 0; j < num_points; j++) {
			long double y = y_min + j * y_step;
			for (uint32_t k = 0; k < num_points; k++) {
				long double z = z_min + k * z_step;
				grid_xyz.push_back({ x, y, z });
			}
		}
	}
}





void biot_savart_chukman(const vector<long double>& s0,
	const vector<long double>& s1,
	const vector<long double>& r,
	vector<long double>& B) { //field[T] of a wire segment from s0[mm] to s1[mm] with current 1e4 A, at observation point r[mm]
	//field[T] of a wire segment from s0[mm] to s1[mm] with current 1e4 A, at observation point r[mm]
	vector<long double> i, k0, k1, A;
	long double norm_i, norm_k0, norm_k1, norm_A_sqr, a0;

	// Current vector
	i = subtract(s1, s0);
	norm_i = norm(i);
	if (norm_i < TOO_SMALL) return;
	norm_i = 1. / norm_i;
	i = scale(i, norm_i);


	// Position vectors
	k0 = subtract(s0, r);
	k1 = subtract(s1, r);
	norm_k0 = norm(k0);
	norm_k1 = norm(k1);
	if (norm_k0 < TOO_SMALL) return; norm_k0 = 1. / norm_k0;
	if (norm_k1 < TOO_SMALL) return; norm_k1 = 1. / norm_k1;

	// Cross product
	A = cross_product(k0, i);
	norm_A_sqr = (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
	if (norm_A_sqr < TOO_SMALL) return;
	norm_A_sqr = 1. / norm_A_sqr;

	// Normalize
	A = scale(A, norm_A_sqr);
	k0 = scale(k0, norm_k0);
	k1 = scale(k1, norm_k1);


	// B field
	//a0 = (k1[0] - k0[0]) * i[0] + (k1[1] - k0[1]) * i[1] + (k1[2] - k0[2]) * i[2];
	a0 = dot_product(subtract(k1, k0), i);
	B = scale(A, a0);

	return;
};


void biot_savart_se(const vector<long double>& start,
	const vector<long double>& end,
	const vector<long double>& grid,
	vector<long double>& b1)
{
	vector<long double> rf = subtract(end, grid); // r_(i+1)	
	vector<long double> ri = subtract(start, grid); // r_i
	vector<long double> dl = subtract(end, start); // dl

	long double norm_dl = norm(dl);
	long double norm_rf = norm(rf);
	long double norm_ri = norm(ri);

	// distance from grid to the wire segment = rf * cos(phi)
	vector<long double> dist = cross_product(rf, dl); // r cross dl
	dist = scale(dist, 1. / norm_dl);


	long double distance = norm(dist); // r

	long double cos_theta_i = dot_product(rf, dl) / (norm_rf * norm_dl);
	long double cos_theta_f = dot_product(ri, dl) / (norm_ri * norm_dl);


	// absolute value of B1
	long double absB1 = distance > TOO_SMALL ? (cos_theta_i - cos_theta_f) / distance : 0.;

	vector<long double> dir = cross_product(rf, dl);
	long double norm_dir = norm(dir);

	absB1 = norm_dir > TOO_SMALL ? absB1 / norm_dir : 0.;

	b1 = add(b1, scale(dir, absB1));

}


// ---------------------------------------------------------
// C interface
// ---------------------------------------------------------

extern "C" {
	bool calculate_b0(const uint32_t num_seg,						// number of wire segments
		const vector<vector<long double>>& seg_start,	// start of segments
		const vector<vector<long double>>& seg_end,		// end of segments
		const uint32_t num_grid,					// number of grid points
		const vector<vector<long double>>& grid_xyz,     // the grid
		vector<vector<long double>>& b1_xyz              // Output B-field grid
	)

	{
		// reset b1
		for (auto& b : b1_xyz) {
			b = { 0.0, 0.0, 0.0 };
		}


		// loop over wire segments
		auto start = chrono::steady_clock::now();
		for (uint32_t seg_indx = 0; seg_indx < num_seg; seg_indx++) {
			const auto& start_point = seg_start[seg_indx];
			const auto& end_point = seg_end[seg_indx];

			try
			{
				//#pragma omp parallel for
				for (uint32_t grid_indx = 0; grid_indx < num_grid; grid_indx++)
					biot_savart_se(start_point, end_point, grid_xyz[grid_indx], b1_xyz[grid_indx]);
			}
			catch (exception& ex)
			{
				cout << "Simulation failed." << endl;
				cout << ex.what() << endl;
				return false;
			}

			// Print progress
			if (seg_indx % 100 == 0)
				cout << "Progress: " << seg_indx << " / " << num_seg << endl;
		}



		//for (auto& b : b1_xyz) {
		//	b[0] *= MU0 / (4.0 * PI);
		//	b[1] *= MU0 / (4.0 * PI);
		//	b[2] *= MU0 / (4.0 * PI);
		//}


		cout << "Simulation finished in " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count() / 1000. << " second(s)." << endl;

		return true;
	}


}

// Function to generate a circular loop
void generate_circular_loop(vector<vector<long double>>& segment_start, vector<vector<long double>>& segment_end, long double radius, int num_segments) {
	segment_start.clear();
	segment_end.clear();

	long double dtheta = 2.0 * PI / num_segments;

	for (int i = 0; i < num_segments; i++) {
		long double theta1 = i * dtheta;
		long double theta2 = (i + 1) * dtheta;

		vector<long double> s0 = { radius * cos(theta1), radius * sin(theta1), 0.0 }; // Start of segment
		vector<long double> s1 = { radius * cos(theta2), radius * sin(theta2), 0.0 }; // End of segment

		segment_start.push_back(s0);
		segment_end.push_back(s1);
	}
}

int main() {
	cout << fixed << setprecision(16);


	// A simple wire in the z-axis
	//int num_segments = 10;
	//vector<vector<long double>> segment_start, segment_end;
	//long double zmin = -10.0, zmax = 10.0;
	//long double dz = (zmax - zmin) / num_segments;

	//for (int seg = 0; seg < num_segments; seg++) {
	//	segment_start.push_back({ 0.0, 0.0, zmin + seg * dz });
	//	segment_end.push_back({ 0.0, 0.0, zmin + (seg + 1) * dz });
	//}
	long double loop_radius = 2.0;
	int num_segments = 100;
	vector<vector<long double>> segment_start, segment_end;
	generate_circular_loop(segment_start, segment_end, loop_radius, num_segments);


	// Uniform Grid
	vector<vector<long double>> grid;
	generate_uniform_grid(grid, -10.0, 10.0, -10.0, 10.0, -10.0, 10.0, 30);



	ofstream output_file_chukman("B_field_loop_chukman.txt");
	ofstream output_file_se("B_field_loop_se.txt");
	//ofstream output_file_chukman("B_field_wire_chukman.txt");
	//ofstream output_file_se("B_field_wire_se.txt");

	output_file_chukman.precision(20);
	output_file_se.precision(20);

	if (!output_file_chukman.is_open()) {
		cerr << "Error: Could not open file for writing!" << endl;
		return 1;
	}

	if (!output_file_se.is_open()) {
		cerr << "Error: Could not open file for writing!" << endl;
		return 1;
	}

	output_file_chukman << "x y z Bx By Bz\n";
	output_file_se << "x y z Bx By Bz\n";


	size_t total_points = grid.size();
	size_t count = 0;
	size_t update_interval = max(total_points / 10, size_t(1)); // Update every 10% or at least once

	// Call calculate_b0 for the entire grid
	vector<vector<long double>> b1_xyz(grid.size(), vector<long double>(3, 0.0));  // Initialize B-field grid
	calculate_b0(num_segments, segment_start, segment_end, grid.size(), grid, b1_xyz);
	//cout << "First point B-field: " << b1_xyz[0][0] << " " << b1_xyz[0][1] << " " << b1_xyz[0][2] << endl; 
	//cout << "second point B-field: " << b1_xyz[1][0] << " " << b1_xyz[1][1] << " " << b1_xyz[1][2] << endl;


	for (const auto& r : grid) {
		vector<long double> B_chukman = { 0.0, 0.0, 0.0 }; // Initialize B-field for Chukman routine

		for (size_t i = 0; i < segment_start.size(); i++) {
			vector<long double> B_temp(3, 0.0);
			biot_savart_chukman(segment_start[i], segment_end[i], r, B_temp);
			B_chukman = add(B_chukman, B_temp);
		}

		// Extract the corresponding magnetic field from the b1_xyz grid for calculate_b0
		size_t grid_index = &r - &grid[0];  // Calculate the index of the current grid point
		vector<long double> B_se = b1_xyz[grid_index];  // Extract the magnetic field at the current grid point

		// Write to Chukman output file: x y z Bx By Bz
		output_file_chukman << r[0] << " " << r[1] << " " << r[2] << " "
			<< B_chukman[0] << " " << B_chukman[1] << " " << B_chukman[2] << "\n";

		// Write to calculate_b0 output file: x y z Bx By Bz
		output_file_se << r[0] << " " << r[1] << " " << r[2] << " "
			<< B_se[0] << " " << B_se[1] << " " << B_se[2] << "\n";

		// Progress update
		count++;
		if (count % update_interval == 0) {
			cout << "Progress: " << (100 * count / total_points) << "% ("
				<< count << "/" << total_points << " points processed)" << endl;
		}
	}


	output_file_chukman.close();
	output_file_se.close();

	cout << "Done!" << endl;

}


// ---------------------------------------------------------
//void compute_b_field_on_grid(const vector<long double>& grid, std::ofstream& file) {
//	int N = 100;
//	long double R = 1.0;
//	long double current = 1.0;
//
//	for (size_t i = 0; i < grid_xyz.size(); i += 3) {
//		std::vector<long double> r = { grid_xyz[i], grid_xyz[i + 1], grid_xyz[i + 2] };
//		std::vector<long double> B = { 0.0, 0.0, 0.0 };
//
//		for (int j = 0; j < N; j++) {
//			long double theta0 = 2.0 * M_PI * j / N;
//			long double theta_i = 2.0 * M_PI * (j + 1) / N;
//
//			std::vector<long double> s0 = { R * cos(theta0), R * sin(theta0), 0.0 };
//			std::vector<long double> s1 = { R * cos(theta_i), R * sin(theta_i), 0.0 };
//
//			std::vector<long double> dB = biot_savart(r, s0, s1, current);
//			B[0] += dB[0];
//			B[1] += dB[1];
//			B[2] += dB[2];
//		}
//
//		file << r[0] << " " << r[1] << " " << r[2] << " "
//			<< B[0] << " " << B[1] << " " << B[2] << "\n";
//	}
//}

//void compute_b_field_on_grid(const vector<long double>& grid_xyz, ofstream& file) {
//	int N = 100;
//	long double R = 1.0;
//	long double current = 1.0;
//
//	for (size_t i = 0; i < grid_xyz.size(); i += 3) {
//		vector<long double> r = { grid_xyz[i], grid_xyz[i + 1], grid_xyz[i + 2] };
//		vector<long double> B = { 0.0, 0.0, 0.0 };
//
//		for (int j = 0; j < N; j++) {
//			long double theta0 = 2.0 * M_PI * j / N;
//			long double theta_i = 2.0 * M_PI * (j + 1) / N;
//
//			vector<long double> s0 = { R * cos(theta0), R * sin(theta0), 0.0 };
//			vector<long double> s1 = { R * cos(theta_i), R * sin(theta_i), 0.0 };
//
//			vector<long double> dB = biot_savart_chukman(r, s0, s1, current);
//			B[0] += dB[0];
//			B[1] += dB[1];
//			B[2] += dB[2];
//		}
//
//		file << r[0] << " " << r[1] << " " << r[2] << " "
//			<< B[0] << " " << B[1] << " " << B[2] << "\n";
//	}
//}