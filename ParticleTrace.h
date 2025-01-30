#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

class ParticleTrace {
private:
    vector<vector<long double>> segment_start;
    vector<vector<long double>> segment_end;
    size_t num_segments;

public:
    ParticleTrace(const string& filename);
    void readCSV(const string& filename);

    const vector<vector<long double>>& getSegmentStart() const { return segment_start; }
    const vector<vector<long double>>& getSegmentEnd() const { return segment_end; }
    const size_t getNumSegments() const { return num_segments; }
};

ParticleTrace::ParticleTrace(const string& filename) {
    readCSV(filename);
}

void ParticleTrace::readCSV(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;
    vector<vector<long double>> positions;

    // Read the first line (header)
    getline(file, line);

    // Read the data
    while (getline(file, line)) {
        stringstream ss(line);

        vector<long double> row;
        string value;
        int col_indx = 0;

        while (getline(ss, value, ',')) {
            if (col_indx == 3) {
                break;
            }
            row.push_back(stold(value));
            col_indx++;
        }

        if (row.size() == 3) {
            positions.push_back(row);
        }

    }
    file.close();

    for (size_t i = 0; i < positions.size() - 1; i++) {
        segment_start.push_back(positions[i]);
        segment_end.push_back(positions[i + 1]);
    }

    num_segments = segment_start.size();

}

