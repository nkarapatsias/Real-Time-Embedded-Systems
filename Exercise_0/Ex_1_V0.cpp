// g++ Ex_1_Final.cpp -o Ex_1_Final.exe -I"./include" -L"./lib" -lopenblas -O3

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cblas.h>

using namespace std;

// Υπολογισμός του αθροισματος των τετραγωνων
void compute_squared_norms(const double* matrix, int rows, int cols, vector<double>& norms) {
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < cols; ++j) {
            double val = matrix[i * cols + j];
            sum += val * val;
        }
        norms[i] = sum;
    }
}


//1. Ακριβής Υπολογισμός
void knnsearch_exact(const double* C, const double* Q, int N, int M, int d, int k, 
                     vector<int>& idx, vector<double>& dst) {
    
    vector<double> C_sq(N);
    compute_squared_norms(C, N, d, C_sq);

    idx.resize(M * k);
    dst.resize(M * k);

    int BLOCK_SIZE = 1000; 

    for (int b_start = 0; b_start < M; b_start += BLOCK_SIZE) {
        int current_B = min(BLOCK_SIZE, M - b_start);
        
        vector<double> Q_block_sq(current_B);
        compute_squared_norms(Q + b_start * d, current_B, d, Q_block_sq);

        vector<double> D_block(current_B * N, 0.0);

        // OpenBLAS: Υπολογισμός -2 * Q_block * C^T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                    current_B, N, d, 
                    -2.0, Q + b_start * d, d, C, d, 
                    0.0, D_block.data(), N);

        for (int i = 0; i < current_B; ++i) {
            vector<pair<double, int>> row_dists(N);
            for (int j = 0; j < N; ++j) {
                double dist_sq = Q_block_sq[i] + C_sq[j] + D_block[i * N + j];
                dist_sq = max(0.0, dist_sq); // Αποτροπή αρνητικών (precision issues)
                row_dists[j] = {sqrt(dist_sq), j};
            }

            // Quick-select
            nth_element(row_dists.begin(), row_dists.begin() + k, row_dists.end());
            sort(row_dists.begin(), row_dists.begin() + k);

            int global_i = b_start + i;
            for (int j = 0; j < k; ++j) {
                dst[global_i * k + j] = row_dists[j].first;
                // +1 ΓΙΑ ΣΥΜΒΑΤΟΤΗΤΑ ΜΕ MATLAB (1-based indexing)
                idx[global_i * k + j] = row_dists[j].second + 1; 
            }
        }
    }
}

// Βοηθητική συνάρτηση για Ευκλείδεια απόσταση
inline double euclidean_dist(const double* a, const double* b, int d) {
    double sum = 0;
    for(int i = 0; i < d; i++) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrt(max(0.0, sum));
}


//Βοηθητική συνάρτηση ενημέρωσης των top-k γειτόνων ενός σημείου
void update_knn(vector<pair<double, int>>& knn, double dist, int idx, int k) {
    // Αποτροπή διπλοεγγραφών
    for(const auto& p : knn) {
        if(p.second == idx) return;
    }

    if (knn.size() < k) {
        knn.push_back({dist, idx});
        sort(knn.begin(), knn.end());
    } else if (dist < knn.back().first) {
        knn.back() = {dist, idx};
        sort(knn.begin(), knn.end());
    }
}

// 2. Προσεγγιστικός Αναδρομικός Υπολογισμός (KD-Tree Style Partitioning)
void approx_recursive_helper(const double* C_full, vector<int>& subset, int d, int k, int depth,
                             vector<vector<pair<double, int>>>& global_neighbors) {
    
    int N_sub = subset.size();
    
    if (N_sub <= 200) {
        for (int i = 0; i < N_sub; ++i) {
            int q_idx = subset[i];
            for (int j = 0; j < N_sub; ++j) {
                if (i == j) continue; 
                int c_idx = subset[j];
                double dist = euclidean_dist(C_full + q_idx * d, C_full + c_idx * d, d);
                update_knn(global_neighbors[q_idx], dist, c_idx, k);
            }
        }
        return;
    }

    int split_dim = depth % d; 
    
    auto comparator = [&](int a, int b) {
        return C_full[a * d + split_dim] < C_full[b * d + split_dim];
    };

    int mid = N_sub / 2;
    nth_element(subset.begin(), subset.begin() + mid, subset.end(), comparator);
    
    double median_val = C_full[subset[mid] * d + split_dim];

    vector<int> left_subset(subset.begin(), subset.begin() + mid);
    vector<int> right_subset(subset.begin() + mid, subset.end());

    approx_recursive_helper(C_full, left_subset, d, k, depth + 1, global_neighbors);
    approx_recursive_helper(C_full, right_subset, d, k, depth + 1, global_neighbors);
    
    for (int q_idx : left_subset) {
        double dist_to_plane = abs(C_full[q_idx * d + split_dim] - median_val);
        
        // Αν δεν έχουμε βρει k γείτονες ή αν η απόσταση στο σύνορο είναι μικρότερη από τον k-οστό γείτονα
        if (global_neighbors[q_idx].size() < k || dist_to_plane <= global_neighbors[q_idx].back().first) {
            for (int c_idx : right_subset) {
                // Προσεγγιστικό φίλτρο: Ελέγχουμε μόνο σημεία του right που είναι κοντά στο σύνορο
                double dist_to_plane_c = abs(C_full[c_idx * d + split_dim] - median_val);
                if (dist_to_plane_c <= global_neighbors[q_idx].back().first) {
                    double dist = euclidean_dist(C_full + q_idx * d, C_full + c_idx * d, d);
                    update_knn(global_neighbors[q_idx], dist, c_idx, k);
                }
            }
        }
    }

    for (int q_idx : right_subset) {
        double dist_to_plane = abs(C_full[q_idx * d + split_dim] - median_val);
        
        if (global_neighbors[q_idx].size() < k || dist_to_plane <= global_neighbors[q_idx].back().first) {
            for (int c_idx : left_subset) {
                double dist_to_plane_c = abs(C_full[c_idx * d + split_dim] - median_val);
                if (dist_to_plane_c <= global_neighbors[q_idx].back().first) {
                    double dist = euclidean_dist(C_full + q_idx * d, C_full + c_idx * d, d);
                    update_knn(global_neighbors[q_idx], dist, c_idx, k);
                }
            }
        }
    }
}


// Συνάρτηση εκκίνησης για την Προσεγγιστική, Αναδρομική Λύση (C == Q)

void knnsearch_approx_recursive(const double* C, int N, int d, int k, 
                                vector<int>& idx, vector<double>& dst) {
    idx.resize(N * k);
    dst.resize(N * k);
    
    vector<int> initial_indices(N);
    iota(initial_indices.begin(), initial_indices.end(), 0); 
    
    vector<vector<pair<double, int>>> global_neighbors(N);

    approx_recursive_helper(C, initial_indices, d, k, 0, global_neighbors);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < k && j < global_neighbors[i].size(); ++j) {
            dst[i * k + j] = global_neighbors[i][j].first;
            idx[i * k + j] = global_neighbors[i][j].second + 1; // 1-based indexing for MATLAB
        }
    }
}

int main() {
    int N = 5000; // Dataset μέγεθος 
    int d = 64;   // Διαστάσεις
    int k = 5;    // Κοντινότεροι γείτονες

    cout << "Generating Mock Data..." << endl;
    vector<double> C(N * d);
    for (double& val : C) val = (rand() % 100) / 10.0;

    // ----------------------------------------------------
    // EXACT KNN TEST 
    // ----------------------------------------------------
    vector<int> exact_indices;
    vector<double> exact_distances;

    cout << "\nRunning Exact k-NN (Block-wise OpenBLAS)..." << endl;
    // Εδώ περνάμε το C και στις δύο θέσεις C, Q
    knnsearch_exact(C.data(), C.data(), N, N, d, k, exact_indices, exact_distances);

    cout << "Exact - Nearest neighbors for Point 0:\n";
    for(int j = 0; j < k; ++j) {
        cout << "Neighbor " << j+1 << " -> Index: " << exact_indices[j] 
             << " | Distance: " << exact_distances[j] << "\n";
    }

    // ----------------------------------------------------
    // APPROXIMATE KNN TEST
    // ----------------------------------------------------
    vector<int> approx_indices;
    vector<double> approx_distances;

    cout << "\nRunning Approximate Recursive k-NN (Divide & Conquer)..." << endl;
    knnsearch_approx_recursive(C.data(), N, d, k, approx_indices, approx_distances);

    cout << "Approximate - Nearest neighbors for Point 0:\n";
    for(int j = 0; j < k; ++j) {
        cout << "Neighbor " << j+1 << " -> Index: " << approx_indices[j] 
             << " | Distance: " << approx_distances[j] << "\n";
    }

    return 0;
}

