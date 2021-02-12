#include <algorithm>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <random>
#include <vector>
#include "utils/FileUtil.h"
#include "utils/Timer.hpp"

using namespace std;

typedef unsigned long ulong;

/*
 *  Internal implementation of the SMAWK algorithm.
 */
template<typename T>
void _smawk(
        const vector<ulong> &rows,
        const vector<ulong> &cols,
        const function<T(ulong, ulong)> &lookup,
        vector<ulong> *result) {
    // Recursion base case
    if (rows.size() == 0) return;

    // ********************************
    // * REDUCE
    // ********************************

    vector<ulong> _cols;  // Stack of surviving columns
    for (ulong col : cols) {
        while (true) {
            if (_cols.size() == 0) break;
            ulong row = rows[_cols.size() - 1];
            if (lookup(row, col) >= lookup(row, _cols.back()))
                break;
            _cols.pop_back();
        }
        if (_cols.size() < rows.size())
            _cols.push_back(col);
    }

    // Call recursively on odd-indexed rows
    vector<ulong> odd_rows;
    for (ulong i = 1; i < rows.size(); i += 2) {
        odd_rows.push_back(rows[i]);
    }
    _smawk(odd_rows, _cols, lookup, result);

    unordered_map<ulong, ulong> col_idx_lookup;
    for (ulong idx = 0; idx < _cols.size(); ++idx) {
        col_idx_lookup[_cols[idx]] = idx;
    }

    // ********************************
    // * INTERPOLATE
    // ********************************

    // Fill-in even-indexed rows
    ulong start = 0;
    for (ulong r = 0; r < rows.size(); r += 2) {
        ulong row = rows[r];
        ulong stop = _cols.size() - 1;
        if (r < rows.size() - 1)
            stop = col_idx_lookup[(*result)[rows[r + 1]]];
        ulong argmin = _cols[start];
        T min = lookup(row, argmin);
        for (ulong c = start + 1; c <= stop; ++c) {
            T value = lookup(row, _cols[c]);
            if (c == start || value < min) {
                argmin = _cols[c];
                min = value;
            }
        }
        (*result)[row] = argmin;
        start = stop;
    }
}

/*
 *  Interface for the SMAWK algorithm, for finding the minimum value in each row
 *  of an implicitly-defined totally monotone matrix.
 */
template<typename T>
vector<ulong> smawk(
        const ulong num_rows,
        const ulong num_cols,
        const function<T(ulong, ulong)> &lookup) {
    vector<ulong> result;
    result.resize(num_rows);
    vector<ulong> rows(num_rows);
    iota(begin(rows), end(rows), 0);
    vector<ulong> cols(num_cols);
    iota(begin(cols), end(cols), 0);
    _smawk<T>(rows, cols, lookup, &result);
    return result;
}

/*
 *  Calculates cluster costs in O(1) using prefix sum arrays.
 */
template<class DT>
class CostCalculator {
    vector<double> cumsum;
    vector<double> cumsum2;

public:
    CostCalculator(const vector<DT> &vec, ulong n) {
        cumsum.push_back(0.0);
        cumsum2.push_back(0.0);
        for (ulong i = 0; i < n; ++i) {
            double x = vec[i];
            cumsum.push_back(x + cumsum[i]);
            cumsum2.push_back(x * x + cumsum2[i]);
        }
    }

    double calc(ulong i, ulong j) {
        if (j < i) return 0.0;
        double mu = (cumsum[j + 1] - cumsum[i]) / (j - i + 1);
        double result = cumsum2[j + 1] - cumsum2[i];
        result += (j - i + 1) * (mu * mu);
        result -= (2 * mu) * (cumsum[j + 1] - cumsum[i]);
        return result;
    }
};

template<typename T>
class Matrix {
    vector<T> data;
    ulong num_rows;
    ulong num_cols;

public:
    Matrix(ulong num_rows, ulong num_cols) {
        this->num_rows = num_rows;
        this->num_cols = num_cols;
        data.resize(num_rows * num_cols);
    }

    inline T get(ulong i, ulong j) {
        return data[i * num_cols + j];
    }

    inline void set(ulong i, ulong j, T value) {
        data[i * num_cols + j] = value;
    }
};

template<class DT>
void cluster(
        DT *array,
        ulong n,
        ulong k,
        ulong *clusters,
        DT *centroids) {
    // ***************************************************
    // * Sort input array and save info for de-sorting
    // ***************************************************

    vector<ulong> sort_idxs(n);
    iota(sort_idxs.begin(), sort_idxs.end(), 0);
    sort(
            sort_idxs.begin(),
            sort_idxs.end(),
            [&array](ulong a, ulong b) { return array[a] < array[b]; });
    vector<ulong> undo_sort_lookup(n);
    vector<DT> sorted_array(n);
    for (ulong i = 0; i < n; ++i) {
        sorted_array[i] = array[sort_idxs[i]];
        undo_sort_lookup[sort_idxs[i]] = i;
    }

    // ***************************************************
    // * Set D and T using dynamic programming algorithm
    // ***************************************************

    // Algorithm as presented in section 2.2 of (Gronlund et al., 2017).

    CostCalculator cost_calculator(sorted_array, n);
    Matrix<DT> D(k, n);
    Matrix<ulong> T(k, n);

    for (ulong i = 0; i < n; ++i) {
        D.set(0, i, cost_calculator.calc(0, i));
        T.set(0, i, 0);
    }

    for (ulong k_ = 1; k_ < k; ++k_) {
        auto C = [&D, &k_, &cost_calculator](ulong i, ulong j) -> DT {
            ulong col = i < j - 1 ? i : j - 1;
            return D.get(k_ - 1, col) + cost_calculator.calc(j, i);
        };
        vector<ulong> row_argmins = smawk<DT>(n, n, C);
        for (ulong i = 0; i < row_argmins.size(); ++i) {
            ulong argmin = row_argmins[i];
            DT min = C(i, argmin);
            D.set(k_, i, min);
            T.set(k_, i, argmin);
        }
    }

    double ratio_avg = 0;
    bool findk = false;
    for (ulong k_ = 3; k_ < k; ++k_) {
        float ratio = D.get(k_ - 1, n - 1) / D.get(k_, n - 1);
        ratio_avg = (ratio_avg * (k_ - 3) + ratio) / (k_ - 2);
        std::cout << k_ << " , " << D.get(k_, n - 1) << " , " << ratio << " , " << ratio_avg << std::endl;
        if (ratio / ratio_avg > 1.5) {
            k = k_ + 1;
            findk = true;
            break;
        }
    }
    if (!findk) {
        return;
    }
    std::cout << "k=" << k << std::endl;

    // ***************************************************
    // * Extract cluster assignments by backtracking
    // ***************************************************

    // TODO: This step requires O(kn) memory usage due to saving the entire
    //       T matrix. However, it can be modified so that the memory usage is O(n).
    //       D and T would not need to be retained in full (D already doesn't need
    //       to be fully retained, although it currently is).
    //       Details are in section 3 of (Grønlund et al., 2017).

    vector<DT> sorted_clusters(n);

    ulong t = n;
    ulong k_ = k - 1;
    ulong n_ = n - 1;
    // The do/while loop was used in place of:
    //   for (k_ = k - 1; k_ >= 0; --k_)
    // to avoid wraparound of an unsigned type.
    do {
        ulong t_ = t;
        t = T.get(k_, n_);
        DT centroid = 0.0;
        for (ulong i = t; i < t_; ++i) {
            sorted_clusters[i] = k_;
            centroid += (sorted_array[i] - centroid) / (i - t + 1);
        }
        centroids[k_] = centroid;
        k_ -= 1;
        n_ = t - 1;
    } while (t > 0);

    // ***************************************************
    // * Order cluster assignments to match de-sorted
    // * ordering
    // ***************************************************

    for (ulong i = 0; i < n; ++i) {
        clusters[i] = sorted_clusters[undo_sort_lookup[i]];
    }
}

int main(int argc, char **argv) {
    size_t num;

    SZ::Timer timer;
    timer.start();

    auto input = SZ::readfile<float>("/Users/kzhao/data/exaalt/exaalt_nano-230434x146/x.f32.dat", num);
//    auto input = SZ::readfile<float>("/Users/kzhao/data/exaalt/exaalt_trinity111-208745427/x.f32.dat", num);
//    auto input = SZ::readfile<float>("/Users/kzhao/data/exaalt/exaalt_2v18he100n1-5759x1040/x.f32.dat", num);
    timer.stop("read file");

    timer.start();
    int sample_rate = 100000;
    std::vector<float> sample;
    sample.reserve(num / sample_rate);
    std::sample(input.get(), input.get() + num,
                std::back_inserter(sample),
                num / sample_rate,
                std::mt19937{std::random_device{}()});
    timer.stop("random sample");

    timer.start();
    int k = 100;
    std::vector<ulong> idx(num);
    std::vector<float> cents(k);
    cluster(sample.data(), num / sample_rate, k, idx.data(), cents.data());
//    cluster(input.get(), num, 16, idx.data(), cents.data());
    timer.stop("kmeans1d");

    for (auto &cent:cents) {
        std::cout << cent << " ";
    }

}