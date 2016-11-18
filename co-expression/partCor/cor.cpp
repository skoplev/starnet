#include <Rcpp.h>

using namespace Rcpp;

struct Interaction {
	int from;
	int to;
	float weight;
};

// Calcualtes partial correlation matrices for more efficient
// memory usage of large.
// [[Rcpp::export]]
DataFrame partCor(NumericMatrix dmat, float min_cor) {

	// Aligned data structures for each identified interaction
	IntegerVector from;
	IntegerVector to;
	NumericVector weight;

	// Allocation of running sums calculacted for each correlation evaluated
	float sum_x = 0.0;
	float sum_y = 0.0;
	float sum_x2 = 0.0;
	float sum_y2 = 0.0;
	float sum_xy = 0.0;

	// Loop over each combinatorial comparison of rows and columns in entry in matrix
	for (int i = 0; i < dmat.ncol(); i++) {
		for (int j = i + 1; j < dmat.ncol(); j++) {
			// Loop over each sample

			// Reset running sums
			sum_x = 0.0;
			sum_y = 0.0;
			sum_x2 = 0.0;
			sum_y2 = 0.0;
			sum_xy = 0.0;

			for (int k = 0; k < dmat.nrow(); k++) {
				// Get data
				float x = dmat(k, i);
				float y = dmat(k, j);

				// Update running sums
				sum_x += x;
				sum_y += y;
				sum_x2 += pow(x, 2);
				sum_y2 += pow(y, 2);
				sum_xy += x*y;
			}

			// Calculate correlation for features i and j.
			float n = dmat.nrow();  // num of samples
			float r = (n * sum_xy - sum_x * sum_y) / (sqrt(n * sum_x2 - pow(sum_x, 2)) * sqrt(n * sum_y2 - pow(sum_y, 2)));

			// Add if correlation is higher than threshold
			if (std::abs(r) > min_cor) {
				from.push_back(i + 1);  // 0-index to 1-index
				to.push_back(j + 1);
				weight.push_back(r);
			}
		}
	}

	return DataFrame::create(
		Named("from")=from,
		Named("to")=to,
		Named("weight")=weight
	);
}