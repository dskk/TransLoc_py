#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
using namespace std;

pair<vector<double>, vector<double>> RTA_L2(vector<double> eig_vals, vector<vector<double>> eig_vecs, vector<vector<double>> HX_x, vector<vector<double>> HX_y, vector<pair<int, int>> inds, double inv_tau) {
	int nmol = eig_vals.size();
	vector<double> x(nmol), y(nmol);
	for (int i = 0; i < nmol; i++) {
		for (int j = 0; j < nmol; j++) {
			double temp_x = 0, temp_y = 0;
			for (pair<int, int> p : inds) {
				int i1 = p.first;
				int i2 = p.second;
				temp_x += (eig_vecs[i1][i] * eig_vecs[i2][j] - eig_vecs[i1][j] * eig_vecs[i2][i]) * HX_x[i1][i2];
				temp_y += (eig_vecs[i1][i] * eig_vecs[i2][j] - eig_vecs[i1][j] * eig_vecs[i2][i]) * HX_y[i1][i2];
			}
			x[i] += temp_x * temp_x * 2 / (inv_tau * inv_tau + (eig_vals[i] - eig_vals[j]) * (eig_vals[i] - eig_vals[j]));
			y[i] += temp_y * temp_y * 2 / (inv_tau * inv_tau + (eig_vals[i] - eig_vals[j]) * (eig_vals[i] - eig_vals[j]));
		}
	}
	return make_pair(x, y);
}

namespace py = pybind11;

PYBIND11_MODULE(RTA_L2, m) {
	m.def("RTA_L2", &RTA_L2);
}