#include "./optimizer/CCJaDE.h"
#include <iostream>

void ResProcessing(double res[], int func_id);

int main() {
	int dim = 1000;
	double res[30];
	for (int fid = 1; fid <= 3; fid++) {
		if (fid >= 8 && fid <= 14)
			dim = 905;
		else
			dim = 1000;
		JaDE optimizer(dim, fid);
		for (int t = 0; t < 30; t++) {
			res[t] = optimizer.Optimize(t);
		}
		ResProcessing(res, fid);
	}
	for (int fid = 8; fid <= 15; fid++) {
		if (fid >= 8 && fid <= 14)
			dim = 905;
		else
			dim = 1000;
		JaDE optimizer(dim, fid);
		for (int t = 0; t < 30; t++) {
			res[t] = optimizer.Optimize(t);
		}
		ResProcessing(res, fid);
	}
}
void ResProcessing(double res[], int func_id) {
	ofstream outRes;		//	out res in each 1000 gens
	stringstream ss;
	ss << "RES-cec2013/F" << func_id << "/sum.txt";
	outRes.open(ss.str());
	ss.clear();
	ss.str("");
	if (!outRes.is_open()) {
		cout << "can not open result file!" << endl;
	}
	double ave = 0, std = 0;
	double best = 0;
	double worst = 1e30;
	for (int t = 0; t < 30; t++) {
		outRes << "T" << t << ":\t" << res[t] << endl;
		ave += res[t];
		if (best < res[t]) best = res[t];
		if (worst > res[t]) worst = res[t];
	}
	ave = ave / 30.0;
	for (int t = 0; t < 30; t++) {
		std += (res[t] - ave) * (res[t] - ave);
	}
	std = sqrt(std / 30.0);
	outRes << "Best:\t" << best;
	outRes << "Worst:\t" << worst;
	outRes << "Mean:\t" << ave;
	outRes << "Std:\t" << std;
	outRes.close();
}
