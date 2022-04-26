/*This version deal with overlapping function transformed from partially separable function*/
/*                      2018-1-25 by Jason SYSU                                            */
/*******************************************************************************************/
#include "CCJaDE.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>

using namespace std;

double JaDE::random_double(double low, double high) {
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}
int JaDE::random_int(int low, int high) {
	if (high>low)
		return low + rand() % (high - low);
	else
		return low;
}

double JaDE::random_cauchy(double alpha, double beta) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::cauchy_distribution<double> distribution(alpha, beta);
	return distribution(generator);
}
// mu: mean, sigma: standard deviation
double JaDE::random_normal(double mu = 0.0, double sigma = 1.0) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(mu, sigma);
	return distribution(generator);
}
void JaDE::shuffle_array(int ary[], int N) {
	for (int i = N - 1; i >= 0; i--)
		swap(ary[i], ary[rand() % (i + 1)]);
}
void JaDE::shuffle_vector(vector <int> &array) {
	for (int i = array.size() - 1; i >= 0; i--)
		swap(array[i], array[rand() % (i + 1)]);
}
void JaDE::OverlapTransform(double *x, double *y) {
	int m = 5;
	int od = 905;
	int xIndex = 0;
	for (int i = 0; i < all_gp[0].size(); i++) {
		int d = all_gp[0][i];
		y[d] = x[i];
	}
	xIndex += all_gp[0].size() - m;
	for (int i = 1; i < all_gp.size(); i++) {
		for (int j = 0; j < all_gp[i].size(); j++) {
			int d = all_gp[i][j];
			y[d] = x[xIndex + j];
		}
		xIndex += all_gp[i].size() - m;
	}
}
JaDE::~JaDE() {
	for (int i = 0; i < P.size(); i++) {
		P[i].CR.clear();
		P[i].F.clear();
		delete[] P[i].x;
	}
	P.clear();
	for (int i = 0; i < all_gp.size(); i++) all_gp[i].clear();
	all_gp.clear();
	for (int i = 0; i < new_agp.size(); i++) new_agp[i].clear();
	new_agp.clear();
	sep_gp.clear();
	gg.clear();
	
	for (int i = 0; i < arcdFit.size(); i++) arcdFit[i].clear();
	arcdFit.clear();
	for (int i = 0; i < thetaMat.size(); i++) thetaMat[i].clear();
	thetaMat.clear();
	for (int i = 0; i < interMat.size(); i++) interMat[i].clear();
	interMat.clear();

	xhigh.clear();
	xlow.clear();
}
JaDE::JaDE(int dim, int fid) {				// dimension and function id
	p = 0.1;
	c = 0.05;
	fp = NULL;
	func_id = fid;
	adapFlag = false;
	maxgpsize = 0;
	dimension = dim;
	samplingbuf = 200;
	doublemax = numeric_limits<double>::max();
	fp = generateFuncObj(func_id);

	P.resize(NP);
	xlow.resize(dimension);
	xhigh.resize(dimension);
	for (int i = 0; i < P.size(); i++) {
		P[i].x = new double[dimension];
		P[i].f = doublemax;
	}
	PB.x = new double[dimension];
	PB.f = doublemax;
	
	for (int i = 0; i < dimension; i++) {
		Dindex.push_back(i);
	}	
	thetaMat.resize(dimension);
	interMat.resize(dimension);
	for (int i = 0; i < dimension; i++) {
		thetaMat[i].resize(dimension);
	}
	cdim.resize(dimension, 0);
}
// random
void JaDE::Setup_group() {
	Rand_group();
	groupsize = new_agp.size();
	gg.resize(groupsize, 0);
	arcdFit.resize(dimension);				
	for (int n = 0; n < P.size(); n++) {
		P[n].F.resize(groupsize);
		P[n].CR.resize(groupsize);
	}

	uF.resize(groupsize, 0.5);
	uCR.resize(groupsize, 0.5);
}
void JaDE::Rand_group() {
	for (int i = 0; i < new_agp.size(); i++) new_agp[i].clear();
	new_agp.clear();
	shuffle_vector(Dindex);
	
	int subgpDim = 50;		//	the number of sub dimension of sub-group
	int subgpNum = dimension / subgpDim;
	for (int i = 0; i < subgpNum; i++) {
		group tmpgp;
		tmpgp.clear();
		int remaindim = dimension - i * subgpDim;
		if (remaindim >= subgpDim && remaindim < 2 * subgpDim)
			tmpgp.resize(remaindim);
		else
			tmpgp.resize(subgpDim);

		for (int d = 0; d < subgpDim; d++) {
			tmpgp[d] = Dindex[i * subgpDim + d];
		}
		new_agp.push_back(tmpgp);
	}
}
void JaDE::Dygroup() {
	Compare compare;
	vector <SortArr> sa;	
	double_archive tmpdFit;

	sa.resize(dimension);
	tmpdFit.resize(dimension, 0);
	for (int d = 0; d < dimension; d++) {
		if (arcdFit[d].size() > 0) {
			for (int i = 0; i < arcdFit[d].size(); i++) {
				tmpdFit[d] += arcdFit[d][i];
			}
			tmpdFit[d] = tmpdFit[d] / arcdFit[d].size();
		}
		else {
			tmpdFit[d] = 0;
		}
	}
	group dimIndex;
	for (int i = 0; i < dimension; i++) {
		dimIndex.push_back(i);
	}
	// shuffle dimension index before sorting
	shuffle_vector(dimIndex);
	for (int d = 0; d < dimension; d++) {
		sa[d].index = dimIndex[d];
		sa[d].val = tmpdFit[dimIndex[d]];
	}
	sort(sa.begin(), sa.end(), compare);
	new_agp[0].clear();
	int size = 0.02 * dimension;
	maxgpsize = 200;
	vector <int> hashArr;
	hashArr.resize(dimension, -1);
	int maxloop = 0;										//in case of endless loop
	while (new_agp[0].size() < maxgpsize && maxloop < 20) {
		for (int i = 0; i < size; i++) {
			int d = sa[i].index;
			if (hashArr[d] == -1) {
				new_agp[0].push_back(d);
				hashArr[d] = d;
			}
			if (new_agp[0].size() < maxgpsize) {
				for (int j = 0; j < interMat[d].size(); j++) {
					int innerDim = interMat[d][j];
					if (new_agp[0].size() < maxgpsize && hashArr[innerDim] == -1) {
						new_agp[0].push_back(innerDim);
						hashArr[innerDim] = innerDim;
						break;
					}
				}
			}
		}
		maxloop++;
	}
	hashArr.clear();
	sa.clear();
	tmpdFit.clear();
}
void JaDE::Input_group() {
	Setup_group();	// only for random grouping 
					// input theta matrix
	if (func_id >= 8 && func_id <= 11) {
		ifstream infile;
		stringstream ss;
		ss << "GP_2013/F" << func_id << ".txt";
		infile.open(ss.str());
		if (!infile.is_open()) {
			cout << "Open group results file error!" << endl;
			return;
		}
		ss.clear();
		ss.str("");

		string str;
		int agp_size, sep_size;
		infile >> str >> agp_size;
		for (int i = 0; i < agp_size; i++) {
			int gpsize;
			infile >> str >> gpsize;
			if (gpsize > 0) {
				group tmp_ngp;
				tmp_ngp.resize(gpsize);
				for (int j = 0; j < gpsize; j++) {
					infile >> tmp_ngp[j];
				}
				all_gp.push_back(tmp_ngp);
				tmp_ngp.clear();
			}
		}
		infile >> str >> sep_size;
		if (sep_size > 0) {
			sep_gp.resize(sep_size);
			for (int i = 0; i < sep_size; i++)
				infile >> sep_gp[i];
		}
		infile.close();
	}
	ifstream intheta;
	stringstream ifss;
	ifss << "Theta/theta" << func_id << ".txt";
	intheta.open(ifss.str());
	if (intheta.is_open()) {
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				intheta >> thetaMat[i][j];
				if (thetaMat[i][j] == 1 && (i != j)) {
					interMat[i].push_back(j);
				}
			}
			if (maxgpsize < interMat[i].size()) {
				maxgpsize = interMat[i].size();
			}
		}
	}
	ifss.str("");
//	cout << "Setup group" << endl;	
}
void JaDE::Evaluation(double *x, double &f) {
	if (func_id >= 12 || func_id <= 3)
		f = fp->compute(x);
	else {
		double *y = new double[1000];
		OverlapTransform(x, y);
		f = fp->compute(y);
		delete[] y;
	}
}
void JaDE::Initializ() {
	srand(time(0));
	PB.f = doublemax;	
	for (int d = 0; d < dimension; d++) {
		xlow[d] = fp->getMinX();
		xhigh[d] = fp->getMaxX();
	}
	for (int i = 0; i < P.size(); i++) {
		for (int d = 0; d < dimension; d++)
			P[i].x[d] = random_double(xlow[d], xhigh[d]);
		Evaluation(P[i].x, P[i].f);
		if (PB.f > P[i].f) {
			PB.f = P[i].f;
			PB.best_index = i;
		}
	}
	for (int d = 0; d < dimension; d++) {
		PB.x[d] = P[PB.best_index].x[d];
	}

	PB.f = P[PB.best_index].f;
	PrePB.f = PB.f;
}

double JaDE::Optimize(int timeIndex) {
	srand(time(0));
	Input_group();
	Initializ();

	/*ofstream out_res;
	stringstream ss;
	ss << "R2013/R" << func_id << ".txt";
	out_res.open(ss.str());
	ss.clear();
	ss.str("");*/

	ofstream outRes;		//	out res in each 1000 gens
	stringstream ss;
	ss << "RES-cec2013/F" << func_id  << "/T" << timeIndex << ".txt";
	outRes.open(ss.str());
	ss.clear();
	ss.str("");

	if (!outRes.is_open()) {
		cout << "can not open result file!" << endl;
	}
	int *index = new int[NP];
	double *cbx = new double[dimension];
	double_archive SF, SCR;
	bool success_flag = false;
	g = 0;
	while (g <= 60000) {
		if (g < 800) {
			cgIndex = g % groupsize;
			if (g <= 700)
				cdim[cgIndex] = 0;
			if (g % groupsize == 0 && g != 0)
				Rand_group();
		}
		else {
			if (g % 100 == 0) {
				Dygroup();
				adapFlag = true;
				uF.clear();					//	reset parameters for new dynamic group
				uCR.clear();
				uF.resize(groupsize, 0.5);
				uCR.resize(groupsize, 0.5);
				A.clear();
			}
			else {
				adapFlag = false;
			}
			cgIndex = 0;					//	new_agp[0] is the dynamic group after 800th generation
		}
		SF.clear();
		SCR.clear();
		success_flag = false;
		for (int n = 0; n < NP; n++) {
			for (int n = 0; n < NP; n++) {
				index[n] = n;
			}
			for (int d = 0; d < dimension; d++) {
				cbx[d] = PB.x[d];
			}

			double *v = new double[dimension];
			double *u = new double[dimension];

			//	choose F
			do {
				P[n].F[cgIndex] = random_cauchy(uF[cgIndex], 0.1);
			} while (P[n].F[cgIndex] < 0);
			if (P[n].F[cgIndex] > 1) {
				P[n].F[cgIndex] = 1;
			}

			//	choose CR
			P[n].CR[cgIndex] = random_normal(uCR[cgIndex], 0.1);
			if (P[n].CR[cgIndex] > 1) P[n].CR[cgIndex] = 1;
			if (P[n].CR[cgIndex] < 0) P[n].CR[cgIndex] = 0;

			//	mutation 
			int idx, idxpbest, r1, r2;	//	idx: record the index, idxpbest: p% best of P, r1: randomly selected from P, r2: randomly selected from {A, P}
			vector <SortArr> sapop;
			sapop.resize(NP);
			for (int np = 0; np < NP; np++) {
				sapop[np].index = np;
				sapop[np].val = P[np].f;
			}
			Compare cp;
			sort(sapop.begin(), sapop.end(), cp);
			idx = random_int(0, random_double(0, 0.2) * NP);
			idxpbest = sapop[NP - 1 - idx].index;
			sapop.clear();
			swap(index[idxpbest], index[NP - 1]);
			idx = random_int(0, NP - 1);
			r1 = index[idx];
			swap(index[idx], index[NP - 2]);
			int temp_PAsize = P.size() - 2 + A.size();
			r2 = random_int(0, temp_PAsize);

			for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
				int d = new_agp[cgIndex][innd];
				if (r2 < (P.size() - 2))
					v[d] = P[n].x[d] + P[n].F[cgIndex] * (P[idxpbest].x[d] - P[n].x[d]) + P[n].F[cgIndex] * (P[r1].x[d] - P[r2].x[d]);
				else
					v[d] = P[n].x[d] + P[n].F[cgIndex] * (P[idxpbest].x[d] - P[n].x[d]) + P[n].F[cgIndex] * (P[r1].x[d] - A[r2 - (P.size() - 2)].x[d]);
			}

			//	crossover
			int jrand_idx = random_int(0, new_agp[cgIndex].size());
			int jrand = new_agp[cgIndex][jrand_idx];
			for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
				int d = new_agp[cgIndex][innd];
				if (random_double(0, 1) < P[n].CR[cgIndex])
					u[d] = v[d];
				else
					u[d] = P[n].x[d];
			}
			u[jrand] = v[jrand];
			for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
				int d = new_agp[cgIndex][innd];
				if (u[d] < xlow[d]) u[d] = xlow[d];
				if (u[d] > xhigh[d]) u[d] = xhigh[d];
			}
			for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
				int d = new_agp[cgIndex][innd];
				cbx[d] = P[n].x[d];
			}
			
			Evaluation(cbx, P[n].f);
			for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
				int d = new_agp[cgIndex][innd];
				cbx[d] = u[d];
			}

			double uf;	// f
			Evaluation(cbx, uf);

			// selection
			if (P[n].f > uf) {
				success_flag = true;
				//	archive inferior individual, update A			
				A.push_back(P[n]);		
				// update P
				P[n].f = uf;
				for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
					int d = new_agp[cgIndex][innd];
					P[n].x[d] = cbx[d];
				}
				SF.push_back(P[n].F[cgIndex]);
				SCR.push_back(P[n].CR[cgIndex]);
			}
			
			if (PB.f > P[n].f) {
				PB.f = P[n].f;
				PB.best_index = n;
				for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
					int d = new_agp[cgIndex][innd];
					PB.x[d] = P[PB.best_index].x[d];
				}
			}
			delete[] v;
			delete[] u;
		}
		// randomly remove solutions from A so that |A| < NP
		if (A.size() > NP) {
			int_archive popIndex;
			for (int i = 0; i < A.size(); i++) {
				popIndex.push_back(i);
			}
			shuffle_vector(popIndex);
			vector <SOLUTION> A2;
			for (int i = 0; i < P.size(); i++) {
				A2.push_back(A[popIndex[i]]);
			}
			A.clear();
			for (int i = 0; i < P.size(); i++) {
				A.push_back(A2[i]);
			}
			A2.clear();
			popIndex.clear();
		}
		//	update uF and uCR
		if (SF.size() > 0 && adapFlag == false) {
			double uA_SCR = 0, uL_SF = 0;
			double sum_F = 0, sum_FF = 0;
			for (int ii = 0; ii < SCR.size(); ii++) {
				uA_SCR += SCR[ii];
				sum_F += SF[ii];
				sum_FF += SF[ii] * SF[ii];
			}
			uA_SCR = uA_SCR / (SCR.size() * 1.0);
			if (sum_F > 0)
				uL_SF = sum_FF / sum_F;
			else
				uL_SF = uF[cgIndex];
			uCR[cgIndex] = (1 - c) * uCR[cgIndex] + c * uA_SCR;
			uF[cgIndex] = (1 - c) * uF[cgIndex] + c * uL_SF;
		}

		double tmpdFit = fabs(PB.f - PrePB.f);
		for (int innd = 0; innd < new_agp[cgIndex].size(); innd++) {
			int d = new_agp[cgIndex][innd];
			if (arcdFit[d].size() < samplingbuf)
				arcdFit[d].push_back(tmpdFit);
			else
				arcdFit[d][cdim[d] % samplingbuf] = tmpdFit;
		}
		for (int gpid = 0; gpid < new_agp[cgIndex].size(); gpid++) {
			int d = new_agp[cgIndex][gpid];
			cdim[d]++;
		}
		// archive PB
		gg[cgIndex]++;
		if (g % 1000 == 0) {
			cout << endl;
			outRes << endl;
			cout << "gen: " << g << "  " << PB.f << endl;
//			outRes << "gen: " << g << "  " << PB.f << endl;
			outRes << PB.f << endl;
		}
		g++;
		prIndex = cgIndex;
		PrePB.f = PB.f;
	}

	delete[] index;
	delete[] cbx;
	outRes.close();
	return PB.f;
}
