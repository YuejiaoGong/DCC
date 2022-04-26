#ifndef _CCJADE_H_
#define _CCJADE_H_

#include <vector>
#include <string>
#include <fstream>
#include "../cec2013/Header.h"

using namespace std;

typedef vector <int> int_archive;
typedef vector <double> double_archive;

typedef vector <int> group;
typedef vector <group> agroup;

class JaDE {
public:
	struct SOLUTION {
		double *x;								//	dim = dimension					
		double f;								//	f1
		double_archive CR;						//	Crossover
		double_archive F;						//	Mutation	
	};

	struct PBEST {
		int		best_index;
		double	f;
		double	*x;
	};

	struct SortArr {
		double val;
		int index;
	};

	class Compare {
	public:
		bool operator() (const SortArr &x, const SortArr &y) { return x.val > y.val; }	//	descending
	};
	
	JaDE(int dim, int fid);
	~JaDE();
	double Optimize(int timeIndex);

private:
	int g;
	int_archive gg;								//	generation of each group being optimized

	Benchmarks* fp;
	double	doublemax;
	int		func_id;
	int		dimension;

	// CC
	int		cgIndex;							//	Current group for optimization in CC, roundbin mode {}
	int		prIndex;							//	Next group for optimization
	int		groupsize;							//	groupsize = new_agp.size() 	
	group	sep_gp;								//	original group 
	agroup	all_gp;						
	agroup	new_agp;							//	dynamic group

	//	GBCC
	bool	adapFlag;							//	default value is false. set to true when new group is setup. Adaptation flag
	int		samplingbuf;
	int		maxgpsize;							//	maximum group size, maxgpsize = arg max (intermat[i])
	agroup	thetaMat;							//	theta matrix
	agroup	interMat;							//	interMat[i]: the dimensions interact with dimension i
	int_archive Dindex;
	int_archive cdim;
	vector <double_archive> arcdFit;			//	Archive deltafit, arcdFit[groupsize][10], each arcdFit[i] is a queue
												//	SaNSDE
	PBEST	PB;
	PBEST	PrePB;								//	the last PB
	vector <SOLUTION> P;						//	population
	vector <SOLUTION> A;						//	Archive population

	static const int	NP = 50;
	double_archive		uCR, uF;				//	size: groupsize
	double_archive		xlow, xhigh;
	double p, c;

	//	Function of SaNSDE
	void Initializ();
	void Evaluation(double *x, double &f); 
	void Setup_group();
	void Input_group();
	void Rand_group();

	//	Random function
	void	shuffle_array(int ary[], int N);
	void	shuffle_vector(std::vector <int> &array);
	int		random_int(int low, int high);
	double	random_double(double low, double high);
	double	random_cauchy(double alpha, double beta);
	double	random_normal(double mu, double sigma);

	//	Auxiliary function
	void Dygroup();
	void OverlapTransform(double *x, double *y);//	transform non-overlapped diensions into overlapped dimensions
};

#endif

