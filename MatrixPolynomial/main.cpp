#include <cstdio>
#include <cstdlib>
#include <string>
#include <exception>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <complex>
#include <functional>
#include <chrono>  // for high_resolution_clock
#include <iostream>

//#define MKL_Complex8 std::complex<float>
#define USE_MEX
#include <mkl.h>
#include "Matrix.h"
#include "Polynomial.h"

#include "mex.h"

// http://www.cs.rochester.edu/~bh/cs400/using_lapack.html

/*
S	Real, the same as float in C
D	Double precision, the same as double in C
C	Complex
Z	Double complex
*/

/*
N		number of columns of A
nrhs	number of columns of b, usually 1
lda		number of rows (Leading Dimension of A) of A
ipiv	pivot indices
ldb		number of rows of b
*/



// CLAPACK http://icl.cs.utk.edu/lapack-for-windows/clapack/

class Argument
{
public:
	Argument(std::string name, mxClassID type, std::string data) : _name(name), _type(type), _string_data(data) { }
	Argument(std::string name, mxClassID type, double data) : _name(name), _type(type), _numeric_data(data) { }

	std::string name() { return _name; }
	mxClassID type() { return _type; }
	std::string string_data() { return _string_data; }
	double numeric_data() { return _numeric_data; }

private:
	std::string _name;
	mxClassID _type;
	std::string _string_data;
	double _numeric_data;
};

std::map<std::string, std::shared_ptr<Argument>> ParseArgs(const mxArray* params_arg)
{
	std::map<std::string, std::shared_ptr<Argument>> params;
	int nfields = mxGetNumberOfFields(params_arg);
	if (mxGetNumberOfElements(params_arg) != 1)
		mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidArgs", "Parameters structure cannot contain more than one structure");

	// mwSize ndim = mxGetNumberOfDimensions(params_arg);
	// const mwSize* dims = mxGetDimensions(params_arg);

	for (int ifield = 0; ifield < nfields; ifield++)
	{
		mxArray* tmp = mxGetFieldByNumber(params_arg, 0, ifield);
		if (tmp == NULL)
		{
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParamsField", "Parameters field is empty!");
		}
		if ((!mxIsChar(tmp) && !mxIsNumeric(tmp)) || mxIsSparse(tmp))
		{
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParamsField", "Parameters field must have either string or numeric non-sparse data.");
		}
		std::string field_name = mxGetFieldNameByNumber(params_arg, ifield);
		mxClassID class_id = mxGetClassID(tmp);

		if (mxIsChar(tmp))
		{
			std::string data(1024, '\0');
			mxGetString(tmp, &data[0], data.size());
			data.resize(strlen(data.c_str()));
			params[field_name] = std::make_shared<Argument>(field_name, class_id, data);
		}
		else if (mxIsNumeric(tmp))
		{
			params[field_name] = std::make_shared<Argument>(field_name, class_id, mxGetScalar(tmp));
		}
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParamsField", "Unsupported parameters field type");
	}
	return params;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *A, *qA; /* pointers to input & output matrices*/
	size_t m, n;      /* matrix dimensions */
						 /* form of op(A) & op(B) to use in matrix multiplication */

 	if (nrhs < 1 || nlhs < 1 || nlhs > 4)
	{
		mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidArgs", "Invalid number of input/output arguments");
	}

	//param = mxGetScalar(prhs[2]);

	A = mxGetPr(prhs[0]); /* input matrix */
	std::vector<float> co;
	int coefficients_count = nrhs - 1;
	std::map<std::string, std::shared_ptr<Argument>> params;
	if (mxIsStruct(prhs[nrhs - 1]))
	{
		coefficients_count--;
		params = ParseArgs(prhs[nrhs - 1]);
	}

	double delta = 0.1;
	int nworkers = 0;
	SchurDecompositionAlgorithm schur_algorithm = SchurDecompositionAlgorithm::Default;
	EigenvaluesClusterSortAlgorithm eigenvalues_sort_algorithm = EigenvaluesClusterSortAlgorithm::SmartSwapSort;
	PolynomialBlockAlgorithm blocks_algorithm = PolynomialBlockAlgorithm::ParlettRecurrence;
	SylvesterFunction sylvester_function = SylvesterFunction::SylvesterLapack;

	if (params.find("delta") != params.end())
		delta = params["delta"]->numeric_data();
	if (params.find("schur") != params.end())
	{
		std::string schur_algorithm_string = params["schur"]->string_data();
		if (schur_algorithm_string == "default")
			schur_algorithm = SchurDecompositionAlgorithm::Default;
		else if (schur_algorithm_string == "hassenberg")
			schur_algorithm = SchurDecompositionAlgorithm::Hassenberg;
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParameterSchurAlgorithm", schur_algorithm_string.c_str());
	}
	if (params.find("blocks_algorithm") != params.end())
	{
		std::string blocks_algorithm_string = params["blocks_algorithm"]->string_data();
		if (blocks_algorithm_string == "parlett")
			blocks_algorithm = PolynomialBlockAlgorithm::ParlettRecurrence;
		else if (blocks_algorithm_string == "parallel")
			blocks_algorithm = PolynomialBlockAlgorithm::ParallelParlettRecurrence;
		else if (blocks_algorithm_string == "sequential")
			blocks_algorithm = PolynomialBlockAlgorithm::Sequential;
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParameterBlocksAlgorithm", blocks_algorithm_string.c_str());
	}
	if (params.find("sylvester") != params.end())
	{
		std::string sylvester_string = params["sylvester"]->string_data();
		if (sylvester_string == "lapack")
			sylvester_function = SylvesterFunction::SylvesterLapack;
		else if (sylvester_string == "recursive")
			sylvester_function = SylvesterFunction::SylvesterRecusrive;
		else if (sylvester_string == "intel")
			sylvester_function = SylvesterFunction::SylvesterIntel;
		else if (sylvester_string == "recsy")
			sylvester_function = SylvesterFunction::SylvesterRecsy;
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParameterSylvester", sylvester_string.c_str());
	}
	if (params.find("eigenvalues_sort") != params.end())
	{
		std::string sort_string = params["eigenvalues_sort"]->string_data();
		if (sort_string == "bubble")
			eigenvalues_sort_algorithm = EigenvaluesClusterSortAlgorithm::BubbleSort;
		else if (sort_string == "smart")
			eigenvalues_sort_algorithm = EigenvaluesClusterSortAlgorithm::SmartSwapSort;
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParameterEigenvaluesSort", sort_string.c_str());
	}
	if (params.find("nworkers") != params.end())
		nworkers = (int)params["nworkers"]->numeric_data();
	if (nworkers > 0)
	{
		mkl_set_num_threads(nworkers);
#if USE_CILK
		char nworkers_ascii[10] = { 0 };
		sprintf_s(nworkers_ascii, "%d", nworkers);
		__cilkrts_end_cilk();
		if (__cilkrts_set_param("nworkers", nworkers_ascii) != 0)
		{
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:cilkError", "Failed setting nworkers");
		}
#endif
	}

	for (int i = 1; i <= coefficients_count; i++)
		co.push_back(mxGetScalar(prhs[i]));

	if (co.empty())
	{
		mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidArgs", "Invalid input polynomial factors - cannot be empty");
	}

	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	if (m != n)
	{
		mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidArgs", "Invalid input matrix dimensions - must be square");
	}

	/* create output matrix C */
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	qA = mxGetPr(plhs[0]);

	/* Do polynomial calculation */
#define ElmType double
	Matrix<ElmType> input_matrix(A, m, n);
	Polynomial<ElmType> p(input_matrix);
	Matrix<ElmType> poly = p.calculatePolynomial(co, delta, schur_algorithm, eigenvalues_sort_algorithm, blocks_algorithm, sylvester_function);
	memcpy(qA, poly.data(), m * n * sizeof(*qA));

	if (nlhs > 1)
	{
		double* timings = nullptr;
		double* properties = nullptr;
		std::chrono::duration<double> total;
		size_t checkpoints_count = (size_t)Benchmark::CheckpointType::Checkpoint_MAX;
		size_t properties_count = (size_t)Benchmark::PropertyType::Property_MAX;
		if (nlhs > 2)
		{
			plhs[2] = mxCreateDoubleMatrix(checkpoints_count, 3, mxREAL);
			timings = mxGetPr(plhs[2]);
		}
		if (nlhs > 3)
		{
			plhs[3] = mxCreateDoubleMatrix(properties_count, 1, mxREAL);
			properties = mxGetPr(plhs[3]);
		}
		for (size_t i = 0; i < checkpoints_count; i++)
		{
			auto time = Benchmark::Instance().Get((Benchmark::CheckpointType)i);
			auto cpu = Benchmark::Instance().GetCpu((Benchmark::CheckpointType)i);
			total += time;
			if (timings != nullptr)
			{
				timings[i] = time.count();
				timings[i + 1 * checkpoints_count] = Benchmark::Instance().GetCount((Benchmark::CheckpointType)i);
				//double time_ticks = time.count() / 100.0; // ticks are 100ns
				//double cpu_percentage = cpu / time_ticks;
				timings[i + 2 * checkpoints_count] = cpu;
			}
		}
		if (properties != nullptr)
		{
			for (size_t i = 0; i < properties_count; i++)
			{
				properties[i] = Benchmark::Instance().Get((Benchmark::PropertyType)i);
			}
		}
		plhs[1] = mxCreateDoubleScalar(std::chrono::duration_cast<std::chrono::milliseconds>(total).count());
	}
}

// Way 1: Layers + small SYL equations. Paralleled by layers
// Way 2: From left to right. Fewer SYLs. Parallelism ??

/*
benchmarking:
Automate the benchmarking for easy execution/reexecution. Generate graphs in matlab based on the exported results
1. Each scaling experiment has to have X - size of matrix (up to 8K*8K), Y - time (logarithmic scale or divided by n^3)
- sequential sizes - 1 to 8K, with 24 cores
- Show X datasets on each graph (all different ways to calculate and the naive)
- Compare with baseline of a known 3rd party algo - recsy / schur decomp.
	- OR breakdown the execution time to the different components - schur, my code, sylvester - and show the percentage of the CPU time in the polynomial calculation
2. Scaling #2 - X axis is #cores, on fixes matrix size

3 different parameters in benchmark:
- size of clusters on the diagonal (clusters with 50 values, clusters with 1K values)
- The variance of the cluster sizes (all are the same size, different sizes) (half 50-size blocks, half 1K-size blocks)
- The rank of the polynomial (rank 4, rank n-1, rank 1K)

Input generation:
Generate clustered matrix manually in matlab
generate random matrix which will contain the "eigenvectors"
L = diag(1,20)
V = randn(20,20)
A = V\(L*V)
spy(schur(A))
- make sure during the algo execution that the clusters are as expected.
 */

/* TODO:
0. Convert to real values - done.
1. Implement way#2 in addition to BlockParlettRecurrence (way#1)
2. Parallel
- matrix multiply (split to several concurrent LAPACK calls)
- Patterson-Stockmeyer (Van Loan - PS-MV)
- way#1 - BlockParlettRecurrence by layers
- way#2 - ??
3. Benchmarking - way#1 vs way#2 vs P-way#1 vs P-way#2

way#1 - parallel in layers + parallel matrix multiplication + recsy non-parallel
way#2 - parallel recsy, sequential blocks









histograms + graphs

1. 1-1 blocks - eigenvalues very different. schur-parlett on the whole matrix. Work of n^3.
2. clusters of 100-100 out of 8000-8000 matrix with high rank (1000). SYL on bigger chunks.
	Try to prove that higham-davis algo from sec. 9 gives a boost to the performance, where schur-parlett can't be run on (very similar eigenvalues, high error rate)
3. Few clusters (naive polynomial calculation)

Things we want to prove:
1. small blocks give perf. boost over the naive approach (A^300 will be calculated in sqrt(d) steps, not so naive)
2. diff between way#1 and way#2.


*/
