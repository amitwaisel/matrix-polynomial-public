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
#include <chrono>
#include <iostream>

#include <mkl.h>
#include "Matrix.h"
#include "Polynomial.h"

#include "mex.h"

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
	PatersonStockmeyerType ps_type = PatersonStockmeyerType::PatersonStockmeyerSequential;

	if (params.find("type") != params.end())
	{
		std::string polynomial_type = params["type"]->string_data();
		if (polynomial_type == "sequential")
			ps_type = PatersonStockmeyerType::PatersonStockmeyerSequential;
		else if (polynomial_type == "parallel")
			ps_type = PatersonStockmeyerType::PatersonStockmeyerParallel;
		else if (polynomial_type == "vanloan")
			ps_type = PatersonStockmeyerType::VanLoan;
		else
			mexErrMsgIdAndTxt("MATLAB:MatrixPolynomial:invalidParameterPSType", polynomial_type.c_str());
	}
	if (params.find("nworkers") != params.end())
		nworkers = (int)params["nworkers"]->numeric_data();
	if (nworkers > 0)
	{
		mkl_set_num_threads(nworkers);
#if USE_CILK
		char nworkers_ascii[10] = { 0 };
		sprintf_s(nworkers_ascii, "%d", nworkers);
		__cilkrts_set_param("nworkers", nworkers_ascii);
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
	Submatrix<ElmType> mat(input_matrix);
	Benchmark::Instance().Start();
	Matrix<ElmType> poly = mat.polynomial(co, origin, m, ps_type);
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculation);


	memcpy(qA, poly.data(), m * n * sizeof(*qA));

	if (nlhs > 1)
	{
		std::chrono::duration<double> total = Benchmark::Instance().Get(Benchmark::CheckpointType::PolynomialCalculation);
		plhs[1] = mxCreateDoubleScalar(std::chrono::duration_cast<std::chrono::milliseconds>(total).count());
	}
}
