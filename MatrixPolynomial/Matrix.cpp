#include <mkl.h>
#include "Matrix.h"
#include "mkl_scalapack.h"
#include <algorithm>

template<typename ElementType>
Matrix<ElementType>::Matrix(int nrows, int ncols) : _data(nrows * ncols), _nrows(nrows), _ncols(ncols)
{
	memset(_data.data(), 0, _data.size() * sizeof(ElementType));
}

template<typename ElementType>
Matrix<ElementType>::Matrix(const ElementType* data, int nrows, int ncols) : Matrix(nrows, ncols)
{
	memcpy(_data.data(), data, _data.size() * sizeof(ElementType));
}

template<typename ElementType>
Matrix<ElementType>::Matrix(const Matrix<ElementType>& other) : Matrix(other.nrows(), other.ncols())
{
	reset(other);
}

template<typename ElementType>
Matrix<ElementType>::Matrix(Matrix<ElementType>&& other) : _data(other._data), _nrows(other._nrows), _ncols(other._ncols)
{
	other._data = std::vector<ElementType>();
}

template<typename ElementType>
Matrix<ElementType>& Matrix<ElementType>::operator=(Matrix&& other)
{
	std::swap(_data, other._data);
	_nrows = other._nrows;
	_ncols = other._ncols;
	other._data = std::vector<ElementType>();
	return *this;
}

template<typename ElementType>
Matrix<ElementType>::~Matrix()
{
}

template<typename ElementType>
void Matrix<ElementType>::reset()
{
	memset(_data.data(), 0, _data.size() * sizeof(ElementType));
}

template<typename ElementType>
void Matrix<ElementType>::reset(const Matrix& other)
{
	memcpy(_data.data(), other._data.data(), other._data.size() * sizeof(ElementType));
}

template<typename ElementType>
void Matrix<ElementType>::add(const Matrix& M, float scalar, MatrixIndex from)
{

	ElementType alpha = { scalar };
	if (this->nrows() < from.first + M.nrows() || this->ncols() < from.second + M.ncols())
		throw std::runtime_error("Dimensions mismatch in Matrix::add");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	cilk_for(int row = 0; row < M.nrows(); row++)
	{
#if USE_CILK
		this->data(from.first + row, from.second)[0:M.ncols()] += M.data(row, 0)[0:M.ncols()] * scalar;
#else
		ElementType* start = (ElementType*)this->data(from.first + row, from.second);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.ncols(), alpha, M.data(row, 0), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.ncols(), &alpha, M.data(row, 0), 1, start, 1);
#endif
#endif

#else
	cilk_for(int col = 0; col < M.ncols(); col++)
	{
#if USE_CILK
		this->data(from.first, from.second + col)[0:M.nrows()] += M.data(0, col)[0:M.nrows()] * scalar;
#else
		ElementType* start = (ElementType*)this->data(from.first, from.second + col);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.nrows(), alpha, M.data(0, col), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.nrows(), &alpha, M.data(0, col), 1, start, 1);
#endif
#endif
#endif
	}
}
template<typename ElementType>
void Submatrix<ElementType>::add(const Submatrix& M, float scalar, MatrixIndex from)
{
	ElementType alpha = { scalar };
	if (this->nrows() < from.first + M.nrows() || this->ncols() < from.second + M.ncols())
		throw std::runtime_error("Dimensions mismatch in Matrix::add");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	cilk_for (int row = 0; row < M.nrows(); row++)
	{
#if USE_CILK
		this->data(from.first + row, from.second)[0:M.ncols()] += M.data(row, 0)[0:M.ncols()] * scalar;
#else
		ElementType* start = (ElementType*)this->data(from.first + row, from.second);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.ncols(), alpha, M.data(row, 0), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.ncols(), &alpha, M.data(row, 0), 1, start, 1);
#endif
#endif

#else
	cilk_for(int col = 0; col < M.ncols(); col++)
	{
#if USE_CILK
		this->data(from.first, from.second + col)[0:M.nrows()] += M.data(0, col)[0:M.nrows()] * scalar;
#else
		ElementType* start = (ElementType*)this->data(from.first, from.second + col);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.nrows(), alpha, M.data(0, col), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.nrows(), &alpha, M.data(0, col), 1, start, 1);
#endif
#endif
#endif
	}
}

template<typename ElementType>
void Matrix<ElementType>::assign(const Matrix& M, MatrixIndex from)
{
	int row = from.first;
	int col = from.second;
	if (col + M.ncols() > this->ncols() || row + M.nrows() > this->nrows())
		throw std::runtime_error("Cannot assign matrix, out of bounds");
	
	int res = LAPACKE_FUNC(ElementType, lacpy, MATRIX_ORDER, 'A', M.nrows(), M.ncols(), M.data(), M.dim(), this->data(row, col), dim());
	if (res != 0)
		throw std::runtime_error("Failed LAPACKE_clacpy");
	//clacpy("A", &m_rows, &m_cols, M.data(), &m_dim, this->data(row, col), &this_dim);
}

template<typename ElementType>
void Matrix<ElementType>::multiply(const Matrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB)
{
	Matrix<ElementType> temp = multiply(*this, B, m, n, k, fromA, fromB);
	assign(temp, fromA);
	//int from_temp = 0;
	//pcgemr2d(&temp_dim, &temp_dim, temp.data(), &from_temp, &from_temp,)
}

template<typename ElementType>
void Matrix<ElementType>::multiply(const Matrix& B)
{
	if (this->ncols() != B.nrows())
		throw std::runtime_error("Dimensions mismatch");
	return this->multiply(B, this->nrows(), this->ncols(), B.ncols(), { 0, 0 }, {0, 0});
}

template<typename ElementType>
Matrix<ElementType> Matrix<ElementType>::multiply(const Matrix& A, const Matrix& B)
{
	if (A.ncols() != B.nrows())
		throw std::runtime_error("Dimensions mismatch");
	return multiply(A, B, A.nrows(), B.ncols(), A.ncols(), MatrixIndex(0, 0), MatrixIndex(0, 0));
}

template<typename ElementType>
Matrix<ElementType> Matrix<ElementType>::multiply(const Matrix& A, const Matrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB)
{
	ElementType alpha = { 1 };
	ElementType beta = { 0 };

	if (fromA.first + m > A._nrows || fromA.second + k > A._ncols)
		throw std::runtime_error("A submatrix out of bounds");
	if (fromB.first + k > B._nrows || fromB.second + n > B._ncols)
		throw std::runtime_error("B submatrix out of bounds");


#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	const ElementType* Ap_start = A._data.data() + A.dim() * fromA.first + fromA.second;
	const ElementType* Bp_start = B._data.data() + B.dim() * fromB.first + fromB.second;
#else
	const ElementType* Ap_start = A._data.data() + fromA.first + A.dim() * fromA.second;
	const ElementType* Bp_start = B._data.data() + fromB.first + B.dim() * fromB.second;
#endif
	Matrix<ElementType> temp(m, n);
	ElementType* temp_start = temp._data.data();
#if CALCULATE_REAL
	CBLAS_FUNC(ElementType, gemm, CBLAS_MATRIX_ORDER, CblasNoTrans, CblasNoTrans, m, n, k, alpha, Ap_start, A.dim(), Bp_start, B.dim(), beta, temp_start, temp.dim());
#endif
#if CALCULATE_COMPLEX
	CBLAS_FUNC(ElementType, gemm3m, CBLAS_MATRIX_ORDER, CblasNoTrans, CblasNoTrans, m, n, k, &alpha, Ap_start, A.dim(), Bp_start, B.dim(), &beta, temp_start, temp.dim());
#endif
	return temp;
}

template<typename ElementType>
void Submatrix<ElementType>::initialize_I()
{
	memset(data(), 0, n*n);
#if USE_CILK
	data()[0:n: n + 1] = 1;
#else
	cilk_for(int i = 0; i < n; i++)
	{
		*(data(i, i)) = { 1 };
	}
#endif
}

//// [in] A:m*k, [in] B:k*n, [ret] C:m*n [C=A*B]
//template<typename ElementType>
//Matrix<ElementType> Submatrix<ElementType>::multiply(const Submatrix& A, const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, bool transA, bool transB)
//{
//	Matrix<ElementType> temp(m, n);
//	Submatrix sTemp(temp);
//	multiply_ex(A, B, sTemp, m, n, k, fromA, fromB, origin, transA, transB);
//	return temp;
//}

template<typename ElementType>
ElementType Submatrix<ElementType>::GetMultiplier(Submatrix<ElementType>::MultiplyOption option)
	{
		switch (option)
		{
		case MultiplyOption::Plus:
			return 1;
		case MultiplyOption::Minus:
			return -1;
		case MultiplyOption::None:
		default:
			return 0;
		}
	}

// [in] A:m*k, [in] B:k*n, [out] C:m*n [C+=A*B]
template<typename ElementType>
void Submatrix<ElementType>::multiply_ex(const Submatrix& A, const Submatrix& B, Submatrix& C, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, MatrixIndex fromC, bool transA, bool transB, MultiplyOption option_AB, MultiplyOption option_C)
{
	ElementType alpha = GetMultiplier(option_AB);
	ElementType beta = GetMultiplier(option_C);

	if (fromA.first + m > A._nrows || fromA.second + k > A._ncols)
		throw std::runtime_error("A submatrix out of bounds");
	if (fromB.first + k > B._nrows || fromB.second + n > B._ncols)
		throw std::runtime_error("B submatrix out of bounds");
	if (fromC.first + m > C._nrows || fromC.second + n > C._ncols)
		throw std::runtime_error("C submatrix out of bounds");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	const ElementType* Ap_start = A._data + A._dim * fromA.first + fromA.second;
	const ElementType* Bp_start = B._data + B._dim * fromB.first + fromB.second;
	ElementType* Cp_start = C._data + C._dim * fromC.first + fromC.second;
#else
	const ElementType* Ap_start = A._data + fromA.first + A._dim * fromA.second;
	const ElementType* Bp_start = B._data + fromB.first + B._dim * fromB.second;
	ElementType* Cp_start = C._data + fromC.first + C._dim * fromC.second;
#endif
	CBLAS_TRANSPOSE opA = transA ? CblasConjTrans : CblasNoTrans;
	CBLAS_TRANSPOSE opB = transB ? CblasConjTrans : CblasNoTrans;
#if CALCULATE_REAL
	CBLAS_FUNC(ElementType, gemm, CBLAS_MATRIX_ORDER, opA, opB, m, n, k, alpha, Ap_start, A._dim, Bp_start, B._dim, beta, Cp_start, C._dim);
#endif
#if CALCULATE_COMPLEX
	CBLAS_FUNC(ElementType, gemm, CBLAS_MATRIX_ORDER, opA, opB, m, n, k, &alpha, Ap_start, A._dim, Bp_start, B._dim, &beta, Cp_start, C._dim);
#endif
}

template<typename ElementType>
void Submatrix<ElementType>::multiply3_ex(
	const Submatrix& A, const Submatrix& B, const Submatrix& C, Submatrix& D,
	int m, int n, int k, int l,
	MatrixIndex fromA, MatrixIndex fromB, MatrixIndex fromC, MatrixIndex fromD,
	bool transA, bool transB, bool transC,
	MultiplyOption option_ABC, MultiplyOption option_D,
	std::shared_ptr<Submatrix<ElementType>> temp)
{
	ElementType alpha = GetMultiplier(option_ABC);
	ElementType one = GetMultiplier(Plus);
	ElementType beta = GetMultiplier(option_D);

	if (fromA.first + m > A._nrows || fromA.second + k > A._ncols)
		throw std::runtime_error("A submatrix out of bounds");
	if (fromB.first + k > B._nrows || fromB.second + l > B._ncols)
		throw std::runtime_error("B submatrix out of bounds");
	if (fromC.first + l > C._nrows || fromC.second + n > C._ncols)
		throw std::runtime_error("C submatrix out of bounds");
	if (fromD.first + m > D._nrows || fromD.second + n > D._ncols)
		throw std::runtime_error("D submatrix out of bounds");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	const ElementType* Ap_start = A._data + A._dim * fromA.first + fromA.second;
	const ElementType* Bp_start = B._data + B._dim * fromB.first + fromB.second;
	const ElementType* Cp_start = C._data + C._dim * fromC.first + fromC.second;
	ElementType* Dp_start = D._data + D._dim * fromD.first + fromD.second;
#else
	const ElementType* Ap_start = A._data + fromA.first + A._dim * fromA.second;
	const ElementType* Bp_start = B._data + fromB.first + B._dim * fromB.second;
	const ElementType* Cp_start = C._data + fromC.first + C._dim * fromC.second;
	ElementType* Dp_start = D._data + fromD.first + D._dim * fromD.second;
#endif
	CBLAS_TRANSPOSE opA = transA ? CblasConjTrans : CblasNoTrans;
	CBLAS_TRANSPOSE opB = transB ? CblasConjTrans : CblasNoTrans;
	CBLAS_TRANSPOSE opC = transC ? CblasConjTrans : CblasNoTrans;
	std::shared_ptr<Matrix<ElementType>> temp_matrix;
	if (temp == nullptr)
	{
		temp_matrix = std::make_shared<Matrix<ElementType>>(m, l);
		temp = std::make_shared<Submatrix<ElementType>>(*temp_matrix);
	}
#if CALCULATE_REAL
	CBLAS_FUNC(ElementType, gemm, CBLAS_MATRIX_ORDER, opA, opB, m, l, k, alpha, Ap_start, A._dim, Bp_start, B._dim, one, temp->data(), temp->dim());
	CBLAS_FUNC(ElementType, gemm, CBLAS_MATRIX_ORDER, CblasNoTrans, opC, m, n, l, one, temp->data(), temp->dim(), Cp_start, C._dim, beta, Dp_start, D._dim);
#endif
#if CALCULATE_COMPLEX
	CBLAS_FUNC(ElementType, gemm3m, CBLAS_MATRIX_ORDER, opA, opB, m, l, k, &alpha, Ap_start, A._dim, Bp_start, B._dim, &one, temp->data(), temp->dim());
	CBLAS_FUNC(ElementType, gemm3m, CBLAS_MATRIX_ORDER, CblasNoTrans, opC, m, n, l, &one, temp->data(), temp->dim(), Cp_start, C._dim, &beta, Dp_start, D._dim);
#endif
}

template<typename ElementType>
void Submatrix<ElementType>::assign(const Submatrix& M, MatrixIndex from)
{
	int row = from.first;
	int col = from.second;
	if (col + M._ncols > this->_ncols || row + M._nrows > this->_nrows)
		throw std::runtime_error("Cannot assign matrix, out of bounds");
	
	//clacpy("A", &m_rows, &m_cols, M._data, &m_dim, this->data(row, col), &this_dim);
	int res = LAPACKE_FUNC(ElementType, lacpy, MATRIX_ORDER, 'A', M.nrows(), M.ncols(), M.data(), M.dim(), this->data(from), dim());
	if (res != 0)
		throw std::runtime_error("Failed LAPACKE_clacpy");
}

template<typename ElementType>
void Submatrix<ElementType>::assign(const ElementType& val_diagonal, const ElementType& val_offdiagonal, MatrixIndex from, int nrows, int ncols)
{
	int row = from.first;
	int col = from.second;
	if (col + ncols > this->_ncols || row + nrows > this->_nrows)
		throw std::runtime_error("Cannot assign matrix, out of bounds");

	LAPACKE_FUNC(ElementType, laset, MATRIX_ORDER, 'A', nrows, ncols, val_offdiagonal, val_diagonal, this->data(from), dim());
}

template<typename ElementType>
void Submatrix<ElementType>::multiply(const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, std::shared_ptr<Submatrix> temp)
{
	std::shared_ptr<Matrix<ElementType>> temp_matrix = nullptr;
	if (temp == nullptr)
	{
		temp_matrix = std::make_shared<Matrix<ElementType>>(n, n);
		temp = std::make_shared<Submatrix>(*temp_matrix);
	}

	temp->assign(0);
	Submatrix<ElementType>::multiply_ex(*this, B, *temp, m, n, k, fromA, fromB, origin, false, false);
	this->assign(*temp, fromA);
}

template<typename ElementType>
void Submatrix<ElementType>::multiply_left(const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, std::shared_ptr<Submatrix> temp)
{
	std::shared_ptr<Matrix<ElementType>> temp_matrix = nullptr;
	if (temp == nullptr)
	{
		temp_matrix = std::make_shared<Matrix<ElementType>>(n, n);
		temp = std::make_shared<Submatrix>(*temp_matrix);
	}

	temp->assign(0);
	Submatrix<ElementType>::multiply_ex(B, *this, *temp, m, n, k, fromB, fromA, origin, false, false);
	this->assign(*temp, fromA);
}

template<typename ElementType>
void Submatrix<ElementType>::multiply(const float scalar, int m, int n, MatrixIndex from)
{
	// TODO: 'G' - full matrix; 'U' - Upper triangular; 'H' - Upper Hessenberg

	ElementType alpha = { scalar };
	if (this->nrows() > from.first + m || this->ncols() > from.second + n)
		throw std::runtime_error("Dimensions mismatch in Submatrix::multiply");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	cilk_for (int row = 0; row < m; row++)
	{
#if USE_CILK
		this->data(from.first + row, from.second)[0:n] *= alpha;
#else
		ElementType* start = (ElementType*)this->data(from.first + row, from.second);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, scal, n, alpha, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, scal, n, &alpha, start, 1);
#endif
#endif
	}

#else
	cilk_for(int col = 0; col < n; col++)
	{
#if USE_CILK
		this->data(from.first, from.second + col)[0:m] *= alpha;
#else
		ElementType* start = (ElementType*)this->data(from.first, from.second + col);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, scal, n, alpha, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, scal, n, &alpha, start, 1);
#endif
#endif
	}
#endif

	//int result = LAPACKE_clascl(MATRIX_ORDER, 'G', 0, 0, 1, scalar, m, n, data(from), _dim);
	//if (result != 0)
	//	throw std::runtime_error("Failed scaling matrix");
}

template<typename ElementType>
void Submatrix<ElementType>::get_scaled_column(const float scalar, int col, int m, int n, MatrixIndex from, std::shared_ptr<Submatrix<ElementType>>& output, MultiplyOption option)
{
	ElementType alpha = { scalar };
	ElementType beta = GetMultiplier(option);
	if (this->nrows() > from.first + m || this->ncols() > from.second + n)
		throw std::runtime_error("Dimensions mismatch in Submatrix::multiply_column");

	if (output->nrows() != m || output->ncols() != 1)
		throw std::runtime_error("Output column has to be m*1 matrix");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	throw std::runtime_error("Not implemented")

#else
#if USE_CILK
	output->data()[0:m] = beta * output->data()[0:m] + this->data(from.first, from.second + col)[0:m] * alpha;
#else
	if (option != None)
		throw new std::runtime_error("It is not yet supported to use option != None without Cilk usage");
	ElementType* start = output->data();
	memcpy(start, this->data(from.first, from.second + col), m * sizeof(ElementType));
#if CALCULATE_REAL
	CBLAS_FUNC(ElementType, scal, n, alpha, start, 1);
#endif
#if CALCULATE_COMPLEX
	CBLAS_FUNC(ElementType, scal, n, &alpha, start, 1);
#endif
#endif
#endif

	//int result = LAPACKE_clascl(MATRIX_ORDER, 'G', 0, 0, 1, scalar, m, n, data(from), _dim);
	//if (result != 0)
	//	throw std::runtime_error("Failed scaling matrix");
}

template<typename ElementType>
void Submatrix<ElementType>::add(const Submatrix& M) { add(M, origin); }

template<typename ElementType>
void Submatrix<ElementType>::add(const Submatrix& M, MatrixIndex from)
{
	ElementType alpha = { 1 };
	if (this->nrows() < from.first + M.nrows() || this->ncols() < from.second + M.ncols())
		throw std::runtime_error("Dimensions mismatch in Submatrix::add");

#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	cilk_for (int row = 0; row < M.nrows(); row++)
	{
#if USE_CILK
		this->data(from.first + row, from.second)[0:M.ncols()] += M.data(row, 0)[0:M.ncols()];
#else
		ElementType* start = (ElementType*)this->data(from.first + row, from.second);
		//cblas_caxpy(M.ncols(), &alpha, M.data(row, 0), 1, start, 1);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.ncols(), alpha, M.data(row, 0), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.ncols(), &alpha, M.data(row, 0), 1, start, 1);
#endif
#endif
	}

#else
	cilk_for(int col = 0; col < M.ncols(); col++)
	{
#if USE_CILK
		this->data(from.first, from.second + col)[0:M.nrows()] += M.data(0, col)[0:M.nrows()];
#else
		ElementType* start = (ElementType*)this->data(from.first, from.second + col);
		//cblas_caxpy(M.ncols(), &alpha, M.data(row, 0), 1, start, 1);
#if CALCULATE_REAL
		CBLAS_FUNC(ElementType, axpy, M.nrows(), alpha, M.data(0, col), 1, start, 1);
#endif
#if CALCULATE_COMPLEX
		CBLAS_FUNC(ElementType, axpy, M.nrows(), &alpha, M.data(0, col), 1, start, 1);
#endif
#endif
	}
#endif
}


template<typename ElementType>
Matrix<ElementType> Matrix<ElementType>::eye(int n)
{
	Matrix I(n, n);
#if USE_CILK
	I.data()[0:n : n + 1] = 1;
#else
	cilk_for (int i = 0; i < n; i++)
	{
		*(I.data(i, i)) = { 1 };
	}
#endif
	return I;
}

//template<typename ElementType>
//Submatrix<ElementType>::operator Matrix<ElementType>() const
//{
//	Matrix<ElementType> I = Matrix<ElementType>::eye(_nrows);
//	Submatrix sI(I);
//
//	return Submatrix::multiply(sI, *this, _nrows, _ncols, _nrows, origin, origin, false, false);
//}

template<typename ElementType>
Matrix<ElementType> Submatrix<ElementType>::polynomial(std::vector<float>& coefficients, MatrixIndex from, int n, PatersonStockmeyerType ps_type)
{
	if (ps_type == VanLoan)
		return polynomial_vanloan(coefficients, from, n);

	int d = coefficients.size(); // Polynomial degree, including I (degree 0)
	size_t p = (size_t)sqrt(d);
	size_t s = (d + 1) / p; // d+1 = ps

	if (coefficients.size() < p*s)
		coefficients.resize(p*s);

	std::vector<ElementType> Aps_buffer((p + 1) * n*n);
	Submatrix<ElementType> I(Aps_buffer.data(), n, n, n);
	I.assign(1, 0); //I.initialize_I();

	//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationInit);
	std::vector<Submatrix<ElementType>> Aps;
	Aps.push_back(I);
	for (size_t i = 1; i < p; i++)
	{
		Submatrix<ElementType> current(Aps_buffer.data() + i * n*n, n, n, n);
		current.assign(0);
		Submatrix::multiply_ex(*this, Aps.back(), current, n, n, n, from, origin, origin, false, false);
		Aps.push_back(current);
	}
	Submatrix<ElementType> Ap(Aps_buffer.data() + p * n*n, n, n, n);
	Submatrix::multiply_ex(*this, Aps.back(), Ap, n, n, n, from, origin, origin, false, false);
	//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationAps);

	Matrix<ElementType> qA(n, n); Submatrix<ElementType> sqA(qA);
	Matrix<ElementType> pPolynomial(n, n); Submatrix<ElementType> sPolynomial(pPolynomial);
	Matrix<ElementType> temp_matrix(n, n); std::shared_ptr<Submatrix<ElementType>> temp = std::make_shared<Submatrix<ElementType>>(temp_matrix);
	for (int si = s - 1; si >= 0; si--) // Reverse order - Horner Rule
	{
		//memset(sPolynomial.data(), 0, n * n * sizeof(ElementType));
		sPolynomial.assign(0);
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationInnerInit);
		if (ps_type == PatersonStockmeyerType::PatersonStockmeyerParallel)
		{
			cilk_for(size_t i = 0; i < Aps.size(); i++)
			{
				int pi = si * p + i;
				sPolynomial.add(Aps[i], coefficients[pi], origin);
				//__sec_reduce_add();
			}
		}
		else if (ps_type == PatersonStockmeyerType::PatersonStockmeyerSequential)
		{
			for (size_t i = 0; i < Aps.size(); i++)
			{
				int pi = si * p + i;
				sPolynomial.add(Aps[i], coefficients[pi], origin);
			}
		}
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationInner);
		// qA = qA * Ap + pPolynomial
		Submatrix<ElementType>::multiply_ex(sqA, Ap, sPolynomial, Ap.nrows(), sqA.ncols(), sqA.nrows(), origin, origin, origin, false, false, Plus, Plus); // pPolynomial += qA*Ap
		sqA.assign(sPolynomial, origin);
		//temp->assign(0);
		//sqA.multiply(Ap, temp);
		//sqA.add(pPolynomial, 1, origin);
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationHorner);
	}
	return qA;
}

template<typename ElementType>
Matrix<ElementType> Submatrix<ElementType>::polynomial_vanloan(std::vector<float>& coefficients, MatrixIndex from, int n)
{
	int d = coefficients.size(); // Polynomial degree, including I (degree 0)
	size_t p = (size_t)sqrt(d);
	size_t s = (d + 1) / p; // d+1 = ps

	if (coefficients.size() < p*s)
		coefficients.resize(p*s);

	std::vector<ElementType> Aps_buffer((p + 1) * n*n);
	Submatrix<ElementType> I(Aps_buffer.data(), n, n, n);
	I.assign(1, 0); //I.initialize_I();

	//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculationInit);
	std::vector<Submatrix<ElementType>> Aps;
	Aps.push_back(I);
	for (size_t i = 1; i < p; i++)
	{
		Submatrix<ElementType> current(Aps_buffer.data() + i * n*n, n, n, n);
		current.assign(0);
		Submatrix::multiply_ex(*this, Aps.back(), current, n, n, n, from, origin, origin, false, false);
		Aps.push_back(current);
	}
	Submatrix<ElementType> Ap(Aps_buffer.data() + p * n*n, n, n, n);
	Submatrix::multiply_ex(*this, Aps.back(), Ap, n, n, n, from, origin, origin, false, false);

	Matrix<ElementType> qA(n, n); Submatrix<ElementType> sqA(qA);
	Matrix<ElementType> pPolynomial(n, n); Submatrix<ElementType> sPolynomial(pPolynomial);
	Matrix<ElementType> Qjs(n, n);
	Matrix<ElementType> temp_matrix(n, n); std::shared_ptr<Submatrix<ElementType>> temp = std::make_shared<Submatrix<ElementType>>(temp_matrix);
	cilk_for(size_t col = 0; col < n; col++)
	{
		Submatrix<ElementType> qAj(qA.data(0, col), n, 1, n);
		Submatrix<ElementType> Apj(Ap.data(0, col), n, 1, n);
		int si = s - 1;
		std::shared_ptr<Submatrix<ElementType>> Qj = std::make_shared<Submatrix<ElementType>>(Qjs.data(1,col), n, 1, n);
		for (; si >= 0; si--) // Reverse order - Horner Rule
		{
			Qj->assign(0);
			cilk_for(size_t i = 0; i < Aps.size(); i++)
			{
				int pi = si * p + i;
				Aps[i].get_scaled_column(coefficients[pi], col, Qj, Plus);
				//__sec_reduce_add();
			}
			// qAj = Ap * qAj + Qj
			Submatrix<ElementType>::multiply_ex(Ap, qAj, *Qj, Ap.nrows(), Qj->ncols(), Qj->nrows(), origin, origin, origin, false, false, Plus, Plus); // Qj += Ap*qAj
			qAj.assign(*Qj, origin);
			//temp->assign(0);
			//qAj.multiply_left(Ap, temp);
			//qAj.add(*Qj, 1, origin);
		}
	}
	return qA;
}

//template<typename ElementType>
//Matrix<ElementType> Matrix<ElementType>::concat_rows(const Matrix<ElementType>& above, const Matrix<ElementType>& below)
//{
//	if (above.ncols() != below.ncols())
//		throw std::runtime_error("Can't concat matrices with incompatible columns");
//
//	Matrix<ElementType> result(above.nrows() + below.nrows(), above.ncols());
//	result.assign(above, origin);
//	result.assign(below, { above.nrows(), 0 });
//	return result;
//}
//
//template<typename ElementType>
//Matrix<ElementType> Matrix<ElementType>::concat_cols(const Matrix<ElementType>& left, const Matrix<ElementType>& right)
//{
//	if (left.nrows() != right.nrows())
//		throw std::runtime_error("Can't concat matrices with incompatible columns");
//
//	Matrix<ElementType> result(left.nrows() , left.ncols() + right.ncols());
//	result.assign(left, origin);
//	result.assign(right, { 0, left.ncols() });
//	return result;
//}


template<typename ElementType>
Matrix<ElementType> Matrix<ElementType>::transpose()
{
	Matrix<ElementType> tr(ncols(), nrows());
	cilk_for(int r = 0; r < nrows(); r++)
	{
		cilk_for(int c = r + 1; c < ncols(); c++)
		{
			*tr.data(c, r) = *data(r, c);
		}
	}
	return tr;
}

extern "C" void RECSY_DTRSYL(char* trana, char* tranb,
	lapack_int* isgn, lapack_int* m, lapack_int* n,
	const double* a, lapack_int* lda, const double* b,
	lapack_int* ldb, double* c, lapack_int* ldc,
	double* scale, int* info);

template<typename ElementType>
void Submatrix<ElementType>::recsy_trsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C)
{
#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	throw std::runtime_error("Should be COL_MAJOR");
#endif
	double scale = 1;
	int info = 0;
	int m = A.nrows();
	int n = B.ncols();
	int lda = A.dim();
	int ldb = B.dim();
	int ldc = C.dim();
	int sign = -1;

	RECSY_DTRSYL("N", "N", &sign, &m, &n, A.data(), &lda, B.data(), &ldb, C.data(), &ldc, &scale, &info);
	if (info < 0)
		throw MatrixException("recsy dtrsyl", info);
}

template<typename ElementType>
void Submatrix<ElementType>::mkl_trsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C)
{
#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	throw std::runtime_error("Should be COL_MAJOR, not supported yet");
#endif
	double scale = 1;
	int info = 0;
	int m = A.nrows();
	int n = B.ncols();
	int lda = A.dim();
	int ldb = B.dim();
	int ldc = C.dim();
	int sign = -1;
	dtrsyl("N", "N", &sign, &m, &n, A.data(), &lda, B.data(), &ldb, C.data(), &ldc, &scale, &info);
	if (info < 0)
		throw MatrixException("mkl dtrsyl", info);
}

template<typename ElementType>
void Submatrix<ElementType>::sylvester(SylvesterFunction sylvester_function, const Submatrix<ElementType>& A, const Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C)
{
	switch (sylvester_function)
	{
	case SylvesterFunction::SylvesterLapack:
		trsylv(A, B, C);
		break;
	case SylvesterFunction::SylvesterRecusrive:
		rtrsylv(A, B, C, 10);
		break;
	case SylvesterFunction::SylvesterRecsy:
		recsy_trsylv(A, B, C);
		break;
	case SylvesterFunction::SylvesterIntel:
		mkl_trsylv(A, B, C);
		break;
	default:
		throw std::runtime_error("Unsupported sylvester solver type");
	}
}

template<typename ElementType>
void Submatrix<ElementType>::trsylv(const Submatrix<ElementType>& A, const Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C)
{
	double scale = 1;
	//auto Ct = ((Matrix<ElementType>)C).transpose();
	int res = LAPACKE_FUNC(ElementType, trsyl, MATRIX_ORDER, 'N', 'N', -1, A.nrows(), B.ncols(), A.data(), A.dim(), B.data(), B.dim(), C.data(), C.dim(), &scale);
	if (res != 0)
		throw MatrixException("ctrsyl", res);
}

template<typename ElementType>
void Submatrix<ElementType>::rtrsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C, int block_size)
{
	int m = A.nrows();
	int n = B.nrows();

	if (A.nrows() != A.ncols())
		throw std::runtime_error("A is not rectangular");
	if (B.nrows() != B.ncols())
		throw std::runtime_error("B is not rectangular");

	if (C.nrows() != m || C.ncols() != n)
		throw std::runtime_error("C has wrong size");

	int half_m = A.get_half_rows();
	int half_n = B.get_half_rows();

	if (n >= 1 && n <= block_size &&
		m >= 1 && m <= block_size)
	{
		double scale = 0;
		int res = LAPACKE_FUNC(ElementType, trsyl, MATRIX_ORDER, 'N', 'N', -1, m, n, A.data(), A.dim(), B.data(), B.dim(), C.data(), C.dim(), &scale);
		//int res = LAPACKE_ctrsyl(MATRIX_ORDER, 'N', 'N', -1, block_m, block_n, Tii.data(), Tii.dim(), Tjj.data(), Tjj.dim(), sF.data(ij), sF.dim(), &scale);
		if (res != 0)
			throw MatrixException("ctrsyl", res);

		return;
	}
	else if (n >= 1 && n <= half_m)
	{
		// Case 1: Split A (by rows and columns), C (by rows only)
		/*
		----| A11 : A12 |	| axa : axb |
		A = |  -  :  -  | = |  -  :  -  |
		----| A21 : A22 |	| bxa : bxb |
		*/
		auto A11 = Submatrix<ElementType>(A.data(0, 0), half_m, half_m, A.dim());
		auto A12 = Submatrix<ElementType>(A.data(0, half_m), half_m, m - half_m, A.dim());
		auto A21 = Submatrix<ElementType>(A.data(half_m, 0), m - half_m, half_m, A.dim());
		auto A22 = Submatrix<ElementType>(A.data(half_m, half_m), m - half_m, m - half_m, A.dim());

		auto C1 = Submatrix<ElementType>(C.data(0, 0), half_m, n, C.dim());
		auto C2 = Submatrix<ElementType>(C.data(half_m, 0), m - half_m, n, C.dim());

		rtrsylv(A22, B, C2, block_size); auto X2 = C2;
		A12.multiply(-1);
		/*C1 = */Submatrix<ElementType>::multiply_ex(A12, X2, C1, A12.nrows(), X2.ncols(), A12.ncols());
		A12.multiply(-1);
		rtrsylv(A11, B, C1, block_size); auto X1 = C1;
		//Submatrix<ElementType>::concat_rows(X1, X2, C);
		return;
	}
	else if (m >= 1 && m <= half_n)
	{
		// Case 2: Split B (by rows and columns), C (by cols only)
		/*
		----| B11 : B12 |	| axa : axb |
		B = |  -  :  -  | = |  -  :  -  |
		----| B21 : B22 |	| bxa : bxb |
		*/
		auto B11 = Submatrix<ElementType>(B.data(0, 0), half_n, half_n, B.dim());
		auto B12 = Submatrix<ElementType>(B.data(0, half_n), half_n, n - half_n, B.dim());
		auto B21 = Submatrix<ElementType>(B.data(half_n, 0), n - half_n, half_n, B.dim());
		auto B22 = Submatrix<ElementType>(B.data(half_n, half_n), n - half_n, n - half_n, B.dim());

		auto C1 = Submatrix<ElementType>(C.data(0, 0), m, half_n, C.dim());
		auto C2 = Submatrix<ElementType>(C.data(0, half_n), m, n - half_n, C.dim());

		rtrsylv(A, B11, C1, block_size); auto X1 = C1;
		/*C2 = */Submatrix<ElementType>::multiply_ex(X1, B12, C2, X1.nrows(), B12.ncols(), X1.ncols());
		rtrsylv(A, B22, C2, block_size); auto X2 = C2;
		//Submatrix<ElementType>::concat_cols(X1, X2, C);
		return;
	}
	else
	{
		// Case 3: Split A, B, and C (all by rows and cols)
		auto A11 = Submatrix<ElementType>(A.data(0, 0), half_m, half_m, A.dim());
		auto A12 = Submatrix<ElementType>(A.data(0, half_m), half_m, m - half_m, A.dim());
		auto A21 = Submatrix<ElementType>(A.data(half_m, 0), m - half_m, half_m, A.dim());
		auto A22 = Submatrix<ElementType>(A.data(half_m, half_m), m - half_m, m - half_m, A.dim());

		auto B11 = Submatrix<ElementType>(B.data(0, 0), half_n, half_n, B.dim());
		auto B12 = Submatrix<ElementType>(B.data(0, half_n), half_n, n - half_n, B.dim());
		auto B21 = Submatrix<ElementType>(B.data(half_n, 0), n - half_n, half_n, B.dim());
		auto B22 = Submatrix<ElementType>(B.data(half_n, half_n), n - half_n, n - half_n, B.dim());

		auto C11 = Submatrix<ElementType>(C.data(0, 0), half_m, half_n, C.dim());
		auto C12 = Submatrix<ElementType>(C.data(0, half_n), half_m, n - half_n, C.dim());
		auto C21 = Submatrix<ElementType>(C.data(half_m, 0), m - half_m, half_n, C.dim());
		auto C22 = Submatrix<ElementType>(C.data(half_m, half_n), m - half_m, n - half_n, C.dim());

		auto C1 = Submatrix<ElementType>(C.data(0, 0), half_m, n, C.dim());
		auto C2 = Submatrix<ElementType>(C.data(half_m, 0), m - half_m, n, C.dim());

		rtrsylv(A22, B11, C21, block_size); auto X21 = C21;
		/*C22 = */cilk_spawn Submatrix<ElementType>::multiply_ex(X21, B12, C22, X21.nrows(), B12.ncols(), X21.ncols());
		A12.multiply(-1);
		/*C11 = */Submatrix<ElementType>::multiply_ex(A12, X21, C11, A12.nrows(), X21.ncols(), A12.ncols());
		A12.multiply(-1);
		cilk_sync;
		cilk_spawn rtrsylv(A22, B22, C22, block_size); auto X22 = C22;
		rtrsylv(A11, B11, C11, block_size); auto X11 = C11;
		cilk_sync;
		/*C12 = */Submatrix<ElementType>::multiply_ex(A12, X22, C12, A12.nrows(), X22.ncols(), A12.ncols());
		/*C12 = */Submatrix<ElementType>::multiply_ex(X11, B12, C12, X11.nrows(), B12.ncols(), X11.ncols());
		rtrsylv(A11, B22, C12, block_size); auto X12 = C12;

		//auto X1 = Matrix<ElementType>::concat_cols(X11, X12);
		//auto X2 = Matrix<ElementType>::concat_cols(X21, X22);
		//return Matrix<ElementType>::concat_rows(X1, X2);
		return;
	}
}

template<typename ElementType>
int Submatrix<ElementType>::get_half_rows() const
{
	return get_half(nrows());
}
template<typename ElementType>
int Submatrix<ElementType>::get_half_cols() const
{
	return get_half(ncols());
}

template<typename ElementType>
int Submatrix<ElementType>::get_half(int n) const
{
	int half = n / 2;
	if (half > 0 && *(data(half, half - 1)) != 0) // check if we are in a middle of 2x2 diagonal block
		half--;
	return half;
}

