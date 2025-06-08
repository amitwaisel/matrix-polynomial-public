#pragma once

#include <vector>
#include <memory>
#include <type_traits>
#include "defs.h"

typedef std::pair<int, int> MatrixIndex;
static const MatrixIndex origin(0, 0);


//#define MKL_FUNC(ElementType, library, func, ...) do {\
//		if (std::is_same<ElementType,MKL_Complex8>::value) library##_c##func(__VA_ARGS__); \
//		if (std::is_same<ElementType,double>::value) library##_d##func(__VA_ARGS__); \
//		if (std::is_same<ElementType,float>::value) library##_s##func(__VA_ARGS__); \
//	} while(false)
#if CALCULATE_REAL
#define NATIVE_FUNC(func, ...) d##func(__VA_ARGS__);
#define MKL_FUNC(library, func, ...) library##_d##func(__VA_ARGS__);
#define VECTOR_FUNC(func, ...) vd##func(__VA_ARGS__)
#endif
#if CALCULATE_COMPLEX
#define NATIVE_FUNC(func, ...) c##func(__VA_ARGS__);
#define MKL_FUNC(library, func, ...) library##_c##func(__VA_ARGS__);
#define VECTOR_FUNC(func, ...) vc##func(__VA_ARGS__)
#endif

#define CBLAS_FUNC(ElementType, func, ...) MKL_FUNC(cblas, func, __VA_ARGS__)
//
//#if USE_LAPACKE && MATRIX_ORDER == LAPACK_COL_MAJOR
//#define ARG(x) x
#define LAPACKE_FUNC(ElementType, func, ...) MKL_FUNC(LAPACKE, func, __VA_ARGS__)
//#else
//#define ARG(x) &x
//#define LAPACKE_FUNC(ElementType, func, ...) do {\
//	int result = 0;\
//	std::vector<ElementType> work(n);\
//	NATIVE_FUNC(func, __VA_ARGS__, &work, &info); \
//} while(false)
//#endif

#define CONVERT_REAL(ElementType, output, data, count) for(unsigned int i = 0; i < count; i++) output[i] = data[i].real

class MatrixException : public std::runtime_error
{
public:
	MatrixException(char* func_name, int lapack_result) : std::runtime_error(func_name), _failed_argument(lapack_result * -1), _function(func_name) {}
	virtual ~MatrixException() throw() = default;


private:
	int _failed_argument;
	std::string _function;
};

enum SylvesterFunction
{
	SylvesterLapack = 0,
	SylvesterRecusrive,
	SylvesterIntel,
	SylvesterRecsy,
};

enum PatersonStockmeyerType
{
	PatersonStockmeyerSequential = 0,
	PatersonStockmeyerParallel,
	VanLoan,
};

template<typename ElementType>
class Matrix
{
public:
	explicit Matrix(int nrows, int ncols);
	Matrix(const ElementType* data, int nrows, int ncols);
	//Matrix(const MKL_Complex8* data, int nrows, int ncols, int dim);
	//Matrix(const Matrix& other);
	Matrix(const Matrix& other);
	Matrix& operator=(const Matrix& other) = default;
	Matrix(Matrix&& other);
	Matrix& operator=(Matrix&& other);
	~Matrix();

	const ElementType* data() const { return _data.data(); }
	ElementType* data() { return _data.data(); }
#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	const ElementType* data(int row, int col) const { return _data.data() + row * _ncols + col; }
	const ElementType* data(const MatrixIndex& i) const { return _data.data() + i.first * _ncols + i.second; }
	ElementType* data(int row, int col) { return _data.data() + row * _ncols + col; }
	ElementType* data(const MatrixIndex& i) { return _data.data() + i.first * _ncols + i.second; }
	int dim() const { return ncols(); }
#else
	const ElementType* data(int row, int col) const { return _data.data() + row + col * _nrows; }
	const ElementType* data(const MatrixIndex& i) const { return _data.data() + i.first + i.second * _nrows; }
	ElementType* data(int row, int col) { return _data.data() + row + col * _nrows; }
	ElementType* data(const MatrixIndex& i) { return _data.data() + i.first + i.second * _nrows; }
	int dim() const { return nrows(); }
#endif
	int nrows() const { return _nrows; }
	int ncols() const { return _ncols; }
	size_t size() const { return _data.size() * sizeof(ElementType); }

	void reset();
	void reset(const Matrix& other);

	// Add whole matrix M to location 'from' in 'this'
	void add(const Matrix& M, float scalar, MatrixIndex from);
	void add(const Matrix& M) { add(M, 1, origin); }
	void assign(const Matrix& M, MatrixIndex from);
	Matrix<ElementType> transpose();

	void multiply(const Matrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB);
	void multiply(const Matrix& B);
	static Matrix multiply(const Matrix& A, const Matrix& B);
	static Matrix multiply(const Matrix& A, const Matrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB);

	static Matrix eye(int n);

	//static Matrix<ElementType> concat_rows(const Matrix<ElementType>& above, const Matrix<ElementType>& below);
	//static Matrix<ElementType> concat_cols(const Matrix<ElementType>& left, const Matrix<ElementType>& right);

	void polynomial(std::vector<float>& coefficients, MatrixIndex from, int n);

private:
	std::vector<ElementType> _data;
	int _nrows, _ncols;
	MatrixIndex _from;
	MatrixIndex _to;

};

template<typename ElementType>
class Submatrix
{
public:
	enum MultiplyOption
	{
		None = 0,
		Plus,
		Minus
	};
public:
	Submatrix(Matrix<ElementType>& m) : _data(m.data()), _nrows(m.nrows()), _ncols(m.ncols()), _dim(m.dim()) {}
	Submatrix(const Matrix<ElementType>& m) : Submatrix(const_cast<Matrix<ElementType>&>(m)) {}
	Submatrix(ElementType* data, int nrows, int ncols, int dim) : _data(data), _nrows(nrows), _ncols(ncols), _dim(dim) {}
	//explicit operator Matrix<ElementType>() const;
#if MATRIX_ORDER == LAPACK_ROW_MAJOR
	ElementType* data(MatrixIndex& from) { return &_data[from.first * _dim + from.second]; }
	ElementType* data(int row, int col) { return &_data[row * _dim + col]; }
	const ElementType* data(int row, int col) const { return &_data[row * _dim + col]; }
#else
	ElementType* data(MatrixIndex& from) { return &_data[from.first + from.second * _dim]; }
	ElementType* data(int row, int col) { return &_data[row + col * _dim]; }
	const ElementType* data(int row, int col) const { return &_data[row + col * _dim]; }
#endif
	ElementType* data() { return data(0, 0); }
	const ElementType* data() const { return data(0, 0); }

	void initialize_I();

	// [in] A:m*k, [in] B:k*n, [ret] C:m*n [C=A*B]
	//static Matrix<ElementType> multiply(const Submatrix& A, const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, bool transA, bool transB);
	// A:m*n, B:n*n, out:m*m [A*B*A^T]
	// [in] A:m*k, [in] B:k*n, [out] C:m*n [C = (0+-)A*B (0+-)C]
	static void multiply_ex(const Submatrix& A, const Submatrix& B, Submatrix& C,
		int m, int n, int k,
		MatrixIndex fromA = origin, MatrixIndex fromB = origin, MatrixIndex fromC = origin,
		bool transA = false, bool transB = false,
		MultiplyOption option_AB = Plus, MultiplyOption option_C = Plus);
	// [in] A:m*k, [in] B:k*l, [in] C:l*n, [out] D:m*n [D = (0+-)A*B*C (0+-)D]
	static void multiply3_ex(
		const Submatrix& A, const Submatrix& B,	const Submatrix& C,	Submatrix& D,
		int m, int n, int k, int l,
		MatrixIndex fromA, MatrixIndex fromB, MatrixIndex fromC, MatrixIndex fromD,
		bool transA = false, bool transB = false, bool transC = false,
		MultiplyOption option_ABC = Plus, MultiplyOption option_D = Plus,
		std::shared_ptr<Submatrix<ElementType>> temp = nullptr);
	void multiply(const Submatrix& B, std::shared_ptr<Submatrix> temp = nullptr) { multiply(B, nrows(), B.ncols(), ncols(), origin, origin, temp); }
	void multiply_left(const Submatrix& B, std::shared_ptr<Submatrix> temp = nullptr) { multiply_left(B, B.nrows(), ncols(), B.ncols(), origin, origin, temp); }
	void multiply(const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, std::shared_ptr<Submatrix> temp = nullptr);
	void multiply_left(const Submatrix& B, int m, int n, int k, MatrixIndex fromA, MatrixIndex fromB, std::shared_ptr<Submatrix> temp = nullptr);
	void multiply(const float scalar, int m, int n, MatrixIndex from);
	void multiply(const float scalar) { multiply(scalar, nrows(), ncols(), origin); }
	void get_scaled_column(const float scalar, int col, int m, int n, MatrixIndex from, std::shared_ptr<Submatrix<ElementType>>& output, MultiplyOption option);
	void get_scaled_column(const float scalar, int col, std::shared_ptr<Submatrix<ElementType>>& output, MultiplyOption option) { return get_scaled_column(scalar, col, nrows(), ncols(), origin, output, option); }
	void assign(const Submatrix& M, MatrixIndex from);
	void assign(const ElementType& val) { assign(val, val); }
	void assign(const ElementType& val_diagonal, const ElementType& val_offdiagonal) { assign(val_diagonal, val_offdiagonal, origin, nrows(), ncols()); }
	void assign(const ElementType& val_diagonal, const ElementType& val_offdiagonal, MatrixIndex from, int nrows, int ncols);
	void add(const Submatrix& M);
	void add(const Submatrix& M, MatrixIndex from);
	void add(const Submatrix& M, float scalar, MatrixIndex from);

	int get_half_rows() const;
	int get_half_cols() const;

	// Coefficients in ascending order - c_0, c_1, ..., c_(d-1), c_d
	Matrix<ElementType> polynomial(std::vector<float>& coefficients, MatrixIndex from, int n, PatersonStockmeyerType ps_type = PatersonStockmeyerType::PatersonStockmeyerParallel);
	Matrix<ElementType> polynomial_vanloan(std::vector<float>& coefficients, MatrixIndex from, int n);

	static void trsylv(const Submatrix<ElementType>& A, const Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C);
	static void rtrsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C, int block_size = 10);
	static void recsy_trsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C);
	static void mkl_trsylv(Submatrix<ElementType>& A, Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C);
	static void sylvester(SylvesterFunction sylvester_function, const Submatrix<ElementType>& A, const Submatrix<ElementType>& B, /*in out*/ Submatrix<ElementType>& C);

	int nrows() const { return _nrows; }
	int ncols() const { return _ncols; }
	int dim() const { return _dim; }

private:
	int get_half(int n) const;
	static ElementType GetMultiplier(MultiplyOption option);

private:
	ElementType* _data;
	int _nrows, _ncols, _dim;
};

template<typename ElementType>
class MatrixBlock
{
public:
	MatrixBlock(const Submatrix<ElementType>& block, MatrixIndex from, int rows, int cols) : _block(block), _from(from), _rows(rows), _cols(cols) {}
	MatrixBlock(const Submatrix<ElementType>& block, MatrixIndex from, int n) : MatrixBlock(block, from, n, n) {}
	Submatrix<ElementType>* matrix() { return &this->_block; };
	MatrixIndex from() const { return _from; }
	int dim() const { return _block.dim(); }
	int nrows() const { return _rows; }
	int ncols() const { return _cols; }
	const ElementType* data() const { return matrix().data(from()); }


private:
	Submatrix<ElementType> _block;
	MatrixIndex _from;
	int _rows, _cols;
};

#include "Matrix.cpp"
