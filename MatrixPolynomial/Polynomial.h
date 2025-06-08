#pragma once
#include <mkl.h>
#include "Matrix.h"
#include "SmartSwap.h"
#include <functional>
#include <memory>

enum PolynomialBlockAlgorithm
{
	ParlettRecurrence = 0,
	ParallelParlettRecurrence,
	Sequential,
};

enum SchurDecompositionAlgorithm
{
	Default = 0,
	Hassenberg,
};

enum EigenvaluesClusterSortAlgorithm
{
	BubbleSort = 0,
	SmartSwapSort,
};


template<typename ElementType>
class Polynomial
{
	typedef std::function<void(Submatrix<ElementType>& A, int n, MatrixIndex from, Submatrix<ElementType>& output)> matrix_function;
public:
	Polynomial(const Matrix<ElementType>& A);
	~Polynomial();

	Matrix<ElementType> calculatePolynomial(const std::vector<float>& coefficients,
		double eigenvalues_cluster_delta = 0.1,
		SchurDecompositionAlgorithm schur_algorithm = SchurDecompositionAlgorithm::Default,
		EigenvaluesClusterSortAlgorithm eigenvalues_sort_algorithm = EigenvaluesClusterSortAlgorithm::BubbleSort,
		PolynomialBlockAlgorithm block_algorithm = PolynomialBlockAlgorithm::ParlettRecurrence,
		SylvesterFunction sylvester_function = SylvesterFunction::SylvesterLapack
	);

private:
	typedef std::vector<ElementType> Eigenvalues;

private:
	static void getSchurDecomposition(const Matrix<ElementType>& A, Matrix<ElementType>& schur, Matrix<ElementType>& schur_co, Eigenvalues& eigenvalues);
	static void getSchurDecompositionHassenberg(const Matrix<ElementType>& A, Matrix<ElementType>& schur, Matrix<ElementType>& schur_co, Eigenvalues& eigenvalues);
	
	static std::vector<int> clusterValues(const ElementType* values, int nvalues, double delta);
	static std::vector<int> chooseClusterPermutation(const std::vector<int>& clusters);
	static std::vector<int> applyPermutation(const std::vector<int>& clusters, const std::vector<int>& permutation);
	static void buildClusterBlocks(Submatrix<ElementType> schur, const std::vector<int>& clusters, const std::vector<int>& permutation, SmartSwap::ClustersWithBlocks& cluster_blocks, SmartSwap::ClustersWithBlocks & ordered_blocks);
	static std::vector<int> bubbleSortEigenvalues(Submatrix<ElementType>& schur, Submatrix<ElementType>& schur_co, int n, std::vector<int>& clusters, const std::vector<int>& permutation);
	static std::vector<int> clusterSchurEigenvalues(Submatrix<ElementType>& schur, Submatrix<ElementType>& schur_co, Eigenvalues& eigenvalues, double tolerance, EigenvaluesClusterSortAlgorithm eigenvalues_sort_algorithm);

	static Matrix<ElementType> calculateDiagonalBlocks(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f);
	static std::vector<MatrixBlock<ElementType>> calculateDiagonalBlocks(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, Submatrix<ElementType>& F);
	
	static void calculateBlockInRecurrence(Submatrix<ElementType>& T, int block_row, int block_col, std::vector<MatrixBlock<ElementType>>& blocks, SylvesterFunction sylvester_function, Submatrix<ElementType>& F);
	static std::shared_ptr<Matrix<ElementType>> calculateParallelBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function);
	
	static std::shared_ptr<Matrix<ElementType>> calculateSuperdiagonalBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function);
	static std::shared_ptr<Matrix<ElementType>> calculateBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function);
	static std::shared_ptr<Matrix<ElementType>> calculateBlockedFunction(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function);

private:
	Matrix<ElementType> raw_input;
	Submatrix<ElementType> input;
	Matrix<ElementType> _schur;
	Matrix<ElementType> _schur_co;
};

#include "Polynomial.cpp"

