#include "Polynomial.h"
#include "Benchmark.h"
#include "SmartSwap.h"
#include <map>
#include <algorithm>
#include <functional>

#define UNREFERENCED_PARAMETER(P) (P)

template<typename ElementType>
Polynomial<ElementType>::Polynomial(const Matrix<ElementType>& A) : raw_input(A), input(A), _schur(A.nrows(), A.ncols()), _schur_co(A.nrows(), A.ncols())
{
	if (A.nrows() != A.ncols())
		throw std::runtime_error("Input matrix for polynomial calculations have to be rectangular");
}


template<typename ElementType>
Polynomial<ElementType>::~Polynomial()
{
}

template<typename ElementType>
Matrix<ElementType> Polynomial<ElementType>::calculatePolynomial(const std::vector<float>& coefficients, double eigenvalues_cluster_delta,
	SchurDecompositionAlgorithm schur_algorithm, EigenvaluesClusterSortAlgorithm eigenvalues_sort_algorithm, PolynomialBlockAlgorithm block_algorithm, SylvesterFunction sylvester_function)
{
	Submatrix<ElementType> sub_schur = (Submatrix<ElementType>)_schur;
	Submatrix<ElementType> sub_schur_co = (Submatrix<ElementType>)_schur_co;
	Eigenvalues eigenvalues;
	Benchmark::Instance().Start();

	Benchmark::Instance().TraceCurrentStep(Benchmark::CheckpointType::SchurDecomposition);
	if (schur_algorithm == SchurDecompositionAlgorithm::Default)
		getSchurDecomposition(raw_input, _schur, _schur_co, eigenvalues);
	else if (schur_algorithm == SchurDecompositionAlgorithm::Hassenberg)
		getSchurDecompositionHassenberg(raw_input, _schur, _schur_co, eigenvalues);
	else
		throw std::runtime_error("Invalid schur algorithm");
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::SchurDecomposition);

	std::vector<int> clusters = clusterSchurEigenvalues(sub_schur, sub_schur_co, eigenvalues, eigenvalues_cluster_delta, eigenvalues_sort_algorithm);

	std::vector<float> co(coefficients);
	matrix_function f = [&co](Submatrix<ElementType>& A, int n, MatrixIndex from, Submatrix<ElementType>& output)
	{
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculation,
		//	[&A, &co, &from, &n, &output]() {
			Matrix<ElementType> qA = A.polynomial(co, from, n);
			output.assign(qA, from);
		//});
	};
	std::shared_ptr<Matrix<ElementType>> Ft = nullptr;
	Benchmark::Instance().TraceCurrentStep(Benchmark::CheckpointType::TotalParlettRecurrence);
	if (block_algorithm == ParlettRecurrence)
		Ft = calculateBlockParlettRecurrence(sub_schur, clusters, f, sylvester_function);
	else if (block_algorithm == ParallelParlettRecurrence)
		Ft = calculateSuperdiagonalBlockParlettRecurrence(sub_schur, clusters, f, sylvester_function);
	else if (block_algorithm == Sequential)
		Ft = calculateBlockedFunction(sub_schur, clusters, f, sylvester_function);

	// Z*Ft*Z'
	Benchmark::Instance().TraceCurrentStep(Benchmark::CheckpointType::FinalMultiplication);
	int n = input.nrows();
	Matrix<ElementType> F(n, n);
	Submatrix<ElementType> sub_F(F);
	Submatrix<ElementType>::multiply3_ex(sub_schur_co, *Ft, sub_schur_co, sub_F, n, n, n, n, origin, origin, origin, origin, false, false, true, Submatrix<ElementType>::MultiplyOption::Plus);
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::FinalMultiplication);
	return F;
}

template<typename ElementType>
void Polynomial<ElementType>::getSchurDecomposition(const Matrix<ElementType>& A, Matrix<ElementType>& schur, Matrix<ElementType>& schur_co, Eigenvalues& eigenvalues)
{
	// Other things to check:
	// * ?gees
	// * Order schur - http://www.netlib.org/lapack/lawnspdf/lawn171.pdf, https://nickhigham.files.wordpress.com/2017/03/kjel_cse17.pdf, http://www8.cs.umu.se/research/uminf/reports/2017/011/part1.pdf

	int n = A.dim();
	int leading_dimension = n > 1 ? n : 1;
	int result = 0;

	UNREFERENCED_PARAMETER(eigenvalues);

	int sdim = 0;
	eigenvalues.resize(n);
	schur.reset(A);
	// Ignore the imaginary part of the eigenvalues for now

	// Real flavor
#if CALCULATE_REAL
	Eigenvalues imaginary_eigenvalues(n);
#if USE_LAPACKE
	result = LAPACKE_FUNC(ElementType, gees, MATRIX_ORDER, 'V', 'N', nullptr, n, schur.data(), leading_dimension, &sdim, eigenvalues.data(), imaginary_eigenvalues.data(), schur_co.data(), leading_dimension);
#else
	int lwork = -1;
	ElementType work_size = 0;
	dgees("V", "N", nullptr, &n, schur.data(), &leading_dimension, &sdim, eigenvalues.data(), imaginary_eigenvalues.data(), schur_co.data(), &leading_dimension, &work_size, &lwork, nullptr, &result);
	lwork = (int)work_size;
	std::vector<ElementType> work(lwork);
	dgees("V", "N", nullptr, &n, schur.data(), &leading_dimension, &sdim, eigenvalues.data(), imaginary_eigenvalues.data(), schur_co.data(), &leading_dimension, work.data(), &lwork, nullptr, &result);
#endif
#endif

#if CALCULATE_COMPLEX
#if USE_LAPACKE
	result = LAPACKE_FUNC(ElementType, gees, MATRIX_ORDER, 'V', 'N', nullptr, n, schur.data(), leading_dimension, &sdim, eigenvalues.data(), schur_co.data(), leading_dimension);
#else
	int lwork = -1;
	ElementType work_size = 0;
	cgees("V", "N", nullptr, &n, schur.data(), &leading_dimension, &sdim, eigenvalues.data(), schur_co.data(), &leading_dimension, &work_size, &lwork, rwork.data(), nullptr, &result);
	lwork = (int)work_size;
	std::vector<ElementType> work(lwork);
	std::vector<float> rwork(n);
	cgees("V", "N", nullptr, &n, schur.data(), &leading_dimension, &sdim, eigenvalues.data(), schur_co.data(), &leading_dimension, work.data(), &lwork, rwork.data(), nullptr, &result);
#endif
#endif

	if (result != 0)
		throw MatrixException("cgees", result);
	
	//// Sanity test
	//Matrix<ElementType> FF(n, n);
	//Submatrix<ElementType> sub_FF(FF);
	//Submatrix<ElementType>::multiply3_ex(schur_co, schur, schur_co, sub_FF, n, n, n, n, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, false, false, true, true);
}
// @param T		IN		input matrix T
// @param n		IN		The order of the matrix T
// @param Z		OUT		unitary matrix of Schur vectors Z*schur*Z^H
// @param schur	OUT		SCHUR form Z*schur*Z^H
template<typename ElementType>
void Polynomial<ElementType>::getSchurDecompositionHassenberg(const Matrix<ElementType>& A, Matrix<ElementType>& schur, Matrix<ElementType>& schur_co, Eigenvalues& eigenvalues)
{
	// http://www.netlib.org/lapack/lug/node50.html
	/*
	1.	A general matrix A is reduced to upper Hessenberg form H which is zero below the first subdiagonal.
	The reduction may be written A=QHQT with Q orthogonal if A is real, or A=QHQH with Q unitary if A is complex.
	The reduction is performed by subroutine xGEHRD, which represents Q in a factored form, as described in section 5.4.
	The routine xORGHR (or in the complex case xUNGHR) is provided to form Q explicitly.
	The routine xORMHR (or in the complex case xUNMHR) is provided to multiply another matrix by Q without forming Q explicitly.
	2.	The upper Hessenberg matrix H is reduced to Schur form T, giving the Schur factorization H=STST (for H real) or H=STSH (for H complex).
	The matrix S (the Schur vectors of H) may optionally be computed as well.
	Alternatively S may be postmultiplied into the matrix Q determined in stage 1, to give the matrix Z = Q S, the Schur vectors of A.
	The eigenvalues are obtained from the diagonal of T.
	All this is done by subroutine xHSEQR.
	*/
	int result = 0;

	// ?gebal	Balances a general matrix to improve the accuracy of computed eigenvalues and eigenvectors.
	// ?gehrd	Reduces a general matrix to upper Hessenberg form.


	// TODO: LDA? ROW/COL order [ the distance in elements between two neighbor elements in a line of minor dimension ]
	int ilo = 0, ihi = 0;
	int n = A.dim();
	std::vector<double> scale(n);
	int leading_dimension = n > 1 ? n : 1;

	Matrix<ElementType> H(A);

#define DONT_PERMUTE_DONT_SCALE	'N' // None
#define PERMUTE_DONT_SCALE		'P' // Permute
#define DONT_PERMUTE_SCALE		'S' // Scale
#define PERMUTE_SCALE			'B' // Both
	result = LAPACKE_FUNC(ElementType, gebal, MATRIX_ORDER, PERMUTE_SCALE, n, H.data(), leading_dimension, &ilo, &ihi, scale.data());
	if (result != 0)
		throw MatrixException("cgebal", result);

	// The routine does not form the matrix Q explicitly. Instead, Q is represented as a product of elementary reflectors.
	std::vector<ElementType> Q(n);
	result = LAPACKE_FUNC(ElementType, gehrd, MATRIX_ORDER, n, ilo, ihi, H.data(), leading_dimension, Q.data());
	if (result != 0)
		throw MatrixException("cgehrd", result);

	schur.reset(H);

#define EIGENVALUES_ONLY		'E'
#define EIGENVALUES_AND_SCHUR	'S'
#define NO_SCHUR_VECTORS			'N'
#define SCHUR_VECTORS_HESSENBERG	'I'
#define SCHUR_VECTORS_ORIGINAL		'V'
	// Calculate SCHUR factorization, with Schur vectors for Hessenberg form. Those vectors will be multiplied to Q to form Z [out].
	eigenvalues.resize(n);
	std::vector<ElementType> temp_Z(leading_dimension*n);
#if CALCULATE_REAL
	Eigenvalues imaginary_eigenvalues(n);
	result = LAPACKE_FUNC(ElementType, hseqr, MATRIX_ORDER, EIGENVALUES_AND_SCHUR, SCHUR_VECTORS_HESSENBERG, n, ilo, ihi, schur.data(), leading_dimension, eigenvalues.data(), imaginary_eigenvalues.data(), schur_co.data(), leading_dimension);
#else
	result = LAPACKE_FUNC(ElementType, hseqr, MATRIX_ORDER, EIGENVALUES_AND_SCHUR, SCHUR_VECTORS_HESSENBERG, n, ilo, ihi, schur.data(), leading_dimension, eigenvalues.data(), schur_co.data(), leading_dimension);
#endif
	if (result != 0)
		throw MatrixException("chseqr", result);

	// Calculate the real Z=Q*temp_Z
#define LEFT_SIDE	'L'
#define RIGHT_SIDE	'R'
#define NORMAL		'N'
#define TRANSPOSE	'T'
#if CALCULATE_REAL
	result = LAPACKE_FUNC(ElementType, ormhr, MATRIX_ORDER, LEFT_SIDE, NORMAL, n, n, ilo, ihi, H.data(), leading_dimension, Q.data(), schur_co.data(), leading_dimension);
#else
	result = LAPACKE_FUNC(ElementType, unmhr, MATRIX_ORDER, LEFT_SIDE, NORMAL, n, n, ilo, ihi, H.data(), leading_dimension, Q.data(), schur_co.data(), leading_dimension);
#endif
	if (result != 0)
		throw MatrixException("cunmhr or dormhr", result);
}

// Implemented according to Higham/Davis algorithm 9.5 (page 227)
template<typename ElementType>
std::vector<int> Polynomial<ElementType>::clusterValues(const ElementType* values, int nvalues, double delta)
{
	int n = nvalues;
	int p = 1;
	std::vector<int> clusters(n);

	for (int i = 0; i < n; i++)
	{
		ElementType lambda_i = values[i];

		if (clusters[i] == 0 || clusters[i] >= p) // q == p
		{
			clusters[i] = p;
			p++;
		}

		int qi = clusters[i];
		for (int j = i + 1; j < n; j++)
		{
			ElementType lambda_j = values[j];
			if (clusters[j] != clusters[i]) // if lambda_j \notin S_qi
			{
				double distance = 0;
				//ElementType diff = { lambda_j.real - lambda_i.real, lambda_j.imag - lambda_i.imag }; // lambda_j - lambda_i
				//vcAbs(1, &diff, &distance); // |lambda_j - lambda_i|
				ElementType diff = lambda_j - lambda_i;
				VECTOR_FUNC(Abs, 1, &diff, &distance); // |lambda_j - lambda_i|
				if (distance <= delta)
				{
					qi = clusters[i];
					if (clusters[j] == 0 || clusters[j] >= p)
					{
						clusters[j] = qi;
					}
					else
					{ // TODO: Test this
						int qj = clusters[j];

						int min_qiqj = (qi < qj) ? qi : qj;
						int max_qiqj = (qi > qj) ? qi : qj;

						cilk_for (int m = 0; m < n; m++)
						{
							if (clusters[m] == max_qiqj)
								clusters[m] = min_qiqj; // Move the elements of Smax(qi,qj ) to Smin(qi,qj )
							if (clusters[m] > max_qiqj)
								clusters[m]--; // Reduce by 1 the indices of sets Sq for q > max(qi, qj ).
						}
						p--;
					}
				}
			}
		}
	}
	return clusters;
}

template<typename ElementType>
std::vector<int> Polynomial<ElementType>::chooseClusterPermutation(const std::vector<int>& clusters)
{
	int n = clusters.size();
	std::vector<int> permutation(n);
	std::map<int, std::pair<int, int>> clusters_index;
	std::map<int, float> cluster_average;
	for (int i = 0; i < n; i++)
	{
		int cluster = clusters[i];
		if (clusters_index.find(cluster) == clusters_index.end())
			clusters_index[cluster] = std::pair<int, float>(0, 0);
		auto index = &clusters_index[cluster];
		index->first++;
		index->second += (i + 1);

		cluster_average[cluster] = index->second / (float)index->first;
	}
	int nclusters = clusters_index.size();
	auto ordered_clusters = std::vector<int>(nclusters);
	cilk_for (int i = 0; i < nclusters; i++)
		ordered_clusters[i] = i + 1;
	//ordered_clusters.data()[0:nclusters] = 0;
	std::sort(ordered_clusters.begin(), ordered_clusters.end(), [&cluster_average](const int& clus1, const int& clus2) {
		return cluster_average[clus1] < cluster_average[clus2];
	});
	return ordered_clusters;
}

template<typename ElementType>
std::vector<int> Polynomial<ElementType>::applyPermutation(const std::vector<int>& clusters, const std::vector<int>& permutation)
{
	int n = clusters.size();
	std::map<int, int> clusters_counter;
	for (int i = 0; i < n; i++)
	{
		int cluster = clusters[i];
		if (clusters_counter.find(cluster) == clusters_counter.end())
			clusters_counter[cluster] = 0;
		clusters_counter[cluster]++;
	}
	std::vector<int> result;
	for (int cluster : permutation)
	{
		for (int i = 0; i < clusters_counter[cluster]; i++)
			result.push_back(cluster);
	}
	return result;
}

template<typename ElementType>
void Polynomial<ElementType>::buildClusterBlocks(Submatrix<ElementType> schur, const std::vector<int>& clusters, const std::vector<int>& permutation, SmartSwap::ClustersWithBlocks& cluster_blocks, SmartSwap::ClustersWithBlocks & ordered_blocks)
{
	std::vector<int> permutation_indexes(permutation.size());
	for (unsigned int i = 0; i < permutation.size(); i++)
		permutation_indexes[permutation[i] - 1] = i;


	std::map<int, int> blocks_count;
	std::map<int, int> single_count;

	int total_single_count = 0;
	int total_block_count = 0;

	for (int i = 0; i < (int)clusters.size(); i++)
	{
		if (blocks_count.count(clusters[i]) == 0)
			blocks_count[clusters[i]] = 0;
		if (single_count.count(clusters[i]) == 0)
			single_count[clusters[i]] = 0;

		cluster_blocks.clusters.push_back(clusters[i]);
		if (i < schur.nrows() - 1 && *(schur.data(i + 1, i)) != 0) // 2*2 block
		{
			cluster_blocks.blocks.push_back(2);
			blocks_count[clusters[i]]++;
			total_block_count++;
			i++;
		}
		else
		{
			cluster_blocks.blocks.push_back(1);
			single_count[clusters[i]]++;
			total_single_count++;
		}
	}

	Benchmark::Instance().AddProperty(Benchmark::PropertyType::SingleBlocksCount, total_single_count);
	Benchmark::Instance().AddProperty(Benchmark::PropertyType::DoubleBlocksCount, total_block_count);

	ordered_blocks.clusters = SmartSwap::Clusters(cluster_blocks.clusters);
	ordered_blocks.blocks = SmartSwap::Clusters(cluster_blocks.blocks);
	
	// Bubble sort
	int temp_n = cluster_blocks.clusters.size();
	while (temp_n > 0)
	{
		int new_n = 0;
		for (int i = 1; i < temp_n; i++)
		{
			int cluster_i_1 = ordered_blocks.clusters[i - 1];
			int cluster_i = ordered_blocks.clusters[i];
			int block_i_1 = ordered_blocks.blocks[i - 1];
			int block_i = ordered_blocks.blocks[i];

			if (permutation_indexes[cluster_i_1 - 1] > permutation_indexes[cluster_i - 1]) // We have to rearrange by clusters
			{
				ordered_blocks.clusters[i - 1] = cluster_i;
				ordered_blocks.clusters[i] = cluster_i_1;
				ordered_blocks.blocks[i - 1] = block_i;
				ordered_blocks.blocks[i] = block_i_1;

				new_n = i;
			}
		}
		temp_n = new_n;
	}
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesSort);
	Benchmark::Instance().TraceCurrentStep(Benchmark::CheckpointType::EigenvaluesSort);
}

template<typename ElementType>
std::vector<int> Polynomial<ElementType>::bubbleSortEigenvalues(Submatrix<ElementType>& schur, Submatrix<ElementType>& schur_co, int n, std::vector<int>& clusters, const std::vector<int>& permutation)
{
	std::vector<int> permutation_indexes(permutation.size());
	for (unsigned int i = 0; i < permutation.size(); i++)
		permutation_indexes[permutation[i] - 1] = i;

	int leading_dimension = n > 1 ? n : 1;

	// Bubble sort
	int temp_n = n;
	int swaps_count = 0, total_swaps = 0;
	while (temp_n > 0)
	{
		int new_n = 0;
		for (int i = 1; i < temp_n; i++)
		{
			int cluster_i_1 = clusters[i - 1];
			int cluster_i = clusters[i];

			if (permutation_indexes[cluster_i_1 - 1] > permutation_indexes[cluster_i - 1]) // We have to rearrange by clusters
			{
				// Replace in T

				int ifst = i - 1 + 1;
				int ilst = i + 1;
				int orig_ifst = ifst, orig_ilst = ilst;
				int ifst_cluster = clusters[ifst - 1];
				int ilst_cluster = clusters[ilst - 1];
				//dtrsen
				//
				// Complex flavor
				//int result = LAPACKE_FUNC(ElementType, trexc, MATRIX_ORDER, 'V', n, schur.data(), leading_dimension, schur_co.data(), leading_dimension, ifst, ilst);
				// Real flavor
				Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesSort);

#if USE_LAPACKE
				int result = LAPACKE_FUNC(ElementType, trexc, MATRIX_ORDER, 'V', n, schur.data(), leading_dimension, schur_co.data(), leading_dimension, &ifst, &ilst);
#else
				int result = 0;
				std::vector<ElementType> work(n);
				dtrexc("V", &n, schur.data(), &leading_dimension, schur_co.data(), &leading_dimension, &ifst, &ilst, work.data(), &result);
#endif
				Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::dtrexc);
				
				swaps_count++;
				total_swaps += ilst - ifst;

				if (result != 0)
					throw MatrixException("ctrexc", result);
				// Replace clusters registration

				// a,b,b -> b,b,a => 1<->2 -> 1<->3
				// a,a,b -> b,a,a => 2<->3 -> 1<->2
				// a,a,b,b -> b,b,a,a => 2<->3 -> 1<->3
				// If ifst--, first block was 2*2
				// If ilst--, second block was 1*1 (if not changed or ++, second block was 2*2)
				clusters[ifst - 1] = ilst_cluster;
				clusters[ilst - 1] = ifst_cluster;
				if (ifst == orig_ifst - 1) // first block (now in ilst) was 2*2, need to set the cluster for bottom value
					clusters[ilst + 1 - 1] = ifst_cluster;
				if (ilst != orig_ilst - 1) // second block (now in ifst) was 2*2, need to set the cluster for bottom value
					clusters[ifst + 1 - 1] = ilst_cluster;
				new_n = i;
			}
		}
		temp_n = new_n;
	}
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesSort);
	Benchmark::Instance().AddProperty(Benchmark::PropertyType::EigenvaluesSwapsCount, swaps_count);
	Benchmark::Instance().AddProperty(Benchmark::PropertyType::EigenvaluesTotalSwaps, total_swaps);
	return clusters;
}

template<typename ElementType>
std::vector<int> Polynomial<ElementType>::clusterSchurEigenvalues(Submatrix<ElementType>& schur, Submatrix<ElementType>& schur_co, Eigenvalues& eigenvalues, double tolerance, EigenvaluesClusterSortAlgorithm eigenvalues_sort_algorithm)
{
	int n = schur.dim();
	if (eigenvalues.size() != n)
	{
		eigenvalues.resize(n);
		for (int i = 0; i < n; i++)
			eigenvalues[i] = *(schur.data(i, i));
	}
	Benchmark::Instance().TraceCurrentStep(Benchmark::CheckpointType::EigenvaluesCluster);
	//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesInit);
	std::vector<int> clusters = clusterValues(eigenvalues.data(), n, tolerance);
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesCluster);
	std::vector<int> permutation = chooseClusterPermutation(clusters);
	Benchmark::Instance().AddProperty(Benchmark::PropertyType::ClustersCount, permutation.size());
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesPermute);

	int leading_dimension = n > 1 ? n : 1;

	if (eigenvalues_sort_algorithm == EigenvaluesClusterSortAlgorithm::SmartSwapSort)
	{
		SmartSwap::ClustersWithBlocks cluster_blocks;
		SmartSwap::ClustersWithBlocks ordered_blocks;
		buildClusterBlocks(schur, clusters, permutation, cluster_blocks, ordered_blocks);
		SmartSwap swap(cluster_blocks, ordered_blocks, permutation.size());
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesSwapInit);
		auto moves = swap.Solve();
		//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesSwapSolve);
		Benchmark::Instance().AddProperty(Benchmark::PropertyType::EigenvaluesSwapsCount, moves.size());
		Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesPermute);
		int total_swaps = 0;
		for (auto& move : moves)
		{
			int ifst = move.first + 1;
			int ilst = move.second + 1;
			//mexPrintf("Swapping %d <-> %d\n", ifst, ilst);
#if USE_LAPACKE
			int result = LAPACKE_FUNC(ElementType, trexc, MATRIX_ORDER, 'V', n, schur.data(), leading_dimension, schur_co.data(), leading_dimension, &ifst, &ilst);
#else
			int result = 0;
			std::vector<ElementType> work(n);
			dtrexc("V", &n, schur.data(), &leading_dimension, schur_co.data(), &leading_dimension, &ifst, &ilst, work.data(), &result);
#endif

			Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::dtrexc);
			total_swaps += ifst - ilst;
			if (result != 0)
				throw MatrixException("ctrexc", result);
			//mexPrintf("Swapped %d <-> %d\n", ifst, ilst);
			//for (int i = 0; i < schur.nrows(); i++)
			//	mexPrintf("%lf, ", *schur.data(i, i));
			//mexPrintf("\n");
		}
		std::vector<int> permuted_clusters = applyPermutation(clusters, permutation);
		Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::EigenvaluesPermute);
		Benchmark::Instance().AddProperty(Benchmark::PropertyType::EigenvaluesTotalSwaps, total_swaps);
		return permuted_clusters;
	}
	else // BubbleSort
	{
		return bubbleSortEigenvalues(schur, schur_co, n, clusters, permutation);
	}
}


//// Patterson-Stockmeyer Polynomial Evaluation
//void Polynomial::calculatePolynomialForCluster(Submatrix& A, std::vector<float>& coefficients, int from, int to)
//{
//	MKL_Complex8* Ap_start = A.data(from, from);
//	int submatrix_size = to - from;
//	Matrix qA = A.polynomial(coefficients, { from, from }, submatrix_size);
//	A.assign(qA, { from, from });
//}
//void Polynomial::calculatePolynomialForCluster(Submatrix& A, std::vector<float>& coefficients, MatrixIndex from, int submatrix_size, Matrix& output)
//{
//	MKL_Complex8* Ap_start = A.data(from);
//	Matrix qA = A.polynomial(coefficients, from, submatrix_size);
//	output.assign(qA, from);
//}

template<typename ElementType>
static void handle_block(Submatrix<ElementType>& T, int cluster_start, int cluster_end, int n, Polynomial<ElementType>::matrix_function f, Submatrix<ElementType>& F, std::vector<MatrixBlock<ElementType>>& blocks)
{

	int block_size = cluster_end - cluster_start;
	Submatrix<ElementType> block_data((ElementType*)T.data(cluster_start, cluster_start), block_size, block_size, n);
	MatrixBlock<ElementType> block(block_data, { cluster_start, cluster_start }, block_size);
	blocks.push_back(block);
	f(T, cluster_end - cluster_start, { cluster_start, cluster_start }, F); //calculate_polynom_for_cluster(T.data(cluster_start, cluster_start), cluster_end - cluster_start, , f, cluster_start, cluster_end);
}

template<typename ElementType>
std::vector<MatrixBlock<ElementType>> Polynomial<ElementType>::calculateDiagonalBlocks(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, Submatrix<ElementType>& F)
{
	int n = T.nrows();
	int cluster_start = 0;
	std::vector<MatrixBlock<ElementType>> blocks;
	for (unsigned int i = 0; i < clusters.size(); i++)
	{
		if (clusters[cluster_start] != clusters[i])
		{
			int cluster_end = i;
			handle_block(T, cluster_start, cluster_end, n, f, F, blocks);
			cluster_start = cluster_end;
		}
	}
	handle_block(T, cluster_start, clusters.size(), n, f, F, blocks);
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::PolynomialCalculation);
	return blocks;
}

// T is a triangular matrix
template<typename ElementType>
std::shared_ptr<Matrix<ElementType>> Polynomial<ElementType>::calculateBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function)
{
	int n = T.nrows();
	if (n != T.ncols())
		throw std::runtime_error("This algorithm is valid for square matrices only");

	std::shared_ptr<Matrix<ElementType>> F = std::make_shared<Matrix<ElementType>>(n, n);
	Submatrix<ElementType> sF(*F);
	////// Fii = f(Tii), i = 1: n
	std::vector<MatrixBlock<ElementType>> blocks = calculateDiagonalBlocks(T, clusters, f, sF);

	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::TotalParlettRecurrence,
		[sylvester_function, &blocks, &sF, &T]() {
		int nblocks = blocks.size();
		for (int j = 1; j < nblocks; j++)
		{
			for (int i = j - 1; i >= 0; i--)
			{
				// Solve for Fij the Sylvester equasion
				// TiiFij - FijTjj = FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
				// AX - XB = C
				Submatrix<ElementType>* Tii = blocks[i].matrix();
				Submatrix<ElementType>* Tjj = blocks[j].matrix();
				int i_row = blocks[i].from().first;
				//int i_col = blocks[i].from().second;
				//int j_row = blocks[j].from().first;
				int j_col = blocks[j].from().second;

				MatrixIndex ii = blocks[i].from();
				MatrixIndex ij = { i_row, j_col };
				MatrixIndex jj = blocks[j].from();

				int block_m = Tii->nrows(); // The order of A, and the number of rows in X and C (m >= 0).
				int block_n = Tjj->ncols(); // The order of B, and the number of columns in X and C(n >= 0).

				// Fij <- FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)

				Submatrix<ElementType>::multiply_ex(sF, T, sF, block_m, block_n, block_m, ii, ij, ij, false, false, Submatrix<ElementType>::MultiplyOption::Plus); // Fij += FiiTij
				Submatrix<ElementType>::multiply_ex(T, sF, sF, block_m, block_n, block_n, ij, jj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Minus); // Fij -= TijFjj

				// \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
				//Matrix FikTkj(block_m, block_n);
				//Matrix TikFkj(block_m, block_n);
				for (int k = i + 1; k < j; k++)
				{
					int k_row = blocks[k].from().first;
					int k_col = blocks[k].from().second;

					int ik_m = blocks[i].nrows();
					int ik_n = blocks[k].nrows();
					//int kj_m = blocks[k].nrows();
					int kj_n = blocks[j].nrows();

					MatrixIndex ik = { i_row, k_col };
					MatrixIndex kj = { k_row, j_col };

					Submatrix<ElementType>::multiply_ex(sF, T, sF, ik_m, kj_n, ik_n, ik, kj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Plus); // Fij += FikTkj
					Submatrix<ElementType>::multiply_ex(T, sF, sF, ik_m, kj_n, ik_n, ik, kj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Minus); // Fij -= TikFkj
				}

				Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::BlockParlettRecurrence);
				Submatrix<ElementType> sFij(sF.data(ij), block_m, block_n, sF.dim());
				Submatrix<ElementType>::sylvester(sylvester_function, *Tii, *Tjj, sFij);
				Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::SylvesterSolve);
			}
		}
	});
	return F;
}

template<typename ElementType>
void Polynomial<ElementType>::calculateBlockInRecurrence(Submatrix<ElementType>& T, int block_row, int block_col, std::vector<MatrixBlock<ElementType>>& blocks, SylvesterFunction sylvester_function, Submatrix<ElementType>& F)
{
	// i >= blocks.size() || j <= 0
	if (block_row == block_col) // Diagonal block
		return;
	
	//int i_row = blocks[i].from().first;
	//int j_col = blocks[j].from().second;

	//int left_block_row = i_row; // blocks[i].from().first;
	//int left_block_col = blocks[j - 1].from().second;
	//int bottom_block_row = blocks[i + i].from().first;
	//int bottom_block_col = j_col;

	// Left block
	/*cilk_spawn*/ calculateBlockInRecurrence(T, block_row, block_col - 1, blocks, sylvester_function, F);
	// Bottom block
	/*cilk_spawn*/ calculateBlockInRecurrence(T, block_row + 1, block_col, blocks, sylvester_function, F);

	//cilk_sync;

	// Solve for Fij the Sylvester equasion
	// TiiFij - FijTjj = FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
	// AX - XB = C
	int i = block_row;
	int j = block_col;
	Submatrix<ElementType>* Tii = blocks[i].matrix();
	Submatrix<ElementType>* Tjj = blocks[j].matrix();
	int i_row = blocks[i].from().first;
	int j_col = blocks[j].from().second;

	MatrixIndex ii = blocks[i].from();
	MatrixIndex ij = { i_row, j_col };
	MatrixIndex jj = blocks[j].from();

	int block_m = Tii->nrows(); // The order of A, and the number of rows in X and C (m >= 0).
	int block_n = Tjj->ncols(); // The order of B, and the number of columns in X and C(n >= 0).

								// Fij <- FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)

	Submatrix<ElementType>::multiply_ex(F, T, F, block_m, block_n, block_m, ii, ij, ij, false, false, true); // Fij += FiiTij
	Submatrix<ElementType>::multiply_ex(T, F, F, block_m, block_n, block_n, ij, jj, ij, false, false, false); // Fij -= TijFjj

																												// \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
																												//Matrix FikTkj(block_m, block_n);
																												//Matrix TikFkj(block_m, block_n);
	for (int k = i + 1; k < j; k++)
	{
		int k_row = blocks[k].from().first;
		int k_col = blocks[k].from().second;

		int ik_m = blocks[i].nrows();
		int ik_n = blocks[k].nrows();
		//int kj_m = blocks[k].nrows();
		int kj_n = blocks[j].nrows();

		MatrixIndex ik = { i_row, k_col };
		MatrixIndex kj = { k_row, j_col };

		Submatrix<ElementType>::multiply_ex(F, T, F, ik_m, kj_n, ik_n, ik, kj, ij, false, false, true); // Fij += FikTkj
		Submatrix<ElementType>::multiply_ex(T, F, F, ik_m, kj_n, ik_n, ik, kj, ij, false, false, false); // Fij -= TikFkj
	}

	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::BlockParlettRecurrence);
	Submatrix<ElementType> sFij(F.data(ij), block_m, block_n, F.dim());
	Submatrix<ElementType>::sylvester(sylvester_function, *Tii, *Tjj, sFij);
	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::SylvesterSolve);

}

// T is a triangular matrix
template<typename ElementType>
std::shared_ptr<Matrix<ElementType>> Polynomial<ElementType>::calculateParallelBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function)
{
	int n = T.nrows();
	if (n != T.ncols())
		throw std::runtime_error("This algorithm is valid for square matrices only");

	std::shared_ptr<Matrix<ElementType>> F = std::make_shared<Matrix<ElementType>>(n, n);
	Submatrix<ElementType> sF(*F);
	// Fii = f(Tii), i = 1:n
	std::vector<MatrixBlock<ElementType>> blocks = calculateDiagonalBlocks(T, clusters, f, sF);

	calculateBlockInRecurrence(T, 0, blocks.size() - 1, blocks, sylvester_function, sF);

	return F;
}

template<typename ElementType>
void calculateSuperdiagonalBlockParlettRecurrenceImpl(Submatrix<ElementType>& T, std::vector<MatrixBlock<ElementType>>& blocks, Submatrix<ElementType>& sF, SylvesterFunction sylvester_function)
{
	int nblocks = blocks.size();

	for (int d = 1; d < nblocks; d++) // For each superdiagonal
	{
		cilk_for(int i = 0; i < nblocks - d; i++)
		{
			int j = i + d;

			Submatrix<ElementType>* Tii = blocks[i].matrix();
			Submatrix<ElementType>* Tjj = blocks[j].matrix();
			int i_row = blocks[i].from().first;
			//int i_col = blocks[i].from().second;
			//int j_row = blocks[j].from().first;
			int j_col = blocks[j].from().second;

			MatrixIndex ii = blocks[i].from();
			MatrixIndex ij = { i_row, j_col };
			MatrixIndex jj = blocks[j].from();

			int block_m = Tii->nrows(); // The order of A, and the number of rows in X and C (m >= 0).
			int block_n = Tjj->ncols(); // The order of B, and the number of columns in X and C(n >= 0).

			Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::BlockParlettRecurrence,
				[i, j, d, i_row, j_col, &T, &ii, &ij, &jj, &blocks, &Tii, &Tjj, &sF, &block_m, &block_n]() {
				// Solve for Fij the Sylvester equasion
				// TiiFij - FijTjj = FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
				// AX - XB = C

											// Fij <- FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)

				Submatrix<ElementType>::multiply_ex(sF, T, sF, block_m, block_n, block_m, ii, ij, ij, false, false, Submatrix<ElementType>::MultiplyOption::Plus); // Fij += FiiTij
				Submatrix<ElementType>::multiply_ex(T, sF, sF, block_m, block_n, block_n, ij, jj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Minus); // Fij -= TijFjj

				// \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
				//Matrix FikTkj(block_m, block_n);
				//Matrix TikFkj(block_m, block_n);

				for (int k = i + 1; k < j; k++)
				{
					int k_row = blocks[k].from().first;
					int k_col = blocks[k].from().second;

					int ik_m = blocks[i].nrows();
					int ik_n = blocks[k].nrows();
					//int kj_m = blocks[k].nrows();
					int kj_n = blocks[j].nrows();

					MatrixIndex ik = { i_row, k_col };
					MatrixIndex kj = { k_row, j_col };

					Submatrix<ElementType>::multiply_ex(sF, T, sF, ik_m, kj_n, ik_n, ik, kj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Plus); // Fij += FikTkj
					Submatrix<ElementType>::multiply_ex(T, sF, sF, ik_m, kj_n, ik_n, ik, kj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Minus); // Fij -= TikFkj
				}
			});
			Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::SylvesterSolve,
				[sylvester_function, &sF, &ij, &Tii, &Tjj, &block_m, &block_n]() {
				Submatrix<ElementType> sFij(sF.data(ij), block_m, block_n, sF.dim());
				Submatrix<ElementType>::sylvester(sylvester_function, *Tii, *Tjj, sFij);
			});
		}
	}
}

// T is a triangular matrix
template<typename ElementType>
std::shared_ptr<Matrix<ElementType>> Polynomial<ElementType>::calculateSuperdiagonalBlockParlettRecurrence(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function)
{
	int n = T.nrows();
	if (n != T.ncols())
		throw std::runtime_error("This algorithm is valid for square matrices only");

	std::shared_ptr<Matrix<ElementType>> F = std::make_shared<Matrix<ElementType>>(n, n);
	Submatrix<ElementType> sF(*F);
	////// Fii = f(Tii), i = 1: n
	//Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::RecurrenceInit);
	std::vector<MatrixBlock<ElementType>> blocks = calculateDiagonalBlocks(T, clusters, f, sF);

	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::TotalParlettRecurrence,
		[sylvester_function, &blocks, &T, &sF]() {
		calculateSuperdiagonalBlockParlettRecurrenceImpl(T, blocks, sF, sylvester_function);
	});
	return F;
}

// T is a triangular matrix
template<typename ElementType>
std::shared_ptr<Matrix<ElementType>> Polynomial<ElementType>::calculateBlockedFunction(Submatrix<ElementType>& T, std::vector<int>& clusters, matrix_function f, SylvesterFunction sylvester_function)
{
	int n = T.nrows();
	if (n != T.ncols())
		throw std::runtime_error("This algorithm is valid for square matrices only");

	std::shared_ptr<Matrix<ElementType>> F = std::make_shared<Matrix<ElementType>>(n, n);
	Submatrix<ElementType> sF(*F);
	////// Fii = f(Tii), i = 1: n
	std::vector<MatrixBlock<ElementType>> blocks = calculateDiagonalBlocks(T, clusters, f, sF);

	Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::TotalParlettRecurrence,
		[sylvester_function, n, &blocks, &T, &sF]() {
		int nblocks = blocks.size();
		for (int i = 1; i < nblocks; i++)
		{
			// Solve for Fij the Sylvester equasion
			// TiiFij - FijTjj = FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)
			// AX - XB = C
			int left_block_size = blocks[i].from().first;
			Submatrix<ElementType> left_block_data((ElementType*)T.data(), left_block_size, left_block_size, n);
			MatrixBlock<ElementType> left_block(left_block_data, { 0, 0 }, left_block_size);
			MatrixBlock<ElementType>& current_block = blocks[i];

			Submatrix<ElementType>* Tii = left_block.matrix();
			Submatrix<ElementType>* Tjj = current_block.matrix();
			int i_row = left_block.from().first;
			//int i_col = blocks[i].from().second;
			//int j_row = blocks[j].from().first;
			int j_col = current_block.from().second;

			MatrixIndex ii = left_block.from();
			MatrixIndex ij = { i_row, j_col };
			MatrixIndex jj = current_block.from();

			int block_m = Tii->nrows(); // The order of A, and the number of rows in X and C (m >= 0).
			int block_n = Tjj->ncols(); // The order of B, and the number of columns in X and C(n >= 0).

			// Fij <- FiiTij - TijFjj + \sum_(j-1)^(k=i+1) (FikTkj - TikFkj)

			Submatrix<ElementType>::multiply_ex(sF, T, sF, block_m, block_n, block_m, ii, ij, ij, false, false, Submatrix<ElementType>::MultiplyOption::Plus); // Fij += FiiTij
			Submatrix<ElementType>::multiply_ex(T, sF, sF, block_m, block_n, block_n, ij, jj, ij, false, false, Submatrix<ElementType>::MultiplyOption::Minus); // Fij -= TijFjj

			// \sum_(j-1)^(k=i+1) (FikTkj - TikFkj) - irrelevant here

			Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::BlockParlettRecurrence);
			MatrixIndex iblock = { 0, j_col };
			Submatrix<ElementType> sFij(sF.data(iblock), block_m, block_n, sF.dim());
			Submatrix<ElementType>::sylvester(sylvester_function, *Tii, *Tjj, sFij);
			Benchmark::Instance().Checkpoint(Benchmark::CheckpointType::SylvesterSolve);
		}
	});
	return F;
}
