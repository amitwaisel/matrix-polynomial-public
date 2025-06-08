#pragma once
#include <ctime>
#include <memory>
#include <map>
#include <chrono>  // for high_resolution_clock
#include <string>
#include "defs.h"

#include <Windows.h>

class Benchmark
{
public:
#if USE_CILK
	typedef cilk::reducer<cilk::op_add<std::chrono::nanoseconds>> timing;
	typedef cilk::reducer<cilk::op_add<int>> counter;
	typedef cilk::reducer<cilk::op_add<double>> cpu_counter;
#define getvalue(x) x->get_value()
#else
	typedef std::chrono::nanoseconds timing;
	typedef int counter;
	typedef double cpu_counter;
#define getvalue(x) x
#endif

	enum CheckpointType
	{
		SchurDecomposition = 0,
		//EigenvaluesInit,
		EigenvaluesCluster,
		EigenvaluesPermute,
		//EigenvaluesSwapInit,
		//EigenvaluesSwapSolve,
		EigenvaluesSort,
		dtrexc,
		// RecurrenceInit,
		PolynomialCalculation,
		// PolynomialCalculationInit,
		// PolynomialCalculationAps,
		// PolynomialCalculationInnerInit,
		// PolynomialCalculationInner,
		// PolynomialCalculationHorner,
		// PolynomialAssign,
		BlockParlettRecurrence,
		SylvesterSolve,
		TotalParlettRecurrence,
		FinalMultiplication,
		Checkpoint_MAX
	};

	static char* checkpoint_names[];

	enum PropertyType
	{
		ClustersCount = 0,
		SingleBlocksCount,
		DoubleBlocksCount,
		EigenvaluesSwapsCount,
		EigenvaluesTotalSwaps,
		Property_MAX
	};

private:
	Benchmark() = default;
	~Benchmark() = default;

public:
	static Benchmark& Instance()
	{
		static Benchmark instance; // Guaranteed to be destroyed. Instantiated on first use.
		return instance;
	}
	void Start();
	void Checkpoint(CheckpointType type);
	void Checkpoint(CheckpointType type, std::function<void(void)> func);
	void AddProperty(PropertyType type, double value);
	void TraceCurrentStep(CheckpointType type);

	std::chrono::nanoseconds Get(CheckpointType type);
	double Get(PropertyType type);
	int GetCount(CheckpointType type);
	double GetCpu(CheckpointType type);

private:
	void init();
	double getCurrentCpu();

private:
	std::chrono::high_resolution_clock::time_point _start;
	std::map<CheckpointType, std::shared_ptr<timing>> _timings;
	std::map<CheckpointType, std::shared_ptr<counter>> _counts;
	std::map<CheckpointType, std::shared_ptr<cpu_counter>> _cpu;
	std::map<PropertyType, double> _properties;

	ULARGE_INTEGER lastCPU, lastSysCPU, lastUserCPU;
	int numProcessors;
};

