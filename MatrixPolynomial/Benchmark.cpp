#include "Benchmark.h"
#include <numeric>
#include <algorithm>
#include <vector>

#ifdef USE_MEX
#include "mex.h"
#else
#define mexPrintf prinf
#endif

char* Benchmark::checkpoint_names[] = {
	   "Schur decomposition",
	   "Cluster eigenvalues",
	   "Permute eigenvalues",
	   "Sort eigenvalues",
	   "Eigenvalue swap",
	   "Polynomial calculation",
	   "Block parlett recurrence",
	   "Sylvester solve",
	   "Parlett recurrence",
	   "Final multiplication"
};

void Benchmark::Start()
{
	_timings.clear();
	_counts.clear();
	_cpu.clear();
	_properties.clear();
	for (int i = 0; i < CheckpointType::Checkpoint_MAX; i++)
	{
		CheckpointType cp = (CheckpointType)i;
		_timings[cp] = std::make_shared<timing>();
		_cpu[cp] = std::make_shared<cpu_counter>();
		_counts[cp] = std::make_shared<counter>(0);
	}
	for (int i = 0; i < PropertyType::Property_MAX; i++)
	{
		PropertyType pt = (PropertyType)i;
		_properties[pt] = 0;
	}
	init();

	_start = std::chrono::high_resolution_clock::now();
}

void Benchmark::init() {
	SYSTEM_INFO sysInfo;
	FILETIME ftime, fsys, fuser;

	GetSystemInfo(&sysInfo);
	numProcessors = sysInfo.dwNumberOfProcessors;

	GetSystemTimeAsFileTime(&ftime);
	memcpy(&lastCPU, &ftime, sizeof(FILETIME));

	GetProcessTimes(GetCurrentProcess(), &ftime, &ftime, &fsys, &fuser);
	memcpy(&lastSysCPU, &fsys, sizeof(FILETIME));
	memcpy(&lastUserCPU, &fuser, sizeof(FILETIME));
}

double Benchmark::getCurrentCpu() {
	FILETIME ftime, fsys, fuser;
	ULARGE_INTEGER now, sys, user;
	double percent;

	GetSystemTimeAsFileTime(&ftime);
	memcpy(&now, &ftime, sizeof(FILETIME));

	GetProcessTimes(GetCurrentProcess(), &ftime, &ftime, &fsys, &fuser);
	memcpy(&sys, &fsys, sizeof(FILETIME));
	memcpy(&user, &fuser, sizeof(FILETIME));
	percent = (sys.QuadPart - lastSysCPU.QuadPart) +
		(user.QuadPart - lastUserCPU.QuadPart);
	percent /= (now.QuadPart - lastCPU.QuadPart);
	percent /= numProcessors;
	lastCPU = now;
	lastUserCPU = user;
	lastSysCPU = sys;

	return percent * 100;
}

std::time_t GetCurrentTimeStamp()
{
	return std::time(nullptr);
}

int64_t GetCurrentMsTimeStamp()
{
	using namespace std::chrono;
	milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
	return ms.count();
}

std::string GetUtcTimeString(std::time_t time)
{
	char time_buffer[21] = { 0 }; // 2017-12-11T10:49:28Z
	if (strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&time)) == 0)
		return "1900-01-01T00:00:00Z";
	return time_buffer;
}

std::string GetCurrentTimeStampString()
{
	return GetUtcTimeString(GetCurrentTimeStamp());
}


static std::vector<Benchmark::CheckpointType> SKIP_CHECKPOINTS = { Benchmark::CheckpointType::BlockParlettRecurrence, Benchmark::CheckpointType::SylvesterSolve };
void Benchmark::TraceCurrentStep(Benchmark::CheckpointType type)
{
	std::string now = GetCurrentTimeStampString();
	mexPrintf("Starting step %s at %s\n", checkpoint_names[type], now.c_str());
}

void Benchmark::Checkpoint(CheckpointType type, std::function<void(void)> func)
{
	if (std::find(SKIP_CHECKPOINTS.begin(), SKIP_CHECKPOINTS.end(), type) != SKIP_CHECKPOINTS.end())
	{
		func();
		return; // Skip checkpoint;
	}

	//double cpu_before = getCurrentCpu();
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	func();
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	//double cpu = getCurrentCpu();

	std::chrono::nanoseconds diff = end - start;
#if USE_CILK
	(*_timings[type])->operator+=(diff);
	//(*_cpu[type])->operator+=(cpu);
	(*_counts[type])->operator++();
#else
	*_timings[type] += diff;
	//*_cpu[type] += cpu ;
	*_counts[type]++;
#endif
	_start = std::chrono::high_resolution_clock::now();
}

void Benchmark::Checkpoint(CheckpointType type)
{
	if (std::find(SKIP_CHECKPOINTS.begin(), SKIP_CHECKPOINTS.end(), type) != SKIP_CHECKPOINTS.end())
		return; // Skip checkpoint;
	std::chrono::high_resolution_clock::time_point start = _start;
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::nanoseconds diff = end - start;
	//double cpu = getCurrentCpu();

	//if (diff.count() < 0)
	//	__debugbreak();
	// mexPrintf("Diif for %d is %lld, start %lld end %lld", type, diff.count(), _start.time_since_epoch().count(), end.time_since_epoch().count());

	//std::lock_guard<std::mutex> lock(_mutex);
	// Doesn't take time from the measurement
#if USE_CILK
	(*_timings[type])->operator+=(diff);
	//(*_cpu[type])->operator+=(cpu);
	(*_counts[type])->operator++();
#else
	*_timings[type] += diff;
	//*_cpu[type] += cpu;
	*_counts[type]++;
#endif
	_start = std::chrono::high_resolution_clock::now();
}

void Benchmark::AddProperty(PropertyType type, double value)
{
	_properties[type] = value;
}

std::chrono::nanoseconds Benchmark::Get(CheckpointType type)
{
	return getvalue(_timings[type]);
}

double Benchmark::GetCpu(CheckpointType type)
{
	return getvalue(_cpu[type]);
}

double Benchmark::Get(PropertyType type)
{
	return _properties[type];
}

int Benchmark::GetCount(CheckpointType type)
{
	return getvalue(_counts[type]);
}
