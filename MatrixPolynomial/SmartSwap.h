#pragma once
#include "Matrix.h"
#include <set>
#include <map>
#include <functional>
#include <memory>

class BadLocation
{
public:
	BadLocation() = default;
	BadLocation(int current, int current_block_size, int destination, int destination_block_size, int index, int absolute_index) :
		_current_cluster(current), _current_block_size(current_block_size),
		_destination_cluster(destination), _destination_block_size(destination_block_size),
		_index(index), _absolute_index(absolute_index)
	{ }
	~BadLocation() = default;

	bool operator==(const BadLocation& other) const
	{
		return _current_cluster == other._current_cluster &&
			_current_block_size == other._current_block_size &&
			_destination_cluster == other._destination_cluster &&
			_destination_block_size == other._destination_block_size &&
			_index == other._index &&
			_absolute_index == other._absolute_index;
	}
	bool operator< (const BadLocation& other) const
	{
		return _index < other._index;
	}

	int current() const { return _current_cluster; }
	int current_block_size() const { return _current_block_size; }
	int destination() const { return _destination_cluster; }
	int destination_block_size() const { return _destination_block_size; }
	int index() const { return _index; }
	int absolute_index() const { return _absolute_index; }
private:
	int _current_cluster;
	int _current_block_size;
	int _destination_cluster;
	int _destination_block_size;
	int _index;
	int _absolute_index;
};

class SwapState
{
public:
	typedef std::pair<int, int> SwapMove;
	typedef std::pair<SwapMove, SwapMove> SwapCouple;

	SwapState(int k) : _k(k), _state(k, k) { }
	~SwapState() = default;

	void Add(int current, int destination);
	std::vector<SwapCouple> Solve();


private:
	int _k;
	Matrix<int> _state;
};

class SmartSwap
{
public:
	typedef std::vector<int> Clusters;
	typedef std::vector<SwapState::SwapMove> Moves;
	typedef struct
	{
		Clusters clusters;
		Clusters blocks;
	} ClustersWithBlocks;
	SmartSwap(const ClustersWithBlocks& clusters, const ClustersWithBlocks& order, int k);
	~SmartSwap() = default;

	Moves Solve();

private:
	BadLocation GetCandidate(int current, int current_block_size);
	BadLocation GetCandidate(std::function<bool(const BadLocation&)> pred);
	BadLocation GetGoodCandidate(std::function<bool(const BadLocation&)> pred);
	bool TryGetCandidate(std::function<bool(const BadLocation&)> pred, BadLocation* output);
	void DoSwap(const BadLocation& c1, const BadLocation& c2);
	//BadLocation GetCandidate(const SwapState::SwapMove& move);

private:
	Clusters _clusters;
	Clusters _clusters_blocks;
	Clusters _order;
	Clusters _order_blocks;
	int _k;
	SwapState _state;
	//std::set<BadLocation, decltype(BadLocation::operator<)> _bad_locations;
	std::map<int, BadLocation> _bad_locations;
	std::map<int, BadLocation> _good_locations;
	//std::vector<BadLocation> _bad_locations;
};

