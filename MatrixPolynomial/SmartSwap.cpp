#include "SmartSwap.h"
#include <functional>
#include <algorithm>

SmartSwap::SmartSwap(const ClustersWithBlocks& clusters, const ClustersWithBlocks& order, int k) :
	_clusters(clusters.clusters), _clusters_blocks(clusters.blocks), _order(order.clusters), _order_blocks(order.blocks), _k(k), _state(k)
{
	if (_clusters.size() != _order.size() || _clusters_blocks.size() != _order_blocks.size())
		throw std::runtime_error("Clusters and Order must have the same size");

	int absolute_index = 0;
	for (size_t i = 0; i < _clusters.size(); i++)
	{
		int current = _clusters[i];
		int current_block_size = _clusters_blocks[i];
		int destination = _order[i];
		int destination_block_size = _order_blocks[i];
		if (current != destination || current_block_size != destination_block_size)
		{
			//_state.Add(current, destination);
			//_bad_locations.emplace(std::make_pair(i, BadLocation(current, current_block_size, destination, destination_block_size, i, absolute_index));
			_bad_locations[i] = BadLocation(current, current_block_size, destination, destination_block_size, i, absolute_index);
		}
		else
			_good_locations[i] = BadLocation(current, current_block_size, destination, destination_block_size, i, absolute_index);

		absolute_index += current_block_size;
	}
}

SmartSwap::Moves SmartSwap::Solve()
{
	SmartSwap::Moves moves;
	
	// Swap couples
	//auto swaps = _state.Solve();
	//for (auto& swap_couple : swaps)
	//{
	//	auto c1 = GetCandidate(swap_couple.first);
	//	auto c2 = GetCandidate(swap_couple.second);
	//	_bad_locations.erase(c1);
	//	_bad_locations.erase(c2);
	//	//_bad_locations.erase(std::find(_bad_locations.begin(), _bad_locations.end(), c1));
	//	//_bad_locations.erase(std::find(_bad_locations.begin(), _bad_locations.end(), c2));
	//	moves.push_back(SwapState::SwapMove(c1.index(), c2.index()));
	//	if (_order[c2.index()] != c1.current())
	//		_bad_locations.emplace(c1.current(), _order[c2.index()], c2.index());
	//	if (_order[c1.index()] != c2.current())
	//		_bad_locations.emplace(c2.current(), _order[c1.index()], c1.index());
	//}

	while (!_bad_locations.empty())
	{
		// sets are ordered, the first element is the smallest (minimal index, by operator<)
		auto bad_location = _bad_locations.begin()->second;
		auto c = GetCandidate([&bad_location](const BadLocation& l) { return l.current() == bad_location.destination() && l.current_block_size() == bad_location.destination_block_size(); });
		SwapState::SwapMove move(c.absolute_index(), bad_location.absolute_index());
		DoSwap(bad_location, c);
		moves.push_back(move);
	}
	return moves;
	/*
		while len(self.bad_locs) > 0:
			bad_loc = min([l for l in self.bad_locs], key = lambda l : l.index)
			c = self.get_candidate_ex(lambda l : l.src_c == bad_loc.final_c and l.src_size == bad_loc.final_size)
			self.do_swap(bad_loc, c)
			moves.append((bad_loc.real_index, c.real_index))
	*/
}

void SmartSwap::DoSwap(const BadLocation& left, const BadLocation& right) // right goes to left position
{
	_bad_locations.erase(right.index());

	int new_absolute_index = left.absolute_index();
	std::map<int, BadLocation> new_locations;
	int shift_size = right.current_block_size();
	for (int i = left.index(); i < right.index(); i++)
	{
		BadLocation loc; // (_clusters[i], _clusters_blocks[i], _order[i], _order_blocks[i], i, new_absolute_index);
		if (_bad_locations.find(i) != _bad_locations.end())
		{
			loc = _bad_locations[i];
			_bad_locations.erase(i);
		}
		else
		{
			loc = _good_locations[i];
			//_good_locations.erase(loc);
		}
		int new_index = i + 1;
		new_absolute_index += shift_size;
		shift_size = loc.current_block_size();
		if (_order[new_index] != loc.current() || _order_blocks[new_index] != loc.current_block_size())
			new_locations[new_index] = BadLocation(loc.current(), loc.current_block_size(), _order[new_index], _order_blocks[new_index], new_index, new_absolute_index);
		else
			_good_locations[new_index] = BadLocation(loc.current(), loc.current_block_size(), _order[new_index], _order_blocks[new_index], new_index, new_absolute_index);
	}
	_bad_locations.insert(new_locations.begin(), new_locations.end());
	/*
	def do_swap(self, c1, c2): # c2 goes left to its right location
        self.bad_locs.remove(c2)
        
        shift_size = c2.src_size
        new_absolute_index = c1.real_index
        new_locs = set()
        for i in range(c1.index, c2.index):
            i_predicate = lambda l: l.index == i
            try:
                c = self.get_candidate_ex(i_predicate)
                self.bad_locs.remove(c)
            except: # No bad location for this index
                c = self.get_good_candidate(i_predicate)
                self.good_locs.remove(c)
            new_index = i + 1
            new_absolute_index += shift_size #c.src_size
            shift_size = c.src_size # Next block will be after this one...
            final_c, final_size = self.final[new_index]
            new_location = BadLoc(c.src_c, c.src_size, final_c, final_size, new_index, new_absolute_index)
            if final_c != c.src_c or final_size != c.src_size:
                new_locs.add(new_location)
            else:
                self.good_locs.add(new_location)
        self.bad_locs.update(new_locs)
	*/
}

bool SmartSwap::TryGetCandidate(std::function<bool(const BadLocation&)> pred, BadLocation* output)
{
	for (auto& loc : _bad_locations)
	{
		if (pred(loc.second))
		{
			*output = loc.second;
			return true;
		}
	}
	return false;
}
BadLocation SmartSwap::GetCandidate(std::function<bool(const BadLocation&)> pred)
{
	BadLocation loc;
	if (!TryGetCandidate(pred, &loc))
		throw std::runtime_error("Failed finding candidates for swap");
	return loc;
}
BadLocation SmartSwap::GetGoodCandidate(std::function<bool(const BadLocation&)> pred)
{
	for (auto& loc : _good_locations)
	{
		if (pred(loc.second))
		{
			return loc.second;
		}
	}
	throw std::runtime_error("Failed finding good candidates for swap");
}

void SwapState::Add(int current, int destination)
{
	if (current != destination)
	{
		int* counter = _state.data(current, destination);
		(*counter)++;
	}
	// (*(_state.data(current, destination)))++;
}

std::vector<SwapState::SwapCouple> SwapState::Solve()
{
	std::vector<SwapCouple> swaps;
	for (int i = 0; i < _k; i++)
	{
		for (int j = i + 1; j < _k; j++)
		{
			if (*(_state.data(i, j)) > 0 && *(_state.data(j, i)) > 0)
			{
				(*(_state.data(i, j)))--;
				(*(_state.data(j, i)))--;
				swaps.push_back(SwapCouple(SwapMove(i, j), SwapMove(j, i)));
			}
		}
	}
	return swaps;
}
