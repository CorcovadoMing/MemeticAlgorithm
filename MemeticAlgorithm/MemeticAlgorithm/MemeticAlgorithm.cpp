#include "MemeticAlgorithm.h"
#include "RandomRange.h"
#include <iostream>
#include <fstream>
#include <algorithm>

MemeticAlgorithm::MemeticAlgorithm(const int population_size, const double crossover_rate, const double mutation_rate, const int localsearch_looptimes, const std::string& filename)
	: population_size_(population_size), crossover_rate_(crossover_rate), mutation_rate_(mutation_rate), localsearch_looptimes_(localsearch_looptimes), filename_(filename)
{
	readfile();
	if (jobs_ == 0 || machines_ == 0)
	{
		std::cout << "[Error] jobs = 0, machines = 0" << std::endl;
	}
	population_ = Population(population_size_, Chromosome(jobs_, 0));

	initialize_.push_back(std::mem_fn(&MemeticAlgorithm::randomInitialize));
	initialize_.push_back(std::mem_fn(&MemeticAlgorithm::heuristicInitialize));

	//crossover_ 

	mutation_.push_back(std::mem_fn(&MemeticAlgorithm::randomSwap));

	localsearch_.push_back(std::mem_fn(&MemeticAlgorithm::II));
	localsearch_.push_back(std::mem_fn(&MemeticAlgorithm::SA));
	localsearch_.push_back(std::mem_fn(&MemeticAlgorithm::TS));
}

void MemeticAlgorithm::run()
{
	std::cout << "=== Testing ===" << std::endl;
	for (std::size_t i = 0; i < initialize_.size(); i += 1)
	{
		initialize_[i](this);
	}

	for (std::size_t i = 0; i < mutation_.size(); i += 1)
	{
		mutation_[i](this, population_[0]);
	}

	for (std::size_t i = 0; i < localsearch_.size(); i += 1)
	{
		localsearch_[i](this, population_[0]);
	}
}

#pragma region Initialization

void MemeticAlgorithm::randomInitialize()
{
	std::cout << "RandomInitialize" << std::endl;
	// TODO: Implemet a randome initialization
}

void MemeticAlgorithm::heuristicInitialize()
{
	std::cout << "HeuristicInitialize" << std::endl;
	// TODO: Implement a heuristic initialization
}

#pragma endregion

#pragma region Crossover
#pragma endregion

#pragma region Mutation

void MemeticAlgorithm::randomSwap(Chromosome &chromosome)
{
	// Assign to Ming rf37535@gmail.com
	const int firstElement = RandomRange::random<int>(0, chromosome.size() - 1);
	const int secondElement = RandomRange::random<int>(0, chromosome.size() - 1);
	std::swap(chromosome[firstElement], chromosome[secondElement]);
}

#pragma endregion

#pragma region LocalSearch

const Chromosome MemeticAlgorithm::II(const Chromosome &chromosome)
{
	Chromosome result(chromosome);
	// TODO: Implement a II and improve result
	// !!Implement!!: use randomSwap and fitness
	// please use RandomRange::random<int>(min, max) or RandomRange::random<double>(min, max) to generate random number [min, max]
	return result;
}

const Chromosome MemeticAlgorithm::SA(const Chromosome &chromosome)
{
	Chromosome result(chromosome);
	// TODO: Implement a SA and improve result
	// !!Implement!!: use randomSwap and fitness
	// please use RandomRange::random<int>(min, max) or RandomRange::random<double>(min, max) to generate random number [min, max]
	return result;
}

const Chromosome MemeticAlgorithm::TS(const Chromosome &chromosome)
{
	// Assign to Ming rf37535@gmail.com
	Chromosome result(chromosome);
	return result;
}

#pragma endregion

#pragma region Helper-Functions

void MemeticAlgorithm::readfile()
{
	std::string useless;
	std::ifstream fin;
	fin.open(filename_.c_str());

	fin >> jobs_ >> machines_ >> useless;

	Matrix matrix(machines_, std::vector<int>(jobs_, 0));

	for (std::size_t i = 0; i < machines_; i += 1)
	{
		for (std::size_t j = 0; j < jobs_; j += 1)
		{
			fin >> matrix[i][j];
		}
	}
	matrix_ = matrix;
}

const int MemeticAlgorithm::fitness(const Chromosome& chromosome)
{
	std::vector<int> timespan(matrix_.size() + 1, 0);
	for (std::size_t i = 0; i < matrix_[0].size(); i += 1)
	{
		for (std::size_t j = 1; j <= matrix_.size(); j += 1)
		{
			timespan[j] = std::max(timespan[j - 1], timespan[j]) + matrix_[j - 1][chromosome[i]];
		}
	}
	return timespan[matrix_.size()];
}

#pragma endregion