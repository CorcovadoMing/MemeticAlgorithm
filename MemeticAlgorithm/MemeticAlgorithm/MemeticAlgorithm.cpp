#include "MemeticAlgorithm.h"
#include "RandomRange.h"
#include <iostream>
#include <fstream>
#include <algorithm>

MemeticAlgorithm::MemeticAlgorithm(const int population_size, const double crossover_rate, const double mutation_rate, const std::string& filename)
	: population_size_(population_size), crossover_rate_(crossover_rate), mutation_rate_(mutation_rate), filename_(filename)
{
	readfile();
	if (jobs_ == 0 || machines_ == 0)
	{
		std::cout << "[Error] jobs = 0, machines = 0" << std::endl;
	}
	population_ = Population(population_size_, Chromosome(jobs_, 0));

	initialize_.push_back(randomInitialize);
	initialize_.push_back(heuristicInitialize);

	//crossover_ 

	mutation_.push_back(randomSwap);

	localsearch_.push_back(II);
	localsearch_.push_back(SA);
	localsearch_.push_back(TS);
}

void MemeticAlgorithm::run()
{
	std::cout << "=== Testing ===" << std::endl;
	for (std::size_t i = 0; i < initialize_.size(); i += 1)
	{
		initialize_[i]();
	}

	for (std::size_t i = 0; i < mutation_.size(); i += 1)
	{
		mutation_[i](population_[0]);
	}

	for (std::size_t i = 0; i < localsearch_.size(); i += 1)
	{
		localsearch_[i](population_[0]);
	}
}

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

void MemeticAlgorithm::randomSwap(Chromosome &chromosome)
{
	// TODO: Implement a function that randomly swap two elements in chromosome 
}

Chromosome MemeticAlgorithm::II(const Chromosome &chromosome)
{
	Chromosome result(chromosome);
	// TODO: Implement a II and improve result
	// !!Implement!!: use randomSwap and fitness
	// please use RandomRange::random<int>(min, max) or RandomRange::random<double>(min, max) to generate random number [min, max)
	return result;
}

Chromosome MemeticAlgorithm::SA(const Chromosome &chromosome)
{
	Chromosome result(chromosome);
	// TODO: Implement a SA and improve result
	// !!Implement!!: use randomSwap and fitness
	// please use RandomRange::random<int>(min, max) or RandomRange::random<double>(min, max) to generate random number [min, max)
	return result;
}

Chromosome MemeticAlgorithm::TS(const Chromosome &chromosome)
{
	Chromosome result(chromosome);
	// TODO: Implement a TS and improve result
	// !!Implement!!: use randomSwap and fitness
	// please use RandomRange::random<int>(min, max) or RandomRange::random<double>(min, max) to generate random number [min, max)
	return result;
}

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