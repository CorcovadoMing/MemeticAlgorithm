#pragma once

#include <vector>
#include <string>
#include <functional>

typedef std::vector<int> Chromosome;
typedef std::vector<Chromosome> Population;
typedef std::vector< std::vector<int> > Matrix;

class MemeticAlgorithm
{
	typedef std::function<void(MemeticAlgorithm*)> InitializeOperation;
	typedef std::vector<InitializeOperation> Initialize;

	typedef std::function<const Chromosome(MemeticAlgorithm*, const Chromosome &, const Chromosome &)> CrossoverOperation;
	typedef std::vector<CrossoverOperation> Crossover;

	typedef std::function<void(MemeticAlgorithm*, Chromosome &)> MutationOperation;
	typedef std::vector<MutationOperation> Mutation;

	typedef std::function<const Chromosome(MemeticAlgorithm*, const Chromosome &)> LocalSearchOperation;
	typedef std::vector<LocalSearchOperation> LocalSearch;

public:
	MemeticAlgorithm(const int, const double, const double, const int, const std::string &);
	void run();

private:
	Population population_;
	Matrix matrix_;
	Initialize initialize_;
	Crossover crossover_;
	Mutation mutation_;
	LocalSearch localsearch_;

	void readfile();
	const int fitness(const Chromosome &);
	void randomSwap(Chromosome &);

	void randomInitialize();
	void heuristicInitialize();


	const Chromosome II(const Chromosome &);
	const Chromosome SA(const Chromosome &);
	const Chromosome TS(const Chromosome &);
	//void matingSelect();
	//void environmentSelect();

	const double crossover_rate_;
	const double mutation_rate_;
	const int population_size_;
	const int localsearch_looptimes_;
	const std::string filename_;
	unsigned jobs_ = 0;
	unsigned machines_ = 0;
};

