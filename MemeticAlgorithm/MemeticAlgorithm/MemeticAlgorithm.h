#pragma once

#include <vector>
#include <string>
#include <functional>

typedef std::vector<int> Chromosome;
typedef std::vector<Chromosome> Population;
typedef std::vector< std::vector<int> > Matrix;

#pragma region Functions-Definition
typedef std::function<void()> InitializeOperation;
typedef std::vector<InitializeOperation> Initialize;

typedef std::function<const Chromosome(const Chromosome &, const Chromosome &)> CrossoverOperation;
typedef std::vector<CrossoverOperation> Crossover;

typedef std::function<void(Chromosome &)> MutationOperation;
typedef std::vector<MutationOperation> Mutation;

typedef std::function<const Chromosome(const Chromosome &)> LocalSearchOperation;
typedef std::vector<LocalSearchOperation> LocalSearch;
#pragma endregion

class MemeticAlgorithm
{
public:
	MemeticAlgorithm(const int, const double, const double, const std::string &);
	void run();
private:
	Matrix matrix_;
	void readfile();
	const int fitness(const Chromosome &);

	Initialize initialize_;
	Crossover crossover_;
	Mutation mutation_;
	LocalSearch localsearch_;

	static void randomInitialize();
	static void heuristicInitialize();

	static void randomSwap(Chromosome &);

	static Chromosome II(const Chromosome &);
	static Chromosome SA(const Chromosome &);
	static Chromosome TS(const Chromosome &);
	//void matingSelect();
	//void environmentSelect();

	Population population_;
	const double crossover_rate_, mutation_rate_;
	const int population_size_;
	
	const std::string filename_;
	int jobs = 0, machines = 0;
};

