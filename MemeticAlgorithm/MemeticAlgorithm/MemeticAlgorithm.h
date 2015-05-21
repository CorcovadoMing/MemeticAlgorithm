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

    typedef std::function<void(MemeticAlgorithm*, const Chromosome &, const Chromosome &)> CrossoverOperation;
    typedef std::vector<CrossoverOperation> Crossover;

    typedef std::function<void(MemeticAlgorithm*, Chromosome &)> MutationOperation;
    typedef std::vector<MutationOperation> Mutation;

    typedef std::function<const Chromosome(MemeticAlgorithm*, const Chromosome &)> LocalSearchOperation;
    typedef std::vector<LocalSearchOperation> LocalSearch;

	typedef std::function<const int(MemeticAlgorithm*, Chromosome &, const int)> ApplyLocalSearchOperation;
	typedef std::vector<ApplyLocalSearchOperation> ApplyLocalSearch;

public:
    MemeticAlgorithm(const int, const double, const double, const int, const std::string &);
    void run();

private:
	Population population_, offspring_;
    Matrix matrix_;
    Initialize initialize_;
    Crossover crossover_;
    Mutation mutation_;
    LocalSearch localSearch_;
	ApplyLocalSearch applyLocalSearch_;

    // Helper
    void readfile_();
    const int fitness_(const Chromosome &);

    // Initialize
    void randomInitialize();
    void heuristicInitialize();

    // Mating Selection
	void tournament();

    // Crossover
	void OX(const Chromosome &, const Chromosome &);
	void LOX(const Chromosome &, const Chromosome &);
	void PMX(const Chromosome &, const Chromosome &);
	void CX(const Chromosome &, const Chromosome &);

    // Mutation
	void insertion(Chromosome &);
	void swap(Chromosome &);
	void inverse(Chromosome &);

    // Local Search
    const Chromosome II(const Chromosome &);
    const Chromosome SA(const Chromosome &);
    const Chromosome TS(const Chromosome &);

	// Apply Local Search
	const int applyLocalSearchByLamarckian(Chromosome &, const int);
	const int applyLocalSearchByBaldwinian(Chromosome &, const int);

    //void environmentSelect();
	void tophalf();
	void generationModel();

    const double crossover_rate_;
    const double mutation_rate_;
    const int population_size_;
    const int localsearch_looptimes_;
    const std::string filename_;
    unsigned jobs_ = 0;
    unsigned machines_ = 0;

	std::vector<std::size_t> parent_;
	std::vector<int> fitness_table_;
};

