#include "MemeticAlgorithm.h"
#include "RandomRange.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>

#pragma region Constructor

MemeticAlgorithm::MemeticAlgorithm(const int population_size, const double crossover_rate, const double mutation_rate, const int localsearch_looptimes, const std::string& filename)
    : population_size_(population_size), crossover_rate_(crossover_rate), mutation_rate_(mutation_rate), localsearch_looptimes_(localsearch_looptimes), filename_(filename)
{
    readfile_();
    if (jobs_ == 0 || machines_ == 0)
    {
        std::cout << "[Error] jobs = 0, machines = 0" << std::endl;
    }
    population_ = Population(population_size_, Chromosome(jobs_, 0));

    initialize_.push_back(std::mem_fn(&MemeticAlgorithm::randomInitialize));
    initialize_.push_back(std::mem_fn(&MemeticAlgorithm::heuristicInitialize));

    crossover_.push_back(std::mem_fn(&MemeticAlgorithm::OX));
    crossover_.push_back(std::mem_fn(&MemeticAlgorithm::LOX));
    crossover_.push_back(std::mem_fn(&MemeticAlgorithm::PMX));
    crossover_.push_back(std::mem_fn(&MemeticAlgorithm::CX));

    //mutation_

	localSearch_.push_back(std::mem_fn(&MemeticAlgorithm::II));
	localSearch_.push_back(std::mem_fn(&MemeticAlgorithm::SA));
	localSearch_.push_back(std::mem_fn(&MemeticAlgorithm::TS));

	applyLocalSearch_.push_back(std::mem_fn(&MemeticAlgorithm::applyLocalSearchByLamarckian));
	applyLocalSearch_.push_back(std::mem_fn(&MemeticAlgorithm::applyLocalSearchByBaldwinian));
}

#pragma endregion

#pragma region MainLogic

void MemeticAlgorithm::run()
{
    std::cout << "=== Testing ===" << std::endl;
    for (std::size_t i = 0; i < initialize_.size(); i += 1)
    {
        initialize_[i](this);
    }

    for (std::size_t i = 0; i < crossover_.size(); i += 1)
    {
        crossover_[i](this, population_[0], population_[1]);
    }

    for (std::size_t i = 0; i < mutation_.size(); i += 1)
    {
        mutation_[i](this, population_[0]);
    }

    //for (std::size_t i = 2; i < localsearch_.size(); i += 1)
    //{
    //    localSearch_[1](this, population_[0]);
    //}

	for (std::size_t i = 0; i < applyLocalSearch_.size(); i += 1)
	{
		applyLocalSearch_[i](this, population_[0]);
	}
}

#pragma endregion

#pragma region Initialization

void MemeticAlgorithm::randomInitialize()
{
    // Assign to Ming rf37535@gmail.com [Done]
    for (auto &i : population_)
    {
        int count = 0;
        for (auto &j : i)
        {
            j = count;
            count += 1;
        }
        std::shuffle(i.begin(), i.end(), RandomRange::RandomGenerator);
    }
}

void MemeticAlgorithm::heuristicInitialize()
{
    // Assign to Ming rf37535@gmail.com
}

#pragma endregion

#pragma region Crossover

void MemeticAlgorithm::OX(Chromosome &first_parent, Chromosome &second_parent)
{
    // TODO: Implement OX to first_parent and second_parent
}

void MemeticAlgorithm::LOX(Chromosome &first_parent, Chromosome &second_parent)
{
    // TODO: Implement LOX to first_parent and second_parent
}

void MemeticAlgorithm::PMX(Chromosome &first_parent, Chromosome &second_parent)
{
    // TODO: Implement PMX to first_parent and second_parent
}

void MemeticAlgorithm::CX(Chromosome &first_parent, Chromosome &second_parent)
{
    // TODO: Implement CX to first_parent and second_parent
}

#pragma endregion

#pragma region Mutation
#pragma endregion

#pragma region LocalSearch

const Chromosome MemeticAlgorithm::II(const Chromosome &chromosome)
{
    // Assign to Shin bazukaoc@gmail.com [Done]
    Chromosome result(chromosome);
    int best = fitness_(result);
    int looptimes = localsearch_looptimes_;
    while (looptimes -= 1)
    {
        for (std::size_t i = 0; i < jobs_ - 1; i += 1)
        {
            for (std::size_t j = i + 1; j < jobs_; j += 1)
            {
                std::swap(result[i], result[j]);
				int score = fitness_(result);
                if (best > score)
                {
                    best = score;
                    return result;
                }
                else
                {
                    std::swap(result[j], result[i]);
                }
            }
        }
    }
    return result;
}

const Chromosome MemeticAlgorithm::SA(const Chromosome &chromosome)
{
    // Assign to Wei [Need fix]
    Chromosome result(chromosome);
	int best = fitness_(result), score;
    int looptimes = localsearch_looptimes_;
    double temperature = 2000;
    int changefirst, changesecond;
    while (looptimes -= 1 && temperature >= 1)
    {
        changefirst = RandomRange::random<int>(0, jobs_ - 1);
        changesecond = RandomRange::random<int>(0, jobs_ - 1);
        std::swap(result[changefirst], result[changesecond]);
        score = fitness_(result);
        if (best > score)
        {
            best = score;
        }
        else
        {
            if(RandomRange::random<double>(0, 1) < exp((best - score) / temperature))
            {
                best = score;
            }
            else
            {
                std::swap(result[changefirst], result[changesecond]);
            }
        }
        temperature *= 0.99;
    }
    return result;
}

const Chromosome MemeticAlgorithm::TS(const Chromosome &chromosome)
{
	// Assign to Ming rf37535@gmail.com [Done]
	const int tabu_length = 7;
	int tabu_current = 0;
	std::vector<int> tabulist(tabu_length, 0);
	Chromosome result(chromosome);
	int looptimes = localsearch_looptimes_;
	while (looptimes -= 1)
	{
		int best = INT_MAX;
		Chromosome current_best;
		for (std::size_t i = 0; i < jobs_ - 1; i += 1)
		{
			for (std::size_t j = i + 1; j < jobs_; j += 1)
			{
				std::swap(result[i], result[j]);
				int score = fitness_(result);
				bool in_tabu = false;
				for (std::size_t k = 0; k < tabu_length; k += 1)
				{
					if (score == tabulist[k])
					{
						in_tabu = true;
					}
				}
				if (!in_tabu)
				{
					if (score < best)
					{
						best = score;
						current_best = result;
					}
				}
				std::swap(result[i], result[j]);
			}
		}
		tabulist[tabu_current] = best;
		tabu_current += 1;
		tabu_current = tabu_current % tabu_length;
		result = current_best;
	}
	return result;
}

#pragma endregion

#pragma region ApplyLocalSearch

const int MemeticAlgorithm::applyLocalSearchByLamarckian(Chromosome &chromosome)
{
	// Implement the Lamarckian function
	// Apply local search and update the chromosome, return updated fitness
	return fitness_(chromosome);
}

const int MemeticAlgorithm::applyLocalSearchByBaldwinian(Chromosome &chromosome)
{
	// Implement the Baldwinian function
	// Only apply local search, return updated fitness
	return fitness_(chromosome);
}

#pragma endregion

#pragma region Helper-Functions

void MemeticAlgorithm::readfile_()
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

const int MemeticAlgorithm::fitness_(const Chromosome& chromosome)
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
