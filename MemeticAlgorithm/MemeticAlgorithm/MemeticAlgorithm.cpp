#include "MemeticAlgorithm.h"
#include "RandomRange.h"
#include <iostream>
#include <fstream>
#include <algorithm>

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

	mutation_.push_back(std::mem_fn(&MemeticAlgorithm::insertion));
	mutation_.push_back(std::mem_fn(&MemeticAlgorithm::swap));
	mutation_.push_back(std::mem_fn(&MemeticAlgorithm::inverse));

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
    std::cout << "=== Testing Initialize ===" << std::endl;
    for (std::size_t i = 0; i < initialize_.size(); i += 1)
    {
        initialize_[i](this);
    }
	std::cout << "=== Testing Crossover ===" << std::endl;
    for (std::size_t i = 0; i < crossover_.size(); i += 1)
    {
		crossover_[i](this, population_[0], population_[1]);
    }
	std::cout << "=== Testing Mutation ===" << std::endl;
    for (std::size_t i = 0; i < mutation_.size(); i += 1)
    {
		mutation_[i](this, population_[0]);
    }
	std::cout << "=== Testing Localsearch ===" << std::endl;
    for (std::size_t i = 0; i < localSearch_.size(); i += 1)
    {
        localSearch_[i](this, population_[0]);
    }
	std::cout << "=== Testing ApplyLocalsearch ===" << std::endl;
	for (std::size_t i = 0; i < applyLocalSearch_.size(); i += 1)
	{
		applyLocalSearch_[i](this, population_[0], 0);
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
    // Assign to Ming rf37535@gmail.com [Done]
	const std::size_t heuristic_solution_size = population_size_ / 10;
	randomInitialize();
	for (std::size_t i = 0; i < heuristic_solution_size; i += 1)
	{
		localSearch_[0](this, population_[i]); // II
	}
}

#pragma endregion

#pragma region Mating-Selection

void MemeticAlgorithm::tournament()
{

}

#pragma endregion

#pragma region Crossover

void MemeticAlgorithm::OX(const Chromosome &first_parent, const Chromosome &second_parent)
{
    std::cout << "OX" << std::endl;
    std::size_t chromosome_size = first_parent.size();
	std::size_t inherit_index = RandomRange::random<int>(0, chromosome_size - 2);
//	cout << "inherit_index = " << inherit_index << endl;
	std::size_t inherit_length = RandomRange::random<int>(inherit_index + 1, chromosome_size - 1) - inherit_index;
//  cout << "inherit_length = " << inherit_length << endl;

    //first_child inherit first_parent, and other is second_parent, use first_temp
    Chromosome first_temp(second_parent), first_child(first_parent);
    //second_child inherit second_parent, and other is first_parent, use second_temp
    Chromosome second_temp(first_parent), second_child(second_parent);

    //mark which had been inherited
    for(std::size_t i = 0; i < chromosome_size; i += 1)
    {
        for(std::size_t j = 0; j < inherit_length; j += 1)
        {
            if(first_temp[i] == first_parent[inherit_index+j])
            {
                first_temp[i] = -1;
            }
            if(second_temp[i] == second_parent[inherit_index+j])
            {
                second_temp[i] = -1;
            }
        }
    }

    //add temp to child
    for(std::size_t i = inherit_index + inherit_length, first_addindex = 0, second_addindex = 0; i < chromosome_size + inherit_index; i += 1)
    {
        while(first_temp[first_addindex] == -1 && first_addindex < chromosome_size)
        {
            first_addindex += 1;
        }
        while(second_temp[second_addindex] == -1 && second_addindex < chromosome_size)
        {
            second_addindex += 1;
        }

        first_child[i % chromosome_size] = first_temp[first_addindex];
        second_child[i % chromosome_size] = second_temp[second_addindex];

        first_addindex += 1;
        second_addindex += 1;
    }
	offspring_.push_back(first_child);
	offspring_.push_back(second_child);
}

void MemeticAlgorithm::LOX(const Chromosome &first_parent, const Chromosome &second_parent)
{
    std::size_t chromosome_size = first_parent.size();
    std::size_t inherit_index = RandomRange::random<int>(0, chromosome_size - 2);
    std::size_t inherit_length = RandomRange::random<int>(inherit_index + 1, chromosome_size - 1) - inherit_index; 
    Chromosome first_temp(second_parent), first_child(first_parent);
    Chromosome second_temp(first_parent), second_child(second_parent);
    for (std::size_t i = 0; i < chromosome_size; i += 1)
    {
        for (std::size_t j = 0; j < inherit_length; j += 1)
        {
            if (first_temp[i] == first_parent[inherit_index+j])
            {
                first_temp[i] = -1;
            }
            if (second_temp[i] == second_parent[inherit_index+j])
            {
                second_temp[i] = -1;
            }
        }
    }

    for (std::size_t i = 0, k1 = 0, k2 = 0; i < chromosome_size; i += 1)
    {
        if (i == inherit_index)
        {
            i += (inherit_length - 1);
        }
        else
        {
            while (first_temp[k1]  == -1 && k1 < chromosome_size) {
                k1 += 1;
            }
            while (second_temp[k2] == -1 && k2 < chromosome_size) {
                k2 += 1;
            }
            first_child[i]  = first_temp[k1];  if (k1 < chromosome_size) { k1 += 1; }
            second_child[i] = second_temp[k2]; if (k2 < chromosome_size) { k2 += 1; }
        }
    }
    offspring_.push_back(first_child);
    offspring_.push_back(second_child);
}

void MemeticAlgorithm::PMX(const Chromosome &first_parent, const Chromosome &second_parent)
{
    std::size_t chromosome_size = first_parent.size();
    std::size_t inherit_index  = RandomRange::random<int>(0, chromosome_size - 2);
    std::size_t inherit_length = RandomRange::random<int>(inherit_index + 1, chromosome_size - 1) - inherit_index;

    Chromosome first_child(first_parent);
    Chromosome second_child(second_parent);
    for (std::size_t i = inherit_index; i < inherit_length + inherit_index; i += 1)
    {
        if (first_parent[i] != second_parent[i])
        {
            for (std::size_t j = 0; j < chromosome_size; j += 1)
            {
                if (first_child[j] == second_parent[i])
                {
                    std::swap(first_child[i], first_child[j]);
                }
                if (second_child[j] == first_parent[i])
                {
                    std::swap(second_child[i], second_child[j]);
                }
            }
        }
    }
    offspring_.push_back(first_child);
	offspring_.push_back(second_child);
}

void MemeticAlgorithm::CX(const Chromosome &first_parent, const Chromosome &second_parent)
{
	Chromosome first_child(first_parent);
	Chromosome second_child(second_parent);
	std::size_t check_point = RandomRange::random<int>(0, first_parent.size() - 1);
	std::vector<std::size_t> exchange_element;
	const int start_element = first_parent[check_point];
	while (second_parent[check_point] != start_element)
	{
		exchange_element.push_back(check_point);
		for (std::size_t i = 0; i < first_parent.size(); i += 1)
		{
			if (first_parent[i] == second_parent[check_point])
			{
				check_point = i;
				break;
			}
		}
	}
	for (std::size_t i = 0; i < exchange_element.size(); i += 1)
	{
		std::swap(first_child[exchange_element[i]], second_child[exchange_element[i]]);
	}
	offspring_.push_back(first_child);
	offspring_.push_back(second_child);
}

#pragma endregion

#pragma region Mutation

void MemeticAlgorithm::insertion(Chromosome &chromosome)
{
	const std::size_t point = RandomRange::random<int>(1, chromosome.size() - 1);
	const std::size_t insert = RandomRange::random<int>(0, point - 1);

	std::size_t i = point;
	while (i != insert)
	{
		std::swap(chromosome[i], chromosome[i - 1]);
		i -= 1;
	}
}

void MemeticAlgorithm::swap(Chromosome &chromosome)
{
    std::size_t a = RandomRange::random<int>(0, chromosome.size() - 1);
    std::size_t b = RandomRange::random<int>(0, chromosome.size() - 1);
    while (a == b) 
	{
        b = RandomRange::random<int>(0, chromosome.size() - 1);
    }
    std::swap(chromosome[a], chromosome[b]);
}

void MemeticAlgorithm::inverse(Chromosome &chromosome)
{
    std::size_t start = RandomRange::random<int>(0, chromosome.size() - 2);
    std::size_t range = RandomRange::random<int>(2, chromosome.size() - start);
    for(int i = 0; i < range / 2; i += 1)
    {
        std::swap(chromosome[start + i], chromosome[start + range - 1 - i]);
    }
}

#pragma endregion

#pragma region Environment-Selection

void MemeticAlgorithm::tophalf()
{

}

#pragma endregion

#pragma region LocalSearch

const Chromosome MemeticAlgorithm::II(const Chromosome &chromosome)
{
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

const int MemeticAlgorithm::applyLocalSearchByLamarckian(Chromosome &chromosome, const int type)
{
	chromosome = localSearch_[type](this, chromosome);
	return fitness_(chromosome);
}

const int MemeticAlgorithm::applyLocalSearchByBaldwinian(Chromosome &chromosome, const int type)
{
	return fitness_(localSearch_[type](this, chromosome));
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
