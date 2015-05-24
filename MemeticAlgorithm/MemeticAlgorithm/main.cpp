#include "MemeticAlgorithm.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <ctime>

void experiment(const int a, const int b, const int c, const int d, const int e, const std::string& filename)
{
	clock_t start_time = clock();
	std::vector<int> result;
	for (std::size_t loop = 0; loop < 100; loop += 1)
	{
		MemeticAlgorithm ma(100, 0.2, 0.1, filename);
		result.push_back(ma.run(a, b, c, d, e)); //2, 4, 3, 2, 3
		// default: 1, 3, 1, 0, 0
	}
	double sum = std::accumulate(result.begin(), result.end(), 0);
	double mean = sum / result.size();
	double squrtSum = std::inner_product(result.begin(), result.end(), result.begin(), 0.0);
	double var = std::sqrt(squrtSum / result.size() - mean * mean);

	std::cout << std::endl;
	std::cout << "Experiment: " << filename << " " << a << " " << b << " " << c << " " << d << " " << e << std::endl;
	std::cout << "Max: " << *std::max_element(result.begin(), result.end()) << std::endl;
	std::cout << "Min: " << *std::min_element(result.begin(), result.end()) << std::endl;
	std::cout << "Avg: " << mean << std::endl;
	std::cout << "Var: " << var << std::endl;
	std::cout << "Time: " << (clock() - start_time) / (double)CLOCKS_PER_SEC << std::endl;
}

int main() 
{
	const std::string dataset[] = { "tai20_5_1.txt", "tai20_10_1.txt", "tai20_20_1.txt", "tai50_5_1.txt", "tai50_10_1.txt", "tai50_20_1.txt", "tai100_5_1.txt", "tai100_10_1.txt", "tai100_20_1.txt" };
	for (const std::string i : dataset)
	{
		for (std::size_t a = 0; a < 2; a += 1)
		{
			experiment(a, 3, 1, 0, 0, i);
		}
		for (std::size_t b = 0; b < 4; b += 1)
		{
			experiment(1, b, 1, 0, 0, i);
		}
		for (std::size_t c = 0; c < 3; c += 1)
		{
			experiment(1, 3, c, 0, 0, i);
		}
		for (std::size_t d = 0; d < 2; d += 1)
		{
			experiment(1, 3, 1, d, 0, i);
		}
		for (std::size_t e = 0; e < 3; e += 1)
		{
			experiment(1, 3, 1, 0, e, i);
		}
	}
	return 0;
}