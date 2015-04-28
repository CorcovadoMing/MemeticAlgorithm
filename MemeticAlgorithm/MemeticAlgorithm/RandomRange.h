#pragma once
#include <random>

namespace RandomRange
{

	std::random_device rd;
	std::mt19937 rng(rd());

	template <class T>
	struct TypeIsInt
	{
		static const bool value = false;
	};

	template <>
	struct TypeIsInt<int>
	{
		static const bool value = true;
	};

	template <class T>
	T random(T min, T max)
	{
		if (TypeIsInt<T>::value)
		{
			std::uniform_int_distribution<int> uniform_int(min, max);
			return uniform_int(rng);
		}
		else
		{
			std::uniform_real_distribution<double> uniform_real(min, max);
			return uniform_real(rng);
		}
	}
}