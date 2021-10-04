// SLAE_computational_mathematics.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Matrix.h"
#include "Solutions.h"

int main()
{
	//	Для тестов используется система варианта К
	
	const size_t N = 10;

	mtrx::Square_Matrix<double> system_matrix(N);
	std::vector<double> f(N);
	for (int i = 1; i <= system_matrix.size(); i++)
	{
		for (int j = 1; j <= system_matrix.size(); j++)
		{
			if (i == j)
			{
				system_matrix[i - 1][j - 1] = 1;
			}
			else
			{
				system_matrix[i - 1][j - 1] = 1 / static_cast<double>(i + j);
			}
		}
		f[i - 1] = 1 / static_cast<double>(i);
	}
	std::cout << "A:" << std::endl << system_matrix << std::endl;
	
	std::cout << "f:" << std::endl;

	for (auto& elem : f)
	{
		std::cout.width(12);
		std::cout << elem << std::endl;
	}

	std::cout << std::endl;

	//	Тест метода Гаусса
	{
		auto expanded_system_matrix = mtrx::make_expanded_matrix(system_matrix, f);
		auto solution = slt::Gauss_method(expanded_system_matrix);

		std::cout << "Solution (by Gauss):" << std::endl;

		for (auto& x_i : solution)
		{
			std::cout.width(12);
			std::cout << x_i << std::endl;
		}
	}

	std::cout << std::endl;

	//	Тест метода Зейделя
	
	{
		std::vector<double> start_sol(N, 0);
		auto solution = slt::Seidel_method(system_matrix, f, start_sol, 1e-10);	//	Критерий останова - норма невязки < 10^-10

		std::cout << "Solution (by Seidel):" << std::endl;

		for (auto& x_i : solution)
		{
			std::cout.width(12);
			std::cout << x_i << std::endl;
		}

	}
	
	//std::cout << system_matrix.det() << std::endl;
	/*
	mtrx::Square_Matrix<double> m = { {0.333333,  1}, {0.25, 0.2} };
	std::cout << m.det() << std::endl;
	*/
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
