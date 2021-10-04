#pragma once
#include <vector>
#include <exception>

#include "Matrix.h"

namespace slt
{

	template <typename T = double>
	std::vector<T> Gauss_method(const mtrx::Matrix<T>&);
	template <typename T = double>
	std::vector<T> Seidel_method(const mtrx::Square_Matrix<T>&, const std::vector<T>&, const std::vector<T>, const double);

	//	Метод Гаусса

	template <typename T = double>
	std::vector<T> Gauss_method(const mtrx::Matrix<T>& const_input_matrix) throw (std::exception&)
	{
		auto matrix = const_input_matrix;
		if (matrix.column_size() + 1 != matrix.row_size())
			throw std::exception("Error: invalid format of matrix");

		std::vector<std::pair<size_t, T>> solutions(matrix.row_size() - 1);
		for (int i = 0; i < solutions.size(); i++)
			solutions[i].first = i;

		//	Прямой ход

		for (int i = 0; i < matrix.column_size() - 1; i++)
		{
			mtrx::transform_matrix_with_main_element(matrix, solutions, i);
			for (int j = i + 1; j < matrix.column_size(); j++)
			{
				if (std::abs(matrix[i][i]) < 10e-6)
				{
					throw std::exception("Error: matrix doesn't allow to find the only solution");
				}
				double q = static_cast<double>(matrix[j][i]) / matrix[i][i];
				for (int k = i; k < matrix.row_size(); k++)
				{
					matrix[j][k] -= matrix[i][k] * q;
				}
			}
			
		}

		//	Обратный ход

		for (int i = matrix.row_size() - 2; i >= 0; i--)
		{
			solutions[i].second = matrix[i][matrix.row_size() - 1];
			for (int j = i + 1; j < solutions.size(); j++)
			{
				solutions[i].second -= solutions[j].second * matrix[i][j];
			}
			if (std::abs(matrix[i][i]) < 10e-6)
			{
				throw std::exception("Error: matrix doesn't allow to find the only solution");
			}
			solutions[i].second /= matrix[i][i];
		}

		//	Восстановление порядка переменных (т.к. используется метод с поиском главного элемента)

		std::sort(solutions.begin(), solutions.end(), [](auto& pair1, auto& pair2) -> bool { return pair1.first < pair2.first; });

		std::vector<T> solutions_only(solutions.size());
		for (int i = 0; i < solutions_only.size(); i++)
		{
			solutions_only[i] = solutions[i].second;
		}

		//	Расчёт невязки
		std::vector<T> r(solutions.size());
		std::fill(r.begin(), r.end(), 0);
		for (int i = 0; i < matrix.column_size(); i++)
		{
			for (int j = 0; j < matrix.row_size() - 1; j++)
			{
				r[i] += const_input_matrix.get_xy(j, i) * solutions_only[j];
			}
			r[i] -= const_input_matrix.get_xy(matrix.row_size() - 1, i);
		}
		std::cout << "r_norm (in Gauss method):" << std::endl << mtrx::norm(r) << std::endl << std::endl;

		return solutions_only;
		
	}

	//	Метод Зейделя
	template <typename T = double>
	std::vector<T> Seidel_method(const mtrx::Square_Matrix<T>& sqr_mat, const std::vector<T>& f, const std::vector<T> start_sol, const double accuracy) throw (std::exception&)
	{

		//	Разложение матрицы на L, D и U матрицы

		mtrx::Square_Matrix<T> L(sqr_mat.size()), D(sqr_mat.size()), U(sqr_mat.size());
		for (int i = 0; i < L.size(); i++)
		{
			for (int j = 0; j < L.size(); j++)
			{
				if (i > j)
				{
					L[i][j] = sqr_mat.get_xy(j, i);
				}
				if (i == j)
				{
					D[i][j] = sqr_mat.get_xy(j, i);
				}
				if (i < j)
				{
					U[i][j] = sqr_mat.get_xy(j, i);
				}
			}
		}

		//	Расчёт матриц B и F для итерационного метода

		auto B = -(L + D).inverse() * U;
		std::vector<T> F = (L + D).inverse() * f;

		//	Инициализация начального приближения решения

		std::vector<T> solution = start_sol;

		//	Начальное значение невязки и её нормы-3

		auto r = sqr_mat * solution - f;
		auto r_norm = mtrx::norm(r);

		std::cout << "r_norm (in Seidel):" << std::endl;

		std::cout.width(18);
		std::cout << "Iteration number";
		std::cout.width(18);
		std::cout << "r_norm" << std::endl;

		//	Итерационный цикл
		int number = 0;
		while (r_norm > accuracy)
		{
			solution = B * solution + F;
			r = sqr_mat * solution - f;
			r_norm = mtrx::norm(r);

			//	Вывод значения невязки

			std::cout.width(18);
			std::cout << number;
			std::cout.width(18);
			std::cout << r_norm << std::endl;
			
			number++;
		}
		std::cout << std::endl;
		return solution;
	}

}	//	namespace slt