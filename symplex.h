#pragma once
#include "vector.h"
#include "matrix.h"
#include <iomanip>
#include <string>
#include "numeric.h"
using namespace std;
////////////////////
/// Lab. work #5 ///
////////////////////
namespace sm
{
	constexpr int EQUAL = 0;
	constexpr int LESS_EQUAL = -1;
	constexpr int MORE_EQUAL = 1;


	class simplex
	{
	private:
		std::vector<int> _inequations;
		std::vector<int> virtualArgsIds;
		std::vector<int> naturalArgsIds;
		std::vector<int> basisArgsIds;
		mat_mn simplex_t;
		mat_mn bounds_m;
		vec_n bounds_v;
		vec_n prices_v;
		int mode = SIMPLEX_MAX;
		bool isPlanOptimal()const;
		int getMainCol()const;
		int getMainRow(const int simplex_col)const;
		bool buildVirtualBasisCol(const int ineq_id, const int _ineq, int& col_index, int& col_index_aditional);
		///		|u 0 0|	
		/// A = |0 v 0|
		///		|0 0 w|
		/// ������ �����������
		///		|a|	
		/// b = |d|
		///		|f|
		/// � -������������ ������� ������� 
		/// f = (x,c)->extr
		///	|u 0 0|   |x| <= |b|
		/// |0 v 0| * |x| >= |f|
		///	|0 0 w|   |x| =  |d|
		/// 
		///  �� ������� �� A,b,c ����������
		/// </summary>
		/// <param name="A"> Ax <= b   -> (A|I)(x|w) = b </param>
		/// <param name="c"> (c,x) ->((-c|0),(x|w)) </param>
		/// <param name="ineq"> ���� ����������� =, >=, <= </param>
		/// <param name="b"></param>
		///( A|I)  b
		///(-c|0)  F(x,c)
		void buildSimplexTable();

		bool excludeVirtualArgs();

		bool validateSolution()const;

	public:
		int naturalArgsN()const;

		inline const mat_mn& boundsMatrix()const;

		inline const vec_n& boundsCoeffs()const;

		inline const vec_n& pricesCoeffs()const;

		inline const std::vector<int>& inequations()const;

		inline const std::vector<int>& basisArgumentsIds()const;

		inline const mat_mn& table()const;

		inline bool                    isTargetFunctionModified()const;
		vec_n                          currentSimplexSolution(const bool only_natural_args = false)const;

		vec_n                          solve(const int mode = SIMPLEX_MAX);

		simplex(const mat_mn& a, const vec_n& c, const std::vector<int>& _ineq, const vec_n& b);

		simplex(const mat_mn& a, const vec_n& c, const vec_n& b);

		friend std::ostream& operator<<(std::ostream& stream, const simplex& s);
	};
	std::ostream& operator<<(std::ostream& stream, const simplex& s)
	{

		const  mat_mn& A = s.table();

		const bool targeFuncMod = s.isTargetFunctionModified();

		const std::vector<int>& basisArgsIds = s.basisArgumentsIds();

		if (A.size() == 0)
		{
			return stream;
		}

		const char separator = ' ';
		const int colom_title_w = 6;
		const int colom_w = 12;

		stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << "";

		int i = 0;

		for (; i < A[0].size() - 1; i++)
		{
			stream << std::left << std::setw(colom_w) << std::setfill(separator) << "| x " + std::to_string(i + 1);
		}

		stream << std::left << std::setw(colom_w) << std::setfill(separator) << "| b";

		stream << "\n";

		int n_row = -1;

		for (auto const& row : A)
		{
			n_row++;

			if (targeFuncMod)
			{
				if (n_row == A.size() - 2)
				{
					stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << " d0 ";
				}
				else if (n_row == A.size() - 1)
				{
					stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << " d1 ";
				}
				else
				{
					stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << " x " + std::to_string(basisArgsIds[n_row] + 1);
				}
			}
			else
			{
				if (n_row == A.size() - 1)
				{
					stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << " d ";
				}
				else
				{
					stream << std::left << std::setw(colom_title_w) << std::setfill(separator) << " x " + std::to_string(basisArgsIds[n_row] + 1);
				}
			}

			for (int col = 0; col < row.size(); col++)
			{
				if (row[col] >= 0)
				{
					stream << std::left << std::setw(colom_w) << std::setfill(separator) << "| " + str_rational(row[col]);
					continue;
				}
				stream << std::left << std::setw(colom_w) << std::setfill(separator) << "|" + str_rational(row[col]);
			}

			stream << "\n";
		}
		stream << "\n";

		return stream;
	}
	bool simplex::isPlanOptimal()const
	{
		/// <summary>
		/// ��������� �������� ��������� ������ �������-���������
		/// �� ���������������. ���� ��� ������������, �� ���� ���������.
		/// </summary>

		const vec_n& row = simplex_t[simplex_t.size() - 1];

		bool opt = true;

		for (int i = 0; i < row.size() - 1; i++)
		{
			if (row[i] < 0)
			{
				opt = false;
				break;
			}
		}

		/// <summary>
		/// ���� �� �������������� ������� �������, �� ����� ������ �������������
		/// ��������� �������� �� ���������������� ������������ ������ ��������-��������� 
		/// </summary>

		if (isTargetFunctionModified())
		{
			if (!opt) return opt;

			const vec_n& row_ = simplex_t[simplex_t.size() - 2];

			for (const auto& id : naturalArgsIds)
			{
				if (row_[id] < 0)
				{
					opt = false;
					break;
				}
			}
		}

		return opt;
	}
	int simplex::getMainCol()const
	{
		const vec_n& row = simplex_t[simplex_t.size() - 1];

		double delta = 0;

		int    index = -1;

		for (int i = 0; i < row.size() - 1; i++)
		{
			if (row[i] >= delta)continue;
			delta = row[i];
			index = i;
		}

		if (isTargetFunctionModified() && index == -1)
		{
			const vec_n& row_add = simplex_t[simplex_t.size() - 2];

			for (const auto& id : naturalArgsIds)
			{
				if (row_add[id] >= delta)continue;
				delta = row_add[id];
				index = id;
			}
		}
		return index;
	}
	int simplex::getMainRow(const int simplex_col)const
	{
		double delta = 1e12;

		int index = -1;

		double a_ik;

		int b_index = simplex_t[0].size() - 1;

		int cntr = 0;

		int rows_n = isTargetFunctionModified() ? simplex_t.size() - 2 : simplex_t.size() - 1;

		for (int i = 0; i < rows_n; i++)
		{
			a_ik = simplex_t[i][simplex_col];

			if (a_ik < 0)
			{
				cntr++;
				continue;
			}
			if (simplex_t[i][b_index] / a_ik > delta)continue;
			delta = simplex_t[i][b_index] / a_ik;
			index = i;
		}

		return index;
	}
	bool simplex::buildVirtualBasisCol(const int ineq_id, const int _ineq, int& col_index, int& col_index_aditional)
	{
		if (_ineq == EQUAL)
		{
			for (int row = 0; row < simplex_t.size(); row++)
			{
				if (row == ineq_id)
				{
					simplex_t[row].push_back(1.0);
					continue;
				}
				simplex_t[row].push_back(0.0);
			}

			col_index = simplex_t[0].size() - 1;

			col_index_aditional = simplex_t[0].size() - 1;

			return true;
		}

		if (_ineq == MORE_EQUAL)
		{
			for (int row = 0; row < simplex_t.size(); row++)
			{
				if (row == ineq_id)
				{
					simplex_t[row].push_back(-1.0);

					simplex_t[row].push_back(1.0);

					continue;
				}

				simplex_t[row].push_back(0.0);

				simplex_t[row].push_back(0.0);
			}

			col_index_aditional = simplex_t[0].size() - 1;

			col_index = simplex_t[0].size() - 2;

			return false;
		}

		for (int row = 0; row < simplex_t.size(); row++)
		{
			if (row == ineq_id)
			{
				simplex_t[row].push_back(1.0);
				continue;
			}
			simplex_t[row].push_back(0.0);
		}

		col_index_aditional = -1;

		col_index = simplex_t[0].size() - 1;

		return true;
	}
	void simplex::buildSimplexTable()
	{
		simplex_t = bounds_m;
		///
		/// ���� ����� ������� b ���� ������������� ��������, �� ��������������� ������
		/// ������� ����������� �������� �� ���� ���� � ������ ���� ���������
		///
		for (int row = 0; row < simplex_t.size(); row++)
		{
			if (bounds_v[row] >= 0)continue;

			_inequations[row] *= -1;

			bounds_v[row] *= -1;

			simplex_t[row] = simplex_t[row] * (-1.0);
		}


		for (int i = 0; i < prices_v.size(); i++)
		{
			naturalArgsIds.push_back(i);
		}
		/// <summary>
		/// ���������� ������������� ������
		/// </summary>
		int basis_arg_id;
		int basis_arg_id_add;
		for (int ineq_id = 0; ineq_id < _inequations.size(); ineq_id++)
		{
			buildVirtualBasisCol(ineq_id, _inequations[ineq_id], basis_arg_id, basis_arg_id_add);

			naturalArgsIds.push_back(basis_arg_id);

			if (basis_arg_id_add != -1)
			{
				basisArgsIds.push_back(basis_arg_id_add);
				virtualArgsIds.push_back(basis_arg_id_add);
				continue;
			}

			basisArgsIds.push_back(basis_arg_id);
		}

		/// <summary>
		/// ������� ������� �����������
		/// </summary>

		for (int row = 0; row < simplex_t.size(); row++)
		{
			simplex_t[row].push_back(bounds_v[row]);
		}

		/// <summary>
		/// ���������� �������� ���������
		/// </summary>

		vec_n s_deltas(simplex_t[0].size());

		if (mode == SIMPLEX_MAX)
		{
			for (int j = 0; j < s_deltas.size(); j++) s_deltas[j] = j < prices_v.size() ? -prices_v[j] : 0.0;
		}
		else
		{
			for (int j = 0; j < s_deltas.size(); j++) s_deltas[j] = j < prices_v.size() ? prices_v[j] : 0.0;
		}

		simplex_t.push_back(s_deltas);

		/// <summary>
		/// ���� ������� �������� �� ���� ��������������
		/// </summary>

		if (!isTargetFunctionModified()) return;

		/// <summary>
		/// ���� �� �� ����...
		/// </summary>
		vec_n s_deltas_add(simplex_t[0].size());

		for (int j = 0; j < virtualArgsIds.size(); j++) s_deltas_add[virtualArgsIds[j]] = 1.0;

		simplex_t.push_back(s_deltas_add);
	}
	bool simplex::excludeVirtualArgs()
	{
		if (!isTargetFunctionModified()) return false;

		int last_row_id = simplex_t.size() - 1;

		for (int i = 0; i < virtualArgsIds.size(); i++)
		{
			for (int row = 0; row < simplex_t.size(); row++)
			{
				if (simplex_t[row][virtualArgsIds[i]] != 0)
				{
					double arg = simplex_t[last_row_id][virtualArgsIds[i]] / simplex_t[row][virtualArgsIds[i]];

					simplex_t[last_row_id] = simplex_t[last_row_id] - arg * simplex_t[row];

					break;
				}
			}
		}

		return true;
	}
	bool simplex::validateSolution()const
	{
		double val = 0;

		int n_rows = isTargetFunctionModified() ? simplex_t.size() - 2 : simplex_t.size() - 1;

		int n_cols = simplex_t[0].size() - 1;

		for (int i = 0; i < basisArgsIds.size(); i++)
		{
			if (basisArgsIds[i] < naturalArgsN())
			{
				val += simplex_t[i][n_cols] * prices_v[basisArgsIds[i]];
			}
		}
		if (mode == SIMPLEX_MAX)
		{
			if (abs(val - simplex_t[n_rows][n_cols]) < N_DIM_ACCURACY)
			{
				if (isTargetFunctionModified())
				{
					return true & (abs(simplex_t[simplex_t.size() - 1][simplex_t[0].size() - 1]) < N_DIM_ACCURACY);
				}

				return true;
			}
		}

		if (abs(val + simplex_t[n_rows][n_cols]) < N_DIM_ACCURACY)
		{
			if (isTargetFunctionModified())
			{
				return true & (abs(simplex_t[simplex_t.size() - 1][simplex_t[0].size() - 1]) < N_DIM_ACCURACY);
			}

			return true;
		}
		return false;
	}
	int simplex::naturalArgsN()const
	{
		return prices_v.size();
	}
	inline const mat_mn& simplex::boundsMatrix()const
	{
		return bounds_m;
	}
	inline const vec_n& simplex::boundsCoeffs()const
	{
		return bounds_v;
	}
	inline const vec_n& simplex::pricesCoeffs()const
	{
		return prices_v;
	}
	inline const std::vector<int>& simplex::inequations()const
	{
		return _inequations;
	};
	inline const std::vector<int>& simplex::basisArgumentsIds()const
	{
		return basisArgsIds;
	};
	inline const mat_mn& simplex::table() const { return simplex_t; };
	inline bool simplex::isTargetFunctionModified()const
	{
		return virtualArgsIds.size() != 0;
	}
	vec_n simplex::currentSimplexSolution(const bool only_natural_args)const
	{
		vec_n solution(only_natural_args ? naturalArgsN() : simplex_t[0].size() - 1);

		for (int i = 0; i < basisArgsIds.size(); i++)
		{
			if (basisArgsIds[i] >= solution.size()) continue;

			solution[basisArgsIds[i]] = simplex_t[i][simplex_t[0].size() - 1];
		}
		return solution;
	}
	vec_n simplex::solve(const int mode)
	{
		this->mode = mode;

		std::cout << "Simplex problem type: " << ((mode == SIMPLEX_MAX) ? "max\n" : "min\n");

		buildSimplexTable();

		vec_n solution;

		double a_ik;

		int main_row;

		int main_col;

		std::cout << "Start simplex table:" << "\n";

		std::cout << *this;

		if (excludeVirtualArgs())
		{
			// ������ ����, ���� ������ ������ �������� ���� ���������(���� �������) ����������
			std::cout << "Simplex table after args exclude:" << "\n";

			std::cout << *this;
		}

		while (!isPlanOptimal())
		{
			main_col = getMainCol();

			if (main_col == -1) break;

			main_row = getMainRow(main_col);

			if (main_row == -1)
			{
				/// ������������� ���������� ������� ������ ���������������� � ���, ��� ������� ������ ������������
				std::cout << "Unable to get main row. Simplex is probably boundless...\n";
				solution.clear();
				return solution;
			}

			basisArgsIds[main_row] = main_col;

			a_ik = simplex_t[main_row][main_col];

			simplex_t[main_row] = simplex_t[main_row] * (1.0 / a_ik);

			for (int i = 0; i < simplex_t.size(); i++)
			{
				if (i == main_row)continue;

				simplex_t[i] = simplex_t[i] - simplex_t[i][main_col] * simplex_t[main_row];
			}
			solution = currentSimplexSolution();

#if _DEBUG
			std::cout << "a_main { " << main_row + 1 << ", " << main_col + 1 << " } = " << str_rational(a_ik) << "\n";
			std::cout << *this;
			std::cout << "current_solution" << str_rational(solution) << "\n";
			std::cout << "\n";
#endif
		}
		if (validateSolution())
		{
			solution = currentSimplexSolution(true);
			/// ������������ ������
			std::cout << "solution : " << str_rational(solution) << "\n";
			return solution;
		}
		std::cout << "Simplex is unresolvable\n";
		/// �������� ������� ������� �� ����� �� �������� �� ���������� �����
		solution.clear();
		return solution;
	}
	simplex::simplex(const mat_mn& a, const vec_n& c, const std::vector<int>& _ineq, const vec_n& b)
	{
		if (b.size() != _ineq.size())
		{
			throw std::runtime_error("Error simplex creation :: b.size() != inequalities.size()");
		}
		if (a.size() != _ineq.size())
		{
			throw std::runtime_error("Error simplex creation :: A.rows_number() != inequalities.size()");
		}

		if (a[0].size() != c.size())
		{
			throw std::runtime_error("Error simplex creation :: A.cols_number() != price_coeffs.size()");
		}

		bounds_v = b;

		bounds_m = a;

		prices_v = c;

		_inequations = std::vector<int>(_ineq);
	}
	simplex::simplex(const mat_mn& a, const vec_n& c, const vec_n& b)
	{
		if (a.size() != b.size())
		{
			throw std::runtime_error("Error simplex creation :: A.rows_number() != bouns_coeffs.size()");
		}

		if (a[0].size() != c.size())
		{
			throw std::runtime_error("Error simplex creation :: A.cols_number() != price_coeffs.size()");
		}

		std::vector<int> _ineq;

		for (int i = 0; i < b.size(); i++)
		{
			_ineq.push_back(LESS_EQUAL);
		}

		bounds_v = b;

		bounds_m = a;

		prices_v = c;

		_inequations = std::vector<int>(_ineq);
	}

}