#include "vector.h"
#include "matrix.h"
#include "numeric.h"
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
#include "symplex.h"
////////////////////
/// Lab. work #5 ///
////////////////////

static double testFunction1d(const double& x)
{
	return  (x - 2) * (x - 5);
}

static double testFunction2d(const vec_n& args)
{
	return (args[0] - 5) * args[0] + (args[1] - 3) * args[1];
}
static double testFunctionNd(const vec_n& args)
{
	double  val = 0;
	for (auto const& x : args)
	{
		val += testFunction1d(x);
	}
	return val;
}

static void simplexMethodTest(mat_mn price,  vector<int> buyer, vector<int> seller, int m, int n)
{
	int MM = m + n;
	int NN = m * n;
	vec_n coef_func; // коэфициенты целевой функции
	//mat_mn coef = zeros(MM+NN,NN); // коэфициенты уравнений ограничений
	mat_mn coef = zeros(MM, NN); // коэфициенты уравнений ограничений
	vec_n func_b;// коэфы ограничений
	vector<int> equals; // кэфы неравенств
	for (int i = 0; i < m; i++) {

		for ( int j = 0; j < n; j++ )
		{
			coef[i][j + i * n ] = 1;
		
		}
	}
	for (int i = 0; i < n; i++) {

		for (int j = 0; j < m; j++)
		{
			coef[i+m][ i + n*j] = 1;
		}
	}

	/*for (int i = 0; i < n*m; i++) {

			coef[i + m+n][i] = 1;
	}*/
	
	//cout << coef;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
		{
			coef_func.push_back(price[i][j]);
		}
	}
	for (int i = 0; i < m; i++) {
		func_b.push_back(buyer[i]);
		equals.push_back(sm::EQUAL);
	}
	for (int i = 0; i < n; i++)
	{
		func_b.push_back(seller[i]);
		equals.push_back(sm::EQUAL);
	}
	/*for (int i = 0; i < n*m; i++)
	{
		func_b.push_back(0);
		equals.push_back(sm::MORE_EQUAL);
	}*/

	sm::simplex sym_1(coef, coef_func, equals, func_b);
	sym_1.solve(SIMPLEX_MAX);
	std::cout << "\n";


}

static void simplexMethodTest()
{
	std::cout << "\n/////////////////////////////" << std::endl;
	std::cout << "//////// SimplexTest ////////" << std::endl;
	std::cout << "/////////////////////////////\n" << std::endl;
	std::cout << " f(x,c) =  2x1 + 3x2+x3;\n arg_max = {4, 8}, f(arg_max) = 32\n";
	std::cout << " |-2x1 + 6x2 +x3<= 40\n";
	std::cout << " | 3x1 + 2x2 +x3<= 28\n";
	std::cout << " | 2x1 -  x2 +x3<= 14\n\n";

	sm::simplex sym_0({ {2, 12}, { 4,6}, { 3,0}, { 0, 18} }, { 12, 10}, { sm::MORE_EQUAL, sm::MORE_EQUAL, sm::MORE_EQUAL, sm::MORE_EQUAL }, { 20,32,14,42 });
	sym_0.solve(SIMPLEX_MAX);
	std::cout << "\n";

}
/*static void simplexMethodTest()
{
	std::cout << "\n/////////////////////////////" << std::endl;
	std::cout << "//////// SimplexTest ////////" << std::endl;
	std::cout << "/////////////////////////////\n" << std::endl;
	std::cout << " f(x,c) =  2x1 + 3x2+x3;\n arg_max = {4, 8}, f(arg_max) = 32\n";
	std::cout << " |-2x1 + 6x2 +x3<= 40\n";
	std::cout << " | 3x1 + 2x2 +x3<= 28\n";
	std::cout << " | 2x1 -  x2 +x3<= 14\n\n";

	sm::simplex sym_0({ {2, 3,6}, { 4, 2,4}, { 4,6,8}, { 2, 3 ,1 } }, { 5, 2, 1 }, { sm::LESS_EQUAL, sm::LESS_EQUAL, sm::LESS_EQUAL, sm::LESS_EQUAL }, { 40, 28, 14, 20 });
	sym_0.solve(SIMPLEX_MAX);
	std::cout << "\n";

}*/

 int main()
{
	 setlocale(LC_ALL, "Russian");

	 int m, n;
	 vector<int> buyer;
	 vector<int> seller;
	 cout << "Введите количество поставщиков и покупателей>";
	 cin >> m >> n;
	 int sum_buyer = 0, sum_seller = 0;
	 
		  sum_buyer = 0;
		  sum_seller = 0;
		 cout << "Введите запасы поставщиков>";
		 for (int i = 0; i < m; i++)
		 {
			 int temp;
			 cin >> temp;
			 buyer.push_back(temp);
			 sum_buyer += temp;
		 }
		 cout << "Введите заказы покупателей>";
		 for (int i = 0; i < n; i++)
		 {
			 int temp;
			 cin >> temp;
			 seller.push_back(temp);
			 sum_seller += temp;
		 }
		 mat_mn price;

		 if (sum_buyer > sum_seller) {
			 seller.push_back(sum_buyer - sum_seller);
			 n++;
			 price = zeros(m, n);

			 for (int i = 0; i < m; i++) {
				 for (int j = 0; j < n-1; j++)
				 {
					 price[i][j] = rand() % 9 + 1;
				 }
			 }
		 }

		 else if (sum_buyer < sum_seller) {
			  buyer.push_back(sum_seller - sum_buyer);
			  m++;
			  price = zeros(m, n);
			  for (int i = 0; i < m-1; i++) {
				  for (int j = 0; j < n; j++)
				  {
					  price[i][j] = rand() % 9 + 1;
				  }
			  }
		 }
		 else {
			 price = zeros(m, n);
			 for (int i = 0; i < m; i++) {
				 for (int j = 0; j < n; j++)
				 {
					 price[i][j] = rand() % 9 + 1;
				 }
			 }
		 }
	 
	/*mat_mn price = { 
		 {1,     0 ,    3,    4 ,   2},
		 {5,     1 ,    2 ,    3,     3},
		 { 4    , 8,     1,     4 ,    3} 
	 };*/
	 
	 
	 
	 /*for (int i = 0; i < m; i++) {
		 for (int j = 0; j < n; j++)
		 {
			cout <<  price[i][j] << " ";
		 }
		 cout << endl;
	 }*/


	simplexMethodTest(price,buyer,seller, m ,n);
	//simplexMethodTest();

	return 0;
}