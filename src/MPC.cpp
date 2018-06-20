#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// 设置预测状态点的个数和两个相邻状态点之间的时间间隔
size_t N = 20;
double dt = 0.05; // 50ms

// 设置汽车头到车辆重心之间的距离
const double Lf = 2.67;

// 参考速度，为了避免车辆在行驶中停止
double ref_v = 65;

// solver使用的是一个向量存储所有的状态值，所以需要确定每种状态在向量中的开始位置
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
public:

	// 参考路径方程的系数
	Eigen::VectorXd coeffs;
	FG_eval(Eigen::VectorXd coeffs) {
		this->coeffs = coeffs;
	}

	typedef CPPAD_TESTVECTOR(AD<double>)ADvector;

	// 该函数的目的是定义约束，fg向量包含的是总损失和约束，vars向量包含的是状态值和驱动器的输入
	void operator()(ADvector& fg, const ADvector& vars) {
		// 任何cost都会加到fg[0]
		fg[0] = 0;
		// 使车辆轨迹和参考路径的误差最小，且使车辆的速度尽量接近参考速度
		for (int t = 0; t < N; t++) {
			fg[0] += CppAD::pow(vars[cte_start + t], 2);
			fg[0] += CppAD::pow(vars[epsi_start + t], 2);
			fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
		}
		// 使车辆行驶更平稳，尽量减少每一次驱动器的输入大小
		// 注意：驱动器的输入是N-1个.最后一个点没有驱动器输入
		for (int t = 0; t < N - 1; t++) {
			fg[0] += CppAD::pow(vars[delta_start + t], 2);
			fg[0] += CppAD::pow(vars[a_start + t], 2);
		}
		// 为了使车辆运动更平滑，尽量减少相邻两次驱动器输入的差距
		// 注意：这个地方是N-2个
		for (int t = 0; t < N - 2; t++) {
			fg[0] += 600*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
			fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
		}
		for (int t = 0; t < N - 2; t++) {
			fg[0] += 600*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
			fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
		}

		// 设置fg的初始值为状态的初始值，这个地方为初始条件约束
		fg[1 + x_start] = vars[x_start];
		fg[1 + y_start] = vars[y_start];
		fg[1 + psi_start] = vars[psi_start];
		fg[1 + v_start] = vars[v_start];
		fg[1 + cte_start] = vars[cte_start];
		fg[1 + epsi_start] = vars[epsi_start];

		// 因为t=0初始条件约束已经有了，计算其他约束条件
		for (int t = 1; t < N; t++) {
			// t+1时刻的状态
			AD<double> x1 = vars[x_start + t];
			AD<double> y1 = vars[y_start + t];
			AD<double> psi1 = vars[psi_start + t];
			AD<double> v1 = vars[v_start + t];
			AD<double> cte1 = vars[cte_start + t];
			AD<double> epsi1 = vars[epsi_start + t];

			// t时刻的状态
			AD<double> x0 = vars[x_start + t - 1];
			AD<double> y0 = vars[y_start + t - 1];
			AD<double> psi0 = vars[psi_start + t - 1];
			AD<double> v0 = vars[v_start + t - 1];
			AD<double> cte0 = vars[cte_start + t - 1];
			AD<double> epsi0 = vars[epsi_start + t - 1];
			// t时刻的驱动器输入
			AD<double> delta0 = vars[delta_start + t - 1];
			AD<double> a0 = vars[a_start + t - 1];

			// t时刻参考路径的距离和角度值
			AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2]*x0*x0 + coeffs[3]*x0*x0*x0;
			AD<double> psides0 = CppAD::atan(coeffs[1]+2*coeffs[2]*x0 + 3 * coeffs[3]*x0*x0);

			// 根据如上的状态值计算约束条件
			fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
			fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
			fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
			fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
			fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
			fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
		}
	}
};

// MPC class definition
MPC::MPC() {
}
MPC::~MPC() {
}

vector<vector<double>> MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
	size_t i;
	typedef CPPAD_TESTVECTOR(double)Dvector;

	double x = x0[0];
	double y = x0[1];
	double psi = x0[2];
	double v = x0[3];
	double cte = x0[4];
	double epsi = x0[5];

	// 独立状态的个数，注意：驱动器的输入（N - 1）* 2
	size_t n_vars = N * 6 + (N - 1) * 2;
	// 约束条件的个数
	size_t n_constraints = N * 6;

	// 除了初始值，初始化每一个状态为0
	Dvector vars(n_vars);
	for (int i = 0; i < n_vars; i++) {
		vars[i] = 0.0;
	}
	vars[x_start] = x;
	vars[y_start] = y;
	vars[psi_start] = psi;
	vars[v_start] = v;
	vars[cte_start] = cte;
	vars[epsi_start] = epsi;

	// 设置每一个状态变量的最大和最小值
	// 【1】设置非驱动输入的最大和最小值
	// 【2】设置方向盘转动角度范围-25—25度
	// 【3】加速度的范围-1—1
	Dvector vars_lowerbound(n_vars);
	Dvector vars_upperbound(n_vars);
	for (int i = 0; i < delta_start; i++) {
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}
	for (int i = delta_start; i < a_start; i++) {
		vars_lowerbound[i] = -0.436332;
		vars_upperbound[i] = 0.436332;
	}
	for (int i = a_start; i < n_vars; i++) {
		vars_lowerbound[i] = -1.0;
		vars_upperbound[i] = 1.0;
	}
	// 前一次计算的驱动器输入作为约束
	vars_lowerbound[a_start] = a_prev;
	vars_upperbound[a_start] = a_prev;
	vars_lowerbound[delta_start] = delta_prev;
	vars_upperbound[delta_start] = delta_prev;

	// 设置约束条件的的最大和最小值，除了初始状态其他约束都为0
	Dvector constraints_lowerbound(n_constraints);
	Dvector constraints_upperbound(n_constraints);
	for (int i = 0; i < n_constraints; i++) {
		constraints_lowerbound[i] = 0;
		constraints_upperbound[i] = 0;
	}
	constraints_lowerbound[x_start] = x;
	constraints_lowerbound[y_start] = y;
	constraints_lowerbound[psi_start] = psi;
	constraints_lowerbound[v_start] = v;
	constraints_lowerbound[cte_start] = cte;
	constraints_lowerbound[epsi_start] = epsi;

	constraints_upperbound[x_start] = x;
	constraints_upperbound[y_start] = y;
	constraints_upperbound[psi_start] = psi;
	constraints_upperbound[v_start] = v;
	constraints_upperbound[cte_start] = cte;
	constraints_upperbound[epsi_start] = epsi;

	// Object that computes objective and constraints
	FG_eval fg_eval(coeffs);

	// options
	std::string options;
	options += "Integer print_level  0\n";
	options += "Sparse  true        forward\n";
	options += "Sparse  true        reverse\n";

	// 计算这个问题
	CppAD::ipopt::solve_result<Dvector> solution;
	CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound,
			vars_upperbound, constraints_lowerbound, constraints_upperbound,
			fg_eval, solution);

	//
	// Check some of the solution values
	//
	bool ok = true;
	ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

	auto cost = solution.obj_value;
	std::cout << "Cost " << cost << std::endl;

	// 返回预测出的轨迹点状态
	vector<double> X;
	vector<double> Y;
	vector<double> Delta;
	vector<double> A;
	for (auto i = 0; i < N - 1; i++) {
		X.push_back(solution.x[x_start + i]);
		Y.push_back(solution.x[y_start + i]);
		Delta.push_back(solution.x[delta_start + i]);
		A.push_back(solution.x[a_start + i]);
	}
	vector<vector<double>> result;
	result.push_back(X);
	result.push_back(Y);
	result.push_back(Delta);
	result.push_back(A);
	return result;
}

