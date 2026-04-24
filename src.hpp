#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
// #define IGNORE_MATRIX 0
#define IGNORE_MATRIX 1

#if IGNORE_MATRIX

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

public:

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) : m(m_), n(n_) {
        if (m <= 0 || n <= 0) {
            throw matrix_error();
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    // 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) : m(obj.m), n(obj.n) {
        if (obj.data == nullptr) {
            data = nullptr;
            return;
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
    }

    // 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept : m(obj.m), n(obj.n), data(obj.data) {
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        
        // 释放当前数据
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
        
        m = obj.m;
        n = obj.n;
        
        if (obj.data == nullptr) {
            data = nullptr;
            return *this;
        }
        
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
        return *this;
    }

    // 重载括号，返回矩阵的第i行(1-based)、第j列(1-based)的元素的引用。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 1 || j > n) {
            throw matrix_error();
        }
        return data[i-1][j-1];
    }

    // 重载乘号，返回矩阵乘法 lhs * rhs 的结果。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                result.data[i][j] = fraction(0);
                for (int k = 0; k < lhs.n; k++) {
                    result.data[i][j] = result.data[i][j] + lhs.data[i][k] * rhs.data[k][j];
                }
            }
        }
        return result;
    }

    // 返回矩阵的转置。
    matrix transposition() {
        if (data == nullptr || m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // 返回矩阵的行列式。建议用高斯消元实现。
    fraction determination() {
        if (data == nullptr || m == 0 || n == 0 || m != n) {
            throw matrix_error();
        }
        
        // 创建副本进行高斯消元
        matrix temp(*this);
        fraction det(1);
        
        for (int i = 0; i < m; i++) {
            // 找到主元
            int pivot = i;
            for (int j = i + 1; j < m; j++) {
                // 简单比较：找非零元素
                if (temp.data[j][i].operator==(fraction(0)) == false) {
                    if (temp.data[pivot][i].operator==(fraction(0))) {
                        pivot = j;
                    }
                }
            }
            
            // 如果主元为0，行列式为0
            if (temp.data[pivot][i].operator==(fraction(0))) {
                return fraction(0);
            }
            
            // 交换行
            if (pivot != i) {
                for (int j = 0; j < m; j++) {
                    fraction tmp = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot][j];
                    temp.data[pivot][j] = tmp;
                }
                det = det * fraction(-1);
            }
            
            // 消元
            det = det * temp.data[i][i];
            for (int j = i + 1; j < m; j++) {
                fraction factor = temp.data[j][i] / temp.data[i][i];
                for (int k = i; k < m; k++) {
                    temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                }
            }
        }
        
        return det;
    }
    
    // 辅助函数：获取行数和列数
    int rows() const { return m; }
    int cols() const { return n; }
};

#endif

class resistive_network {
private:
    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 邻接矩阵A，电导矩阵C，Laplace矩阵(A^tCA)
    matrix adjacency, conduction, laplace;
    
    // 存储连接信息
    int *from_nodes, *to_nodes;
    fraction *resistances;
    
    // 求解线性方程组 Ax = b (使用高斯消元法)
    void solve_linear_system(matrix &A, fraction *b, fraction *x, int n) {
        // 创建增广矩阵
        matrix aug(n, n + 1);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                aug(i, j) = A(i, j);
            }
            aug(i, n + 1) = b[i - 1];
        }
        
        // 高斯消元
        for (int i = 1; i <= n; i++) {
            // 找主元
            int pivot = i;
            for (int j = i + 1; j <= n; j++) {
                if (!(aug(j, i) == fraction(0))) {
                    if (aug(pivot, i) == fraction(0)) {
                        pivot = j;
                    }
                }
            }
            
            // 交换行
            if (pivot != i) {
                for (int j = 1; j <= n + 1; j++) {
                    fraction tmp = aug(i, j);
                    aug(i, j) = aug(pivot, j);
                    aug(pivot, j) = tmp;
                }
            }
            
            // 消元
            for (int j = i + 1; j <= n; j++) {
                if (!(aug(i, i) == fraction(0))) {
                    fraction factor = aug(j, i) / aug(i, i);
                    for (int k = i; k <= n + 1; k++) {
                        aug(j, k) = aug(j, k) - factor * aug(i, k);
                    }
                }
            }
        }
        
        // 回代
        for (int i = n; i >= 1; i--) {
            x[i - 1] = aug(i, n + 1);
            for (int j = i + 1; j <= n; j++) {
                x[i - 1] = x[i - 1] - aug(i, j) * x[j - 1];
            }
            if (!(aug(i, i) == fraction(0))) {
                x[i - 1] = x[i - 1] / aug(i, i);
            }
        }
    }

public:

    // 设置电阻网络
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
        : interface_size(interface_size_), connection_size(connection_size_),
          adjacency(connection_size_, interface_size_),
          conduction(connection_size_, connection_size_) {
        
        // 保存连接信息
        from_nodes = new int[connection_size];
        to_nodes = new int[connection_size];
        resistances = new fraction[connection_size];
        
        for (int i = 0; i < connection_size; i++) {
            from_nodes[i] = from[i];
            to_nodes[i] = to[i];
            resistances[i] = resistance[i];
        }
        
        // 构建邻接矩阵 A (m x n)
        // A[k][i] = 1 if edge k starts from node i, -1 if ends at node i, 0 otherwise
        for (int k = 1; k <= connection_size; k++) {
            for (int i = 1; i <= interface_size; i++) {
                if (from[k-1] == i) {
                    adjacency(k, i) = fraction(1);
                } else if (to[k-1] == i) {
                    adjacency(k, i) = fraction(-1);
                } else {
                    adjacency(k, i) = fraction(0);
                }
            }
        }
        
        // 构建电导矩阵 C (m x m, diagonal)
        for (int i = 1; i <= connection_size; i++) {
            for (int j = 1; j <= connection_size; j++) {
                if (i == j) {
                    conduction(i, j) = fraction(1) / resistance[i-1];
                } else {
                    conduction(i, j) = fraction(0);
                }
            }
        }
        
        // 计算 Laplace 矩阵 L = A^T * C * A
        matrix temp = conduction * adjacency;
        laplace = adjacency.transposition() * temp;
    }

    ~resistive_network() {
        delete[] from_nodes;
        delete[] to_nodes;
        delete[] resistances;
    }

    // 返回节点 interface_id1 和 interface_id2 之间的等效电阻
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) {
            return fraction(0);
        }
        
        // 构造电流向量：在 id1 注入 1A，在 id2 流出 1A
        fraction *current = new fraction[interface_size];
        for (int i = 0; i < interface_size; i++) {
            current[i] = fraction(0);
        }
        current[interface_id1 - 1] = fraction(1);
        current[interface_id2 - 1] = fraction(-1);
        
        // 求解电压 (去掉最后一个节点的方程，设 U_n = 0)
        int n = interface_size - 1;
        matrix L_reduced(n, n);
        fraction *I_reduced = new fraction[n];
        fraction *U_reduced = new fraction[n];
        
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                L_reduced(i, j) = laplace(i, j);
            }
            I_reduced[i - 1] = current[i - 1];
        }
        
        solve_linear_system(L_reduced, I_reduced, U_reduced, n);
        
        // 计算电压差
        fraction u1 = (interface_id1 == interface_size) ? fraction(0) : U_reduced[interface_id1 - 1];
        fraction u2 = (interface_id2 == interface_size) ? fraction(0) : U_reduced[interface_id2 - 1];
        
        fraction result = u1 - u2;
        
        delete[] current;
        delete[] I_reduced;
        delete[] U_reduced;
        
        return result;
    }

    // 在给定节点电流I的前提下，返回节点id的电压
    fraction get_voltage(int id, fraction current[]) {
        // 求解电压 (去掉最后一个节点的方程，设 U_n = 0)
        int n = interface_size - 1;
        matrix L_reduced(n, n);
        fraction *I_reduced = new fraction[n];
        fraction *U_reduced = new fraction[n];
        
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                L_reduced(i, j) = laplace(i, j);
            }
            I_reduced[i - 1] = current[i - 1];
        }
        
        solve_linear_system(L_reduced, I_reduced, U_reduced, n);
        
        fraction result = (id == interface_size) ? fraction(0) : U_reduced[id - 1];
        
        delete[] I_reduced;
        delete[] U_reduced;
        
        return result;
    }

    // 在给定节点电压U的前提下，返回电阻网络的功率
    fraction get_power(fraction voltage[]) {
        // P = sum over all edges of (U_from - U_to)^2 / R
        fraction power(0);
        
        for (int k = 0; k < connection_size; k++) {
            fraction u_from = voltage[from_nodes[k] - 1];
            fraction u_to = voltage[to_nodes[k] - 1];
            fraction voltage_diff = u_from - u_to;
            power = power + (voltage_diff * voltage_diff) / resistances[k];
        }
        
        return power;
    }
};


#endif //SRC_HPP