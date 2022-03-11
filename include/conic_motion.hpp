#ifndef _CONIC_MOTION_HPP_
#define _CONIC_MOTION_HPP_

#include "conic.hpp"
#include <ceres/ceres.h>

namespace conic::motion {

    /*
     * @brief   运动模型损失函数
     * @details 透视变换后的坐标和运动模型计算的坐标之间的残差
     * @tparam  F 用户自定义运动方程。直接实现MotionEquation或继承于MotionEquation
     */
    template<typename F>
    struct MotionCost {
        /*
         * @param var_h 透视变换矩阵9个参数
         * @param var_m 运动模型N个参数。由用户指定
         * @param var_t 时间差(相位差)1个参数
         * @param residual 坐标残差2个值
         */
        template<typename T>
        bool operator()(const T * const var_h, 
                        const T * const var_m, 
                        const T * const var_t, 
                        T * residual) const { 
            // 计算透视变换坐标
            ConstMap33<T> H(var_h);
            Matrix31<T> pt_w = H * pt_c;
            pt_w /= pt_w[2];
            // 计算运动模型角度
            T theta = func(var_m, var_t[0] + timestamp);
            // 计算残差
            residual[0] = radius * ceres::cos(theta) - pt_w[0];
            residual[1] = radius * ceres::sin(theta) - pt_w[1];
            return true;
        }

        F func;             // 用户自定义运动方程
        double radius;      // 圆周半径
        double timestamp;   // 当前采样时间
        Matrix31d pt_c;     // 当前采样坐标(齐次坐标)


        static void build_problem(ceres::Problem &problem, double r, double t, Matrix21d pt,
                                  F func, double *var_h, double *var_m, double *var_t) {
            auto cost = new ceres::AutoDiffCostFunction<MotionCost, 2, 9, F::N, 1>(
                new MotionCost{std::move(func), r, t, {pt.x(), pt.y(), 1.0}} );
            problem.AddResidualBlock(cost, nullptr, var_h, var_m, var_t);
            for(int i=0; i<F::N; i++) {
                problem.SetParameterLowerBound(var_m, i, F::lower_bound[i]);
                problem.SetParameterUpperBound(var_m, i, F::upper_bound[i]);
            }
        }
    };

    /*
     * @brief 字符串字面量。需要C++20
     */
    template<size_t N>
    struct StringLiteral {
        constexpr StringLiteral(const char (&str)[N]) {
            std::copy_n(str, N, value);
        }

        char value[N];
    };

    /*
     * @brief 运动模型模板
     *        可直接实现或继承后实现
     * @tparam _S 运动方程名称，用于区分不同运动方程
     * @tparam _N 运动方程参数数量
     */
    template<StringLiteral _S, int _N>
    struct MotionEquation {
        static constexpr StringLiteral S = _S;
        static constexpr int N = _N;

        // 运动方程参数范围
        // to be implemented
        static const double upper_bound[N];
        static const double lower_bound[N];

        // 给定运动参数和当前时间，返回圆周角角度(弧度制)
        // to be implemented
        template<typename T>
        T operator()(const T * const var_m, const T &t) const;
    };

    /*
     * @brief   透视变换矩阵约束条件
     * @details h0.norm() == h1.norm();
     *          h0.dot(h1) == 0;
     */
    struct PerspectiveConstraintCost {
        template<typename T>
        bool operator()(const T * const var_h, T * residual) const {
            ConstMap31<T> h0(var_h);
            ConstMap31<T> h1(var_h + 3);
            residual[0] = h0.norm() - h1.norm();
            residual[1] = h0.dot(h1);
            return true;
        }


        static void build_problem(ceres::Problem &problem, double *var_h) {
            auto functor = new PerspectiveConstraintCost;
            auto cost = new ceres::AutoDiffCostFunction<PerspectiveConstraintCost, 2, 9>(functor);
            problem.AddResidualBlock(cost, nullptr, var_h);
        }
    };

    /*
     * @brief 联合优化运动参数和透视变换
     * @details 允许同时使用多组运动观测
     *          每组运动观测可以有不同的运动方程和时间(相位)偏移
     *          常用于同一圆周上有多个圆周运动的目标
     */
    struct MotionAndPerspectiveSolver {
        /*
         * @brief  添加一组运动观测
         * @tparam F 该组运动观测对应的运动方程
         * @param points     观测坐标点(归一化像素平面)
         * @param timestamps 观测时间戳
         * @param r          圆周半径
         * @param func       运动方程
         * @param var_m      运动方程参数
         * @param var_t      运动方程的时间(相位)偏移
         */
        template<typename F>
        void addGroupSamples(ConstMap2Nd points, ConstMap1Nd timestamps, double r, F &&func, 
                             double *var_m, double *var_t) {
            for(int i=0; i<points.cols(); i++) {
                MotionCost<std::remove_reference_t<F>>::build_problem(problem, r, timestamps(0, i), points.col(i), 
                        std::forward<F>(func), var_h, var_m, var_t);
            }
        }

        /*
         * @brief 进行一次优化
         * @param hold_h 是否固定透视变换参数
         */
        void solve(bool hold_h) {
            if(hold_h) {
                problem.SetParameterBlockConstant(var_h);
            } else {
                PerspectiveConstraintCost::build_problem(problem, var_h);
            }
            ceres::Solve(options, &problem, &summary);
        }

        double *var_h;

        ceres::Problem problem;
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
    };
}

#endif /* _CONIC_MOTION_HPP_ */
