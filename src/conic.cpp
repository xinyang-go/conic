#include "conic.hpp"
#include <Eigen/Dense>
#include <ceres/ceres.h>

using namespace conic;

Matrix33 conic::fitEllipse(const CMap2N& pt) {
    Matrix6N D(6, pt.cols());
    D.row(0) = pt.row(0).array() * pt.row(0).array();
    D.row(1) = pt.row(0).array() * pt.row(1).array();
    D.row(2) = pt.row(1).array() * pt.row(1).array();
    D.row(3) = pt.row(0).array();
    D.row(4) = pt.row(1).array();
    D.row(5).fill(1);
    Matrix66 C = Eigen::Matrix<double, 6, 6>::Zero();
    C(0, 2) = 2;
    C(1, 1) = -1;
    C(2, 1) = 2;
    Matrix66 S = D * D.transpose();
    Matrix66 H = S.inverse() * C;
    Eigen::EigenSolver<Matrix66> es(H);
    Matrix66 vec = es.pseudoEigenvectors();
    Matrix66 val = es.pseudoEigenvalueMatrix();
    int valMax;
    Matrix61 s = val.colwise().sum();
    s.maxCoeff(&valMax);
    Matrix61 p = vec.col(valMax);

    Matrix33 ellipse = Matrix33::Zero();
    ellipse(0, 0) = p(0, 0);
    ellipse(0, 1) = ellipse(1, 0) = p(1, 0) / 2;
    ellipse(1, 1) = p(2, 0);
    ellipse(0, 2) = ellipse(2, 0) = p(3, 0) / 2;
    ellipse(1, 2) = ellipse(2, 1) = p(4, 0) / 2;
    ellipse(2, 2) = p(5, 0);
    
    double det = ellipse.determinant();
    ellipse /= (det < 0) ? std::pow(-det, 1.0/3.0) : -std::pow(det, 1.0/3.0);

    return ellipse;
}

Matrix33 conic::perspectiveFromEllipseAndCentre(
        double radius, const Matrix33& ellipse, const CMap21& centre) {
    Matrix33 P0 = Matrix33::Identity();
    P0(2, 2) = radius;

    Eigen::EigenSolver<Matrix33> es(ellipse);
    Matrix33 val = es.pseudoEigenvalueMatrix();
    Matrix33 vec = es.pseudoEigenvectors();
    int minIndex;
    val.colwise().sum().minCoeff(&minIndex);
    int index[3] = {0, 1, 2};
    std::swap(index[minIndex], index[2]);

    Matrix33 P1;
    P1.col(0) = vec.col(index[0]) * sqrt(val(index[0], index[0]));
    P1.col(1) = vec.col(index[1]) * sqrt(val(index[1], index[1]));
    P1.col(2) = vec.col(index[2]) * sqrt(-val(index[2], index[2]));

    Matrix31 X0 = Matrix31::Zero();
    X0[2] = 1;
    X0 = P0.transpose() * X0;

    Matrix31 X1{centre[0], centre[1], 1.0};
    X1 = P1.transpose() * X1;

    Matrix33 I = Matrix33::Identity();
    I(2, 2) = -1;

    Matrix33 A = Matrix33::Identity();

    auto func = [&]<typename T>(const T * const x, T* r) -> bool {
        using MatrixT31 = Eigen::Matrix<T, 3, 1>;
        using MatrixT33 = Eigen::Matrix<T, 3, 3>;

        Eigen::Map<MatrixT33> A((T *) x);
        Eigen::Map<MatrixT33> rE(r);
        Eigen::Map<MatrixT31> rX(r + 9);

        rE = I.cast<T>() - A * I.cast<T>() * A.transpose();
        MatrixT31 X1_ = A.transpose() * X0;
        X1_ *=  X1[2] / X1_[2];
        rX = X1.cast<T>() - X1_;

        return true;
    };
    auto *cost = new ceres::AutoDiffCostFunction<decltype(func), 12, 9>(
        &func, ceres::DO_NOT_TAKE_OWNERSHIP);
    ceres::Problem problem;
    problem.AddResidualBlock(cost, nullptr, A.data());
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    
    return P1.inverse().transpose() * A.transpose() * P0.transpose();
}

Matrix34 conic::poseFromPerspective(const Matrix33& H) {
    Matrix34 T;
    double n = (H.col(0).norm() + H.col(1).norm()) / 2;
    T.col(0) = H.col(0) / H.col(0).norm();
    T.col(1) = H.col(1) / H.col(1).norm();
    T.col(2) = T.col(0).cross(T.col(1));
    T.col(3) = H.col(2) / n;
    return T;
}

Matrix31 conic::pointCoordinateCamera(
        const CMap21 &pt, const Matrix33 &H_c2w, const Matrix34 &T_w2c) {
    Matrix31 pn{pt[0], pt[1], 0.};
    Matrix31 pw = H_c2w * pn;
    pw /= pw[2];
    pw[2] = 0;
    Matrix31 pc = T_w2c.leftCols<3>() * pw + T_w2c.col(3);
    return pc;
}
