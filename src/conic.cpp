#include "conic.hpp"
#include <Eigen/Dense>

using namespace conic;

Matrix33d conic::fitEllipse(ConstMap2Nd pt) {
    Matrix6Nd D(6, pt.cols());
    D.row(0) = pt.row(0).array() * pt.row(0).array();
    D.row(1) = pt.row(0).array() * pt.row(1).array();
    D.row(2) = pt.row(1).array() * pt.row(1).array();
    D.row(3) = pt.row(0).array();
    D.row(4) = pt.row(1).array();
    D.row(5).fill(1);
    Matrix66d C = Matrix66d::Zero();
    C(0, 2) = 2;
    C(1, 1) = -1;
    C(2, 1) = 2;
    Matrix66d S = D * D.transpose();
    Matrix66d H = S.inverse() * C;
    Eigen::EigenSolver<Matrix66d> es(H);
    Matrix66d vec = es.pseudoEigenvectors();
    Matrix66d val = es.pseudoEigenvalueMatrix();
    int valMax;
    Matrix61d s = val.colwise().sum();
    s.maxCoeff(&valMax);
    Matrix61d p = vec.col(valMax);

    Matrix33d ellipse = Matrix33d::Zero();
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

Matrix33d conic::perspectiveFromEllipseAndCentre(
        double radius, const Matrix33d& ellipse, const Matrix21d& centre) {
    using namespace std::complex_literals;
    // 仿射矫正
    Matrix31d x{centre[0], centre[1], 1};
    Matrix31d l = ellipse * x;
    Matrix33d Hp_inv = Matrix33d::Identity();
    Hp_inv.row(2) = l.transpose();
    Matrix33d Hp = Hp_inv.inverse();
    Matrix33d Ea = Hp.transpose() * ellipse * Hp;
    // 度量矫正
    std::complex<double> a = Ea(0, 0);
    std::complex<double> b = 2i*Ea(0,1);
    std::complex<double> c = -Ea(1,1);
    std::complex<double> n = (-b+std::sqrt(b*b-4.0*a*c))/(2.0*a);
    Matrix33d Ha = Matrix33d::Identity();
    Ha(0, 0) = n.real();
    Ha(0, 1) = n.imag();
    Matrix33d Et = Ha.transpose() * Ea * Ha;
    // 平移矫正
    Matrix33d Ht = Matrix33d::Identity();
    Ht(0, 2) = -Et(0, 2) / Et(0, 0);
    Ht(1, 2) = -Et(1, 2) / Et(1, 1);
    Matrix33d Es = Ht.transpose() * Et * Ht;
    // 缩放矫正
    Matrix33d Hs = Matrix33d::Identity();
    Hs(0, 0) = Hs(1, 1) = sqrt(-Es(2, 2) / Es(0, 0) / radius / radius);
    // 合并透视变换矩阵
    Matrix33d H = Hp * Ha * Ht * Hs;
    
    return H;
}

Matrix34d conic::poseFromPerspective(const Matrix33d& H) {    
    Matrix34d T;
    T.col(0) = H.col(0);
    T.col(1) = H.col(1);
    T.col(2) = H.col(0).cross(H.col(1));

    Eigen::JacobiSVD<Matrix33d> svd(
        T.leftCols<3>(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix33d &U = svd.matrixU();
    const Matrix33d &V = svd.matrixV();
    Matrix31d S{1.0, 1.0, 1.0 / (U * V).determinant()};
    T.leftCols<3>() = U * S.asDiagonal() * V.transpose();
    
    double n = (H.col(0).norm() + H.col(1).norm()) / 2.0;
    T.col(3) = H.col(2) / n;
    return T;
}

Matrix31d conic::pointCoordinateCamera(
        const Matrix21d &pt, const Matrix33d &H_c2w, const Matrix34d &T_w2c) {
    Matrix31d pn{pt[0], pt[1], 0.};
    Matrix31d pw = H_c2w * pn;
    pw /= pw[2];
    pw[2] = 0;
    Matrix31d pc = T_w2c.leftCols<3>() * pw + T_w2c.col(3);
    return pc;
}
