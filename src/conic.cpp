#include "conic.hpp"
#include <Eigen/Dense>

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
    using namespace std::complex_literals;
    // 仿射矫正
    Matrix31 x{centre[0], centre[1], 1};
    Matrix31 l = ellipse * x;
    Matrix33 Hp_inv = Matrix33::Identity();
    Hp_inv.row(2) = l.transpose();
    Matrix33 Hp = Hp_inv.inverse();
    Matrix33 Ea = Hp.transpose() * ellipse * Hp;
    // 度量矫正
    std::complex<double> a = Ea(0, 0);
    std::complex<double> b = 2i*Ea(0,1);
    std::complex<double> c = -Ea(1,1);
    std::complex<double> n = (-b+std::sqrt(b*b-4.0*a*c))/(2.0*a);
    Matrix33 Ha = Matrix33::Identity();
    Ha(0, 0) = n.real();
    Ha(0, 1) = n.imag();
    Matrix33 Et = Ha.transpose() * Ea * Ha;
    // 平移矫正
    Matrix33 Ht = Matrix33::Identity();
    Ht(0, 2) = -Et(0, 2) / Et(0, 0);
    Ht(1, 2) = -Et(1, 2) / Et(1, 1);
    Matrix33 Es = Ht.transpose() * Et * Ht;
    // 缩放矫正
    Matrix33 Hs = Matrix33::Identity();
    Hs(0, 0) = Hs(1, 1) = sqrt(-Es(2, 2) / Es(0, 0) / radius / radius);
    // 合并透视变换矩阵
    Matrix33 H = Hp * Ha * Ht * Hs;
    
    return H;
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
