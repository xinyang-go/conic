#ifndef _CONIC_HPP_
#define _CONIC_HPP_

#include <Eigen/Core>

/*
 * 函数参数除透视变换矩阵和椭圆参数矩阵外，采用Map进行传参
 * 便于和OpenCV相关数据结构进行兼容
 * 如拟合椭圆时std::vector<cv::Point2d>会使用Map2N进行传参
 */

namespace conic {

    using Matrix21 = Eigen::Matrix<double, 2, 1>;
    using Matrix2N = Eigen::Matrix<double, 2, Eigen::Dynamic>;
    using Matrix31 = Eigen::Matrix<double, 3, 1>;
    using Matrix33 = Eigen::Matrix<double, 3, 3>;
    using Matrix34 = Eigen::Matrix<double, 3, 4>;
    using Matrix41 = Eigen::Matrix<double, 4, 1>;
    using Matrix44 = Eigen::Matrix<double, 4, 4>;
    using Matrix61 = Eigen::Matrix<double, 6, 1>;
    using Matrix66 = Eigen::Matrix<double, 6, 6>;
    using Matrix6N = Eigen::Matrix<double, 6, Eigen::Dynamic>;

    using Map21 = Eigen::Map<Matrix21>;
    using Map2N = Eigen::Map<Matrix2N>;
    using Map31 = Eigen::Map<Matrix31>;
    using Map33 = Eigen::Map<Matrix33>;
    using Map34 = Eigen::Map<Matrix34>;
    using Map41 = Eigen::Map<Matrix41>;
    using Map44 = Eigen::Map<Matrix44>;
    using Map61 = Eigen::Map<Matrix61>;
    using Map66 = Eigen::Map<Matrix66>;
    using Map6N = Eigen::Map<Matrix6N>;
    
    using ConstMap21 = Eigen::Map<const Matrix21>;
    using ConstMap2N = Eigen::Map<const Matrix2N>;
    using ConstMap31 = Eigen::Map<const Matrix31>;
    using ConstMap33 = Eigen::Map<const Matrix33>;
    using ConstMap34 = Eigen::Map<const Matrix34>;
    using ConstMap41 = Eigen::Map<const Matrix41>;
    using ConstMap44 = Eigen::Map<const Matrix44>;
    using ConstMap61 = Eigen::Map<const Matrix61>;
    using ConstMap66 = Eigen::Map<const Matrix66>;
    using ConstMap6N = Eigen::Map<const Matrix6N>;

    /*
     * @brief  椭圆拟合
     * @param  pt 椭圆圆周采样坐标点
     * @return 二次型表示的椭圆参数，行列式被归一化为-1
     */
    Matrix33 fitEllipse(ConstMap2N pt);

    /*
     * @brief   从椭圆方程和圆心坐标计算透视变换
     * @details 计算得出的透视变换可能会出现镜像翻转以及绕圆心旋转的情况
     * @param   radius  目标圆在空间中的物理尺寸
     * @param   ellipse 在归一化像素平面上的椭圆参数
     * @param   centre  在归一化像素平面上的圆心坐标
     * @return  从空间平面到归一化像素平面的透视变换矩阵
     */
    Matrix33 perspectiveFromEllipseAndCentre(
            double radius, const Matrix33& ellipse, const Matrix21& centre);

    /*
     * @brief  从透视变换中恢复位姿
     * @param  H_w2c 从空间平面到归一化像素平面的透视变换矩阵
     * @return 从空间平面坐标系到相机坐标系的[3x4]位姿矩阵([3x3]的旋转矩阵和[3x1]的平移矩阵)
     */
    Matrix34 poseFromPerspective(const Matrix33& H_w2c);

    /*
     * @brief  计算归一化像素平面上一点所对应相机坐标系下的坐标
     * @param  pt    在归一化像素平面上的坐标
     * @param  H_c2w 从空间平面到归一化像素平面的透视变换矩阵
     * @param  T_w2c 空间平面坐标系到相机坐标系的位姿矩阵
     * @return 该点在相机坐标系下的坐标
     */
    Matrix31 pointCoordinateCamera(
            const Matrix21& pt, const Matrix33 &H_c2w, const Matrix34 &T_w2c);

}




#endif /* _CONIC_HPP_ */
