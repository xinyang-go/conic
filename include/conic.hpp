#ifndef _CONIC_HPP_
#define _CONIC_HPP_

#include <Eigen/Core>

/*
 * 函数参数除透视变换矩阵和椭圆参数矩阵外，采用Map进行传参
 * 便于和OpenCV相关数据结构进行兼容
 * 如拟合椭圆时std::vector<cv::Point2d>会使用Map2N进行传参
 */

#define DEFINE_MATRIX(w, h) template<typename T>                           \
                            using Matrix##w##h = Eigen::Matrix<T, w, h>;   \
                            using Matrix##w##h##d = Matrix##w##h<double>;  \
                            using Matrix##w##h##f = Matrix##w##h<float>;   \
                            template<typename T>                           \
                            using Map##w##h = Eigen::Map<Matrix##w##h<T>>; \
                            using Map##w##h##d = Map##w##h<double>;        \
                            using Map##w##h##f = Map##w##h<float>;         \
                            template<typename T>                           \
                            using ConstMap##w##h = Eigen::Map<const Matrix##w##h<T>>; \
                            using ConstMap##w##h##d = ConstMap##w##h<double>;         \
                            using ConstMap##w##h##f = ConstMap##w##h<float>;     
#define N   Eigen::Dynamic
namespace conic {
    DEFINE_MATRIX(1, N);
    DEFINE_MATRIX(2, 1);
    DEFINE_MATRIX(2, N);
    DEFINE_MATRIX(3, 1);
    DEFINE_MATRIX(3, 3);
    DEFINE_MATRIX(3, 4);
    DEFINE_MATRIX(3, N);
    DEFINE_MATRIX(4, 1);
    DEFINE_MATRIX(4, 4);
    DEFINE_MATRIX(4, N);
    DEFINE_MATRIX(6, 1);
    DEFINE_MATRIX(6, 6);
    DEFINE_MATRIX(6, N);
}
#undef  N

namespace conic {

    /*
     * @brief  椭圆拟合
     * @param  pt 椭圆圆周采样坐标点
     * @return 二次型表示的椭圆参数，行列式被归一化为-1
     */
    Matrix33d fitEllipse(ConstMap2Nd pt);

    /*
     * @brief   从椭圆方程和圆心坐标计算透视变换
     * @details 计算得出的透视变换可能会出现镜像翻转以及绕圆心旋转的情况
     * @param   radius  目标圆在空间中的物理尺寸
     * @param   ellipse 在归一化像素平面上的椭圆参数
     * @param   centre  在归一化像素平面上的圆心坐标
     * @return  从空间平面到归一化像素平面的透视变换矩阵
     */
    Matrix33d perspectiveFromEllipseAndCentre(
            double radius, const Matrix33d& ellipse, const Matrix21d& centre);

    /*
     * @brief  从透视变换中恢复位姿
     * @param  H_w2c 从空间平面到归一化像素平面的透视变换矩阵
     * @return 从空间平面坐标系到相机坐标系的[3x4]位姿矩阵([3x3]的旋转矩阵和[3x1]的平移矩阵)
     */
    Matrix34d poseFromPerspective(const Matrix33d& H_w2c);

    /*
     * @brief  计算归一化像素平面上一点所对应相机坐标系下的坐标
     * @details 透视变换矩阵H有8个自由度，而三位姿态仅6个自由度
     *          缺少的两个自由度要求H的前两列正交且模长相等
     *          故前两列模长会被矫正成当前模长的均值
     *          方向会被矫正成在对角线两侧并且保持正交
     * @param  pt    在归一化像素平面上的坐标
     * @param  H_c2w 归一化像素平面到空间平面的透视变换矩阵
     *               (perspectiveFromEllipseAndCentre返回值的逆矩阵)
     * @param  T_w2c 空间平面坐标系到相机坐标系的位姿矩阵
     * @return 该点在相机坐标系下的坐标
     */
    Matrix31d pointCoordinateCamera(
            const Matrix21d& pt, const Matrix33d &H_c2w, const Matrix34d &T_w2c);

}




#endif /* _CONIC_HPP_ */
