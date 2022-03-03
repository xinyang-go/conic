#ifndef _CONIC_CV_HPP_
#define _CONIC_CV_HPP_

#include "conic.hpp"
#include <opencv2/core.hpp>

/*
 * 本文件用于提供OpenCV兼容API
 */

namespace conic {

    /*
     * @brief  椭圆拟合
     * @param  pt 椭圆圆周采样坐标点
     * @return 二次型表示的椭圆参数，行列式被归一化为-1
     */
    inline Matrix33 fitEllipse(const std::vector<cv::Point2d>& pt) {
        return fitEllipse(CMap2N{(const double*)pt.data(), 2, (int)pt.size()});
    }

    /*
     * @brief   从椭圆方程和圆心坐标计算透视变换
     * @details 计算得出的透视变换可能会出现镜像翻转以及绕圆心旋转的情况
     * @param   radius  目标圆在空间中的物理尺寸
     * @param   ellipse 在归一化像素平面上的椭圆参数
     * @param   centre  在归一化像素平面上的圆心坐标
     * @return  从空间平面到归一化像素平面的透视变换矩阵
     */
    inline Matrix33 perspectiveFromEllipseAndCentre(
            double radius, const Matrix33& ellipse, const cv::Point2d& centre) {
        return perspectiveFromEllipseAndCentre(radius, ellipse, CMap21{(const double*)&centre});
    }

    /*
     * @brief  计算归一化像素平面上一点所对应相机坐标系下的坐标
     * @param  pt    在归一化像素平面上的坐标
     * @param  H_c2w 从空间平面到归一化像素平面的透视变换矩阵
     * @param  T_w2c 空间平面坐标系到相机坐标系的位姿矩阵
     * @return 该点在相机坐标系下的坐标
     */
    inline Matrix31 pointCoordinateCamera(
            const cv::Point2d& pt, const Matrix33 &H_c2w, const Matrix34 &T_w2c) {
        return pointCoordinateCamera(CMap21{(const double*)&pt}, H_c2w, T_w2c);
    }
}

#endif /* _CONIC_CV_HPP_ */
