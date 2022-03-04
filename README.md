# 视觉系统圆形目标姿态求解

实现以下功能

* [X] 椭圆参数拟合
* [ ] Ransac 椭圆参数拟合
* [X] 使用图像中的一个椭圆方程和其圆心坐标求解透视变换
* [ ] 使用图像中一个椭圆方程和其直径方程求解透视变换
* [ ] 使用图像中两个椭圆方程进行求解透视变换
* [X] 分解透视变换得到目标姿态

注意事项：圆形目标姿态解不唯一，空间平面法向量可能发生镜像翻转，同时会绕圆心任意角度旋转。

---

 依赖库：Eigen3.3

---

参考文献

* Zdešar A, Škrjanc I, Klančar G. Homography estimation from circular motion for use in visual control[J]. Robotics and Autonomous Systems, 2014, 62(10): 1486-1496.
* 《计算机视觉中的多视图几何》
