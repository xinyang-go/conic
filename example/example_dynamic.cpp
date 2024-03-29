#include <opencv2/opencv.hpp>
#include "conic_cv.hpp"
#include "conic_motion.hpp"
#include <opencv2/core/eigen.hpp>
#include <ceres/ceres.h>
#include <thread>


using namespace conic;
using namespace conic::motion;

/*
 * @brief 统计函数运行时间(ms)
 */
template<typename F, typename ...Ts>
auto timing(std::string_view name, F &&func, Ts &&...args) {
    using namespace std::chrono;
    if constexpr (std::is_same_v<void, std::result_of_t<F(Ts...)>>) {
        auto t1 = steady_clock::now();
        func(std::forward<Ts>(args)...);
        auto t2 = steady_clock::now();
        std::cout << name << ": " << duration_cast<microseconds>(t2-t1).count() / 1e3 << "ms" << std::endl;
    } else {
        auto t1 = steady_clock::now();
        auto result = func(std::forward<Ts>(args)...);
        auto t2 = steady_clock::now();
        std::cout << name << ": " << duration_cast<microseconds>(t2-t1).count() / 1e3 << "ms" << std::endl;
        return result;
    }
}

/*
 * @brief 显示缩放后的图像
 */
void resizeShow(const char* name, cv::Mat img, double scale) {
    cv::Mat im2show;
    cv::resize(img, im2show, {}, scale, scale);
    cv::imshow(name, im2show);
}

/*
 * @brief   实现运动方程
 * @details theta = -a/w * cos(w*t) + b*t
 *          0.780 < a < 1.045
 *          1.884 < w < 2.000
 *          b = 2.090 - a
 */
using CosineMotion = MotionEquation<"Cosine", 2>;
template<>
template<typename T>
T CosineMotion::operator()(const T * const var_m, const T &t) const {
    const T &a = var_m[0];
    const T &w = var_m[1];
    const T b = 2.090 - a;
    const T theta = -a/w * ceres::cos(w*t) + b*t;
    return theta;
}
template<>
const double CosineMotion::lower_bound[CosineMotion::N] = {
    0.780, 1.884,
};
template<>
const double CosineMotion::upper_bound[CosineMotion::N] = {
    1.045, 2.000,
};


void optimizePersectiveAndMotion(const std::vector<cv::Point2d> &points,
                                 const std::vector<double> &timestamps,
                                 CosineMotion func, double radius, 
                                 double *var_h, double *var_m, double *var_t,
                                 bool hold_h = false) {
    MotionAndPerspectiveSolver solver{var_h};
    // 限制优化时间不超过100ms
    // 避免数据不佳时的长时间计算
    solver.options.max_solver_time_in_seconds = 0.1;
    // 启用多线程并行
    solver.options.num_threads = std::thread::hardware_concurrency();
    solver.addGroupSamples({(const double *)points.data(), 2, (int)points.size()},
                           {timestamps.data(), 1, (int)timestamps.size()},
                           radius, func, var_m, var_t);
    solver.solve(hold_h);
}

// 此处偷懒，使用OpenCV的fitellipse计算椭圆中心
void initializePerspective(double radius, const std::vector<cv::Point2d> &points, Matrix33d &H_w2c) {
    Matrix33d ellipse = fitEllipse(points);
    std::vector<cv::Point2f> pts_f(points.begin(), points.end());
    cv::Point2d centre = cv::fitEllipse(pts_f).center;
    H_w2c = perspectiveFromEllipseAndCentre(radius, ellipse, centre);
}

int main() {
    cv::Mat K = (cv::Mat_<double>(3, 3) << 2160.99424477655, 0, 1171.17682403839, 
                    0, 2158.46880358502, 555.153465396425, 0, 0, 1);
    cv::Mat D = (cv::Mat_<double>(1, 5) << 0.0829684284954291,-0.0751893192581673, 0, 0, 0);
    cv::Mat K_inv;
    cv::invert(K, K_inv);

    cv::VideoCapture cap("cosine.mp4");
    float fps = cap.get(cv::CAP_PROP_FPS);

    std::vector<cv::Point2d> points;
    std::vector<double> timestamps;

    Matrix33d H_w2c, H_c2w;
    double motion_param[2] = {0.9125, 1.942};
    double phase_bias = 0;

    CosineMotion motion_equation;

    cv::Mat img;
    for(double current_time = 0; cap.read(img); current_time += 1.0 / fps) {
        // 识别
        cv::Mat img_gray;
        cv::cvtColor(img, img_gray, cv::COLOR_BGR2GRAY);
        cv::Mat bin;
        cv::threshold(img_gray, bin, 120, 255, cv::THRESH_BINARY);
        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(bin, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_NONE);
        cv::Point centre;
        bool found = false;
        for(int i=0; i<(int)contours.size(); i++) {
            if(contours[i].size() < 100) continue;
            if(contours[i].size() > 1000) continue;
            if(hierarchy[i][2] == -1) continue;
            cv::RotatedRect ellipse = cv::fitEllipse(contours[i]);
            double S1 = M_PI * ellipse.size.area() / 4;
            double S2 = cv::contourArea(contours[i]);
            if(0.98 < S1/S2 && S1/S2 < 1/0.98) {
                int idx = hierarchy[i][2];
                while(contours[idx].size() < 10 && idx != -1) {
                    idx = hierarchy[idx][0];
                }
                centre = cv::fitEllipse(contours[idx]).center;
                found = true;
                break;
            }
        }
        // 可视化识别结果
        cv::circle(bin, centre, 10, 0, 2);
        resizeShow("bin", bin, 0.5);

        if(!found) continue;

        // 保存识别结果
        points.emplace_back(centre);
        timestamps.emplace_back(current_time);
        // 可视化采样点
        cv::Mat sample(img.rows, img.cols, CV_8UC1, 255);
        for(auto &pt: points) {
            cv::circle(sample, pt, 3, 0, -1);
        }
        if(points.size() > 10) {
            std::vector<cv::Point2f> pts_f(points.begin(), points.end());
            cv::ellipse(sample, cv::fitEllipse(pts_f), 0);
        }            
        resizeShow("sample", sample, 0.5);

        if(points.size() == 80) {
            // 采样80点时，开始初始化
            std::vector<cv::Point2d> pts_un;
            cv::undistortPoints(points, pts_un, K, D);
            // 初始化透视变换。
            // 假定椭圆中心为透视圆心
            initializePerspective(0.2, pts_un, H_w2c);
            H_c2w = H_w2c.inverse();
            // 初始化运动参数，固定透视变换不动
            optimizePersectiveAndMotion(pts_un, timestamps, motion_equation, 0.2, 
                    H_c2w.data(), motion_param, &phase_bias, true);
        } else if(points.size() > 80) {
            // 采样大于80点时，进行联合优化
            std::vector<cv::Point2d> pts_un;
            cv::undistortPoints(points, pts_un, K, D);
            // 同时优化透视变换和运动参数
            optimizePersectiveAndMotion(pts_un, timestamps, motion_equation, 0.2, 
                    H_c2w.data(), motion_param, &phase_bias, false);
            // 可视化透视变换
            cv::Mat H;
            cv::eigen2cv(H_c2w, H);
            H = K * H * K_inv;
            cv::Mat im2show;
            cv::warpPerspective(img, im2show, H, {img.cols, img.rows});
            resizeShow("perspective", im2show, 0.5);
        } else if(points.size() > 100) {
            // 采样大于100点时，认为透视变换基本收敛，不再优化透视变换矩阵
            std::vector<cv::Point2d> pts_un;
            cv::undistortPoints(points, pts_un, K, D);
            // 固定透视变换不动
            optimizePersectiveAndMotion(pts_un, timestamps, motion_equation, 0.2, 
                    H_c2w.data(), motion_param, &phase_bias, true);
            // 可视化透视变换
            cv::Mat H;
            cv::eigen2cv(H_c2w, H);
            H = K * H * K_inv;
            cv::Mat im2show;
            cv::warpPerspective(img, im2show, H, {img.cols, img.rows});
            resizeShow("perspective", im2show, 0.5);
        } else {
            // 显示原图
            resizeShow("perspective", img, 0.5);
        }

        std::cout << "estimated motion parameter: [" << motion_param[0] << ", " << motion_param[1] << "]" << std::endl;
        std::cout << "real motion parameter: [0.782185, 1.95261]" << std::endl;
        cv::waitKey(0);
    }

    return 0;
}
