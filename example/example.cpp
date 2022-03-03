#include <opencv2/opencv.hpp>
#include "conic_cv.hpp"
#include <opencv2/core/eigen.hpp>

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

int main() {
    cv::Mat K = (cv::Mat_<double>(3, 3) << 2160.99424477655, 0, 1171.17682403839, 
                    0, 2158.46880358502, 555.153465396425, 0, 0, 1);
    cv::Mat D = (cv::Mat_<double>(1, 5) << 0.0829684284954291,-0.0751893192581673, 0, 0, 0);
    Eigen::Matrix3d K_e;
    cv::cv2eigen(K, K_e);

    cv::VideoCapture cap("test.mp4");
    cv::Mat im0, im2show;
    while(cap.read(im0)) {
        cv::Mat im0_gray;
        cv::cvtColor(im0, im0_gray, cv::COLOR_BGR2GRAY);
        cv::Mat bin;
        cv::threshold(im0_gray, bin, 100, 255, cv::THRESH_BINARY);
        // 识别椭圆及圆心
        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(bin, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_NONE);
        int ellipse_id = -1;
        for(int i=0; i<(int)contours.size(); i++) {
            if(contours[i].size() < 1000) continue;
            cv::RotatedRect ellipse = cv::fitEllipse(contours[i]);
            double S1 = M_PI * ellipse.size.area() / 4;
            double S2 = cv::contourArea(contours[i]);
            if(0.98 < S1/S2 && S1/S2 < 1/0.98) {
                if(ellipse_id == -1 || hierarchy[i][3] == ellipse_id) {
                    ellipse_id = i;
                }
            }
        }
        int centre_id = -1;
        for(int i=0; i<(int)contours.size(); i++) {
            if(hierarchy[i][3] == ellipse_id && contours[i].size() > 10) {
                centre_id = i;
                break;
            }
        }
        // 可视化识别结果
        cv::cvtColor(bin, bin, cv::COLOR_GRAY2BGR);
        cv::drawContours(bin, contours, ellipse_id, {0, 255, 0}, 2);
        cv::drawContours(bin, contours, centre_id, {0, 255, 0}, 2);
        resizeShow("detect", bin, 0.5);
        
        // 计算归一化平面上的圆心坐标
        std::vector<cv::Point2f> pt_c;
        std::copy(contours[centre_id].begin(), contours[centre_id].end(), std::back_inserter(pt_c));
        std::vector<cv::Point2f> pt_c_un;
        cv::undistortPoints(pt_c, pt_c_un, K, D);
        cv::RotatedRect centre = cv::fitEllipse(pt_c_un);
        // 计算归一化平面上的椭圆采样点
        std::vector<cv::Point2d> pt_e;
        std::copy(contours[ellipse_id].begin(), contours[ellipse_id].end(), std::back_inserter(pt_e));
        std::vector<cv::Point2d> pt_e_un;
        cv::undistortPoints(pt_e, pt_e_un, K, D);
        // 椭圆拟合
        auto ellipse = timing("fitEllipse", [&](){
                        return conic::fitEllipse(pt_e_un);});
        // 求解透视变换
        auto H_w2c = timing("perspectiveFromEllipseAndCentre", [&](){
                        return conic::perspectiveFromEllipseAndCentre(0.2, ellipse, centre.center);});
        // 求解姿态
        auto T_w2c = timing("poseFromPerspective", [&](){
                        return conic::poseFromPerspective(H_w2c);});

        // 可视化透视变换结果
        cv::Mat H_inv;
        cv::eigen2cv(H_w2c, H_inv);
        cv::Mat K_inv;
        cv::invert(K, K_inv);
        H_inv = K * H_inv * K_inv;
        cv::invert(H_inv, H_inv);
        cv::Mat im2show;
        cv::warpPerspective(im0, im2show, H_inv, {im0.cols, im0.rows});
        resizeShow("homography", im2show, 0.5);

        // 可视化重投影姿态
        Eigen::Matrix3d rmtx_e = T_w2c.leftCols<3>();
        Eigen::Vector3d tvec_e = T_w2c.col(3);
        cv::Mat tvec, rvec, rmtx;
        cv::eigen2cv(tvec_e, tvec);
        cv::eigen2cv(rmtx_e, rmtx);
        cv::Rodrigues(rmtx, rvec);
        cv::drawFrameAxes(im0, K, D, rvec, tvec, 0.1);
        resizeShow("coordinate", im0, 0.5);

        cv::waitKey(0);
    }
}