#include <iostream>
#include <math.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

int main() {
	double M, N;
	Mat image = imread("chatgpt.jpg");

    cout << "Enter Target resolution (M':N') \n";
    cout << "M':";
    cin >> M;
    cout << "N':";
    cin >> N;

	Mat resample(image.rows * M, image.cols * N, CV_8UC3);

	for (int y = 0; y < resample.rows - 1; y++) {
		for (int x = 0; x < resample.cols - 1; x++) {

			int px = (int)(x / N);
			int py = (int)(y / M);

			if (px >= image.cols - 1 || py >= image.rows - 1) break;

			double fx1 = (double)x / (double)N - (double)px;
			double fx2 = 1 - fx1;
			double fy1 = (double)y / (double)M - (double)py;
			double fy2 = 1 - fy1;

			double w1 = fx2 * fy2;
			double w2 = fx1 * fy2;
			double w3 = fx2 * fy1;
			double w4 = fx1 * fy1;

			Vec3b p1 = image.at<Vec3b>(py, px);
			Vec3b p2 = image.at<Vec3b>(py, px + 1);
			Vec3b p3 = image.at<Vec3b>(py + 1, px);
			Vec3b p4 = image.at<Vec3b>(py + 1, px + 1);

			resample.at<Vec3b>(y, x) = (w1 * p1) + (w2 * p2) + (w3 * p3) + (w4 * p4);
		}
	}

	imshow("before resampling", image);
	cout << "before resampling : " << image.rows << " X " << image.cols << "\n";
	imshow("after resampling", resample);
	cout << "after resampling: " << resample.rows << " X " << resample.cols << "\n";
    
    imwrite("resampled_image.jpg", resample);
	
    waitKey(0);

	return 0;
}