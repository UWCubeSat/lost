//
// Created by bkkos on 2/17/2022.
//

// Step 1: Generate Image
// Step 2: Run the vector of Star objects through our copy-and-pasted undistortion algo
// Step 3: Use opencv's undistort and to undistort generated image and compare with ours.
//        We can compare either by comparing each pixel side by side or just eyeballing each undistorted image printed
//        out with an opencv func.

#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdio.h>
#include <iostream>

// Defining the dimensions of the checkerboard
int CHECKERBOARD[2]{6, 9};

int main() {
    // Creating the vector to store vectors of 3D points for each checkerboard image
    std::vector<std::vector<cv::Point3f>> objpoints;

    // Creating a vector to store vectors of 2D points for each checkerboard image
    std::vector<std::vector<cv::Point2f>> imgpoints;

    // Defining the world coordinates for 3D points
    std::vector<cv::Point3f> obj_pts;

    for(int i = 0; i < CHECKERBOARD[1]; i++) {
        for(int j = 0; j < CHECKERBOARD[0]; j++) {
            obj_pts.push_back(cv::Point3f(j, i, 0));
        }
    }

    // Extracting path of individual image stored in a given directory
    std::vector<cv::String> images;

    // Path of the folder containing the checkerboard images
    std::string path = "./images/*.jpg";

    // Create a glob of all possible image file names
    cv::glob(path, images);

    // Defining two opencv Matrix objects
    cv::Mat frame, gray;

    //Vector to store the pixel coordinates of detected checkerboard corners
    std:vector<cv::Point2f> corner_pts;
    bool success;

    // Looping over all the images in the directory
    for(int i = 0; i < images.size(); i++) {
        frame = cv::imread(images[i]);
        cv::cvtColor(frame, gray, cv::COLOR_BGR2GRAY);

        // Finding checkerboard corners
        // If desired num of corners are found in the image, then success = true
        success = cv::findChessboardCorners(gray, cv::Size(CHECKERBOARD[0], CHECKERBOARD[1]), corner_pts,
                                            CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK |
                                            CV_CALIB_CB_NORMALIZE_IMAGE);

        // If the desired num of corners are detected, then we refine the pixel coordinates
        // and display them on the images of the checkerboard
        if(success) {
            cv::TermCriteria criteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.001);

            // Refining pixel coordinates for given 2D points.
            cv::cornerSubPix(gray, corner_pts, cv::Size(11, 11), cv::Size(-1, -1), criteria);

            // Displaying the detected corner points on the checkerboard
            cv::drawChessboardCorners(frame, cv::Size(CHECKERBOARD[0], CHECKERBOARD[1]), corner_pts, success);

            objpoints.push_back(obj_pts);
            imgpoints.push_back(corner_pts);
        }

        cv::imshow("Image", frame);
        cv::waitKey(0);
    }

    cv::destroyAllWindows();

    cv::Mat cameraMatrix, distCoeffs, R, T;

    // Performing camera calibration by passing the value of known 3D points (objpoints) and
    // corresponding pixel coordinates of the detected corners (imgpoints)
    cv::calibrateCamera(objpoints, imgpoints, cv::Size(gray.rows, gray.cols), cameraMatrix, distCoeffs,
                        R, T);

    std::cout << "cameraMatrix: " << cameraMatrix << std::endl;
    std::cout << "distCoeffs: " << distCoeffs << std::endl;
    std::cout << "Rotation vector: " << R << std::endl;
    std::cout << "Translation vector: " << T << std::endl;

    return 0;


}




