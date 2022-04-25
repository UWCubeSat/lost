"""
Calibrate camera, i.e. calculate camera parameters (esp. camera matrix, distortion coefficients)

Note:
  Images must be of the same pattern, taken from the same camera with unchanged focal length

Givens:
  Let dimensions of the pattern be n x m

Find 3D coordinates for all of chessboard's corners (assume z = 0)
  [
    [0, 0, 0], [0, 1, 0], ... , [0, m-1, 0]
    [1, 0, 0], [1, 1, 0], ... , [1, m-1, 0]
    ...
    [n-1, 0, 0], [n-1, 1, 0], ... , [n-1, m-1, 0]
  ]

Find the pixel coordinates (u, v) for each board corner for each different image
  Use findChessboardCorners(img, patternSize, flags)
    cv.findChessboardCorners(img, (n, m), None)

  Refine corners with cornerSubPix(img, corners, winSize, zeroZone, criteria)
    Assume (-1, -1) for zeroZone

Find camera parameters using calibrateCamera()
  cv.calibrateCamera(objectPoints, imagePoints, imageSize)
  This function returns: retVal, cameraMatrix, distortionCoeffs, rotationVecs, translationVec

Save computed parameters in JSON:
  Distortion coefficients: k1, k2, p1, p2, k3 (k for radial, p for tangential distortion)
  Camera matrix:
    [
      [ fx  0   cx ]
      [ 0   fy  cy ]
      [ 0   0    1 ]
    ]

"""

import cv2 as cv
import glob
import numpy as np
import json

WIN_SIZE = (11, 11)
TERM_CRITERIA = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 30, 0.001)

testingDir = input("Directory for calibration patterns: ")
print("Input dimensions of pattern")
# normally use: 7, corresponds to number of points per column
n = int(input("Number of rows: "))
# normally use: 6, corresponds to number of points per row
m = int(input("Number of columns: "))

v = input("View corners? (y/n): ")

imgPoints = []
allObjPoints = []

objPoints = np.zeros((n*m, 3), np.float32)
objPoints[:, :2] = np.mgrid[0:n, 0:m].T.reshape(-1, 2)

testImgs = glob.glob(testingDir + '/*.jpg')
imgSize = None

for count, image in enumerate(testImgs):
    img = cv.imread(image)
    gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

    res, corners = cv.findChessboardCorners(gray, (n, m), None)  # no flags
    if res:
        imgSize = gray.shape[::-1]
        allObjPoints.append(objPoints)
        cv.cornerSubPix(gray, corners, WIN_SIZE, (-1, -1),
                        TERM_CRITERIA)  # assume no zero zone
        imgPoints.append(corners)

    if v == 'y':
        cv.drawChessboardCorners(img, (n, m), corners, res)
        cv.imshow(f"Calibration Test Image {count}", img)
        cv.waitKey(0)
        cv.destroyWindow(f"Calibration Test Image {count}")

cv.destroyAllWindows()

allObjPoints = np.array(allObjPoints)

res, cameraMatrix, distortionCoeffs, rotVecs, tVec = cv.calibrateCamera(
    allObjPoints, imgPoints, imgSize, None, None)

print("Distortion Coefficients", distortionCoeffs.flatten())
print("Camera Matrix", cameraMatrix)
# print("Translation Vectors", tVec)  # don't need
# print("Rotation Vectors", rotVecs)  # don't need

res = {
    "distCoeffs": distortionCoeffs.flatten().tolist(),
    "cameraMatrix": cameraMatrix.tolist()
}

with open("calib.json", "w+") as f:
    json.dump(res, f, indent=4)
