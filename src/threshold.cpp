// C++ program to compute submatrix query sum in O(1)
// time
#include<iostream>
using namespace std;

// constructor(type of thresholding)
//  which preprocess
//  which query function

class Threshold {       // The class
  public:             // Access specifier
    Threshold() : image(0.0f), imageWidth(0), imageHeight(0) {};
    //Parameterized Constructor
    Threshold(unsigned char *image_, int imageWidth_, int imageHeight_)
    {
        image = image_;
        imageWidth = imageWidth_;
        imageHeight = imageHeight_;
    }
    unsigned char *image;
    int imageWidth;
    int imageHeight;

    int preProcess(int mat[imageWidth][imageHeight], int aux[imageWidth][imageHeight]);

    int sumQuery(int aux[imageWidth][imageWidth], int tli, int tlj, int rbi, int rbj);

    int BadThreshold(unsigned char *image, int imageWidth, int imageHeight);

    int LessBasicThreshold(unsigned char *image, int imageWidth, int imageHeight, int x, int y);

    int OtsusThreshold(unsigned char *image, int imageWidth, int imageHeight);

    int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight);

    int BasicThresholdOnePass(unsigned char *image, int imageWidth, int imageHeight);

};


// Function to preprocess input mat[M][N]. This function
// mainly fills aux[M][N] such that aux[i][j] stores sum
// of elements from (0,0) to (i,j)
int Threshold::preProcess(int mat[M][N], int aux[M][N])
{
  // Copy first row of mat[][] to aux[][]
  for (int i=0; i<N; i++)
    aux[0][i] = mat[0][i];

  // Do column wise sum
  for (int i=1; i<M; i++)
    for (int j=0; j<N; j++)
      aux[i][j] = mat[i][j] + aux[i-1][j];

  // Do row wise sum
  for (int i=0; i<M; i++)
    for (int j=1; j<N; j++)
      aux[i][j] += aux[i][j-1];
}

// A O(1) time function to compute sum of submatrix
// between (tli, tlj) and (rbi, rbj) using aux[][]
// which is built by the preprocess function
int Threshold::sumQuery(int aux[M][N], int tli, int tlj, int rbi, int rbj)
{
  // result is now sum of elements between (0, 0) and
  // (rbi, rbj)
  int res = aux[rbi][rbj];

  // Remove elements between (0, 0) and (tli-1, rbj)
  if (tli > 0)
  res = res - aux[tli-1][rbj];

  // Remove elements between (0, 0) and (rbi, tlj-1)
  if (tlj > 0)
  res = res - aux[rbi][tlj-1];

  // Add aux[tli-1][tlj-1] as elements between (0, 0)
  // and (tli-1, tlj-1) are subtracted twice
  if (tli > 0 && tlj > 0)
  res = res + aux[tli-1][tlj-1];

  return res;
}


// a poorly designed thresholding algorithm
int Threshold::BadThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    //loop through entire array, find sum of magnitudes
    long totalMag = 0;
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    return (((totalMag/(imageHeight * imageWidth)) + 1) * 15) / 10;
}

// a simple, but well tested thresholding algorithm that works well with star images
int Threshold::LessBasicThreshold(unsigned char *image, int imageWidth, int imageHeight, int x, int y) {
    //preprocess
    //sumQuery(img, max(0, i - k), max(0, j - k), min(i+k, m), min(j+k, n))
    unsigned long totalMag = 0;
    float std = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        sumQuery(img, max(0, i - k), max(0, j - k), min(i+k, m), min(j+k, n))
        totalMag += image[i];
    }
    float mean = totalMag / totalPixels;
    for (long i = 0; i < totalPixels; i++) {
        std += std::pow(image[i] - mean, 2);
    }
    std = std::sqrt(std / totalPixels);
    return mean + (std * 5);
}


// a more sophisticated thresholding algorithm, not tailored to star images
int Threshold::OtsusThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    // code here, duh
    long total = imageWidth * imageHeight;
    //float top = 255;
    float sumB = 0;
    float sum1 = 0;
    float wB = 0;
    float maximum = 0;
    int level = 0;
    // make the histogram (array length 256)
    int histogram[256];

    memset(histogram, 0, sizeof(int)*256);

    for (long i = 0; i < total; i++) {
        histogram[image[i]]++;
    }
    for (int i = 0; i < 256; i ++) {
        sum1 += i * histogram[i];
    }
    for (int i = 0; i < 256; i ++) {
        float wF = total - wB;
        //std::cout << "wF\n" << wB << "\n";
        //std::cout << "wB\n" << wF << "\n";
        if (wB > 0 && wF > 0) {
            float mF = (sum1 - sumB) / wF;
            float val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
            //std::cout << val << "\n";
            if (val >= maximum) {
                level = i;
                maximum = val;
            }
        }
        wB = wB + histogram[i];
        sumB = sumB + i * histogram[i];
    }
    return level;
}

// a simple, but well tested thresholding algorithm that works well with star images
int Threshold::BasicThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    unsigned long totalMag = 0;
    float std = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
    }
    float mean = totalMag / totalPixels;
    for (long i = 0; i < totalPixels; i++) {
        std += std::pow(image[i] - mean, 2);
    }
    std = std::sqrt(std / totalPixels);
    return mean + (std * 5);
}

// basic thresholding, but do it faster (trade off of some accuracy?)
int Threshold::BasicThresholdOnePass(unsigned char *image, int imageWidth, int imageHeight) {
    unsigned long totalMag = 0;
    float std = 0;
    float sq_totalMag = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
        sq_totalMag += image[i] * image[i];
    }
    float mean = totalMag / totalPixels;
    float variance = (sq_totalMag / totalPixels) - (mean * mean);
    std = std::sqrt(variance);
    return mean + (std * 5);
}


// Driver program
int main()
{
  int mat[M][N] = {{1, 2, 3, 4, 6},
            {5, 3, 8, 1, 2},
            {4, 6, 7, 5, 5},
            {2, 4, 8, 9, 4} };
  int aux[M][N];

  preProcess(mat, aux);

  int tli = 2, tlj = 2, rbi = 3, rbj = 4;
  cout << "\nQuery1: " << sumQuery(aux, tli, tlj, rbi, rbj);

  tli = 0, tlj = 0, rbi = 1, rbj = 1;
  cout << "\nQuery2: " << sumQuery(aux, tli, tlj, rbi, rbj);

  tli = 1, tlj = 2, rbi = 3, rbj = 3;
  cout << "\nQuery3: " << sumQuery(aux, tli, tlj, rbi, rbj);

  return 0;
}
