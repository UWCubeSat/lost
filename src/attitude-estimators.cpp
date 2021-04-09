#include "attitude-estimators.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace lost {

Quaternion DavenportQAlgorithm::Go(const Camera &camera,
                                   const Stars &stars,
                                   const Catalog &catalog,
                                   const StarIdentifiers &starIdentifiers) {

    //create a vector that'll hold {bi} (Stars in our frame)
    //create a vector that'll hold {ri} (Stars in catalog frame)
    Eigen::Matrix3f B;
    B.setZero();
    for (const StarIdentifier &s: starIdentifiers) {
        Star bStar = stars[s.starIndex];
        Vec3 bStarSpatial = camera.CameraToSpatial({bStar.x, bStar.y});
        Eigen::Vector3f bi;
        bi << bStarSpatial.y, bStarSpatial.x, bStarSpatial.z;

        CatalogStar rStar = catalog[s.catalogIndex];
        Eigen::Vector3f ri;
        ri << cos(rStar.raj2000)*cos(rStar.dej2000),
            sin(rStar.raj2000)*cos(rStar.dej2000),
            sin(rStar.dej2000);

        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        //Calculate matrix B = sum({ai}{bi}{ri}T)
        B += ri * bi.transpose() * s.weight;
    }
    // S = B + Transpose(B)
    Eigen::Matrix3f S = B + B.transpose();
    //sigma = B[0][0] + B[1][1] + B[2][2]
    float sigma = B.trace();
    //Z = [[B[1][2] - B[2][1]], [B[2][0] - B[0][2]], [B[0][1] - B[1][0]]]
    Eigen::Vector3f Z;
    Z << B(1,2) - B(2,1), 
        B(2,0) - B(0,2), 
        B(0,1) - B(1,0);
    //K =  [[[sigma], [Z[0]], [Z[1]], [Z[2]]], [[Z[0]], [S[0][0] - sigma], [S[0][1]], [S[0][2]]], [[Z[1]], [S[1][0]], [S[1][1] - sigma], [S[1][2]]], [[Z[2]], [S[2][0]], [S[2][1]], [S[2][2] - sigma]]]
    Eigen::Matrix4f K;
    K << sigma, Z(0), Z(1), Z(2), 
        Z(0), S(0,0) - sigma, S(0,1), S(0,2),
        Z(1), S(1,0), S(1,1) - sigma, S(1,2),
        Z(2), S(2,0), S(2,1), S(2,2) - sigma;
    //Find eigenvalues of K, store the largest one as lambda 
    //find the maximum index
    Eigen::EigenSolver<Eigen::Matrix4f> solver(K);
    Eigen::Vector4cf values = solver.eigenvalues();
    Eigen::Matrix4cf vectors = solver.eigenvectors();
    int maxIndex = 0;
    std::complex<float> maxIndexValue = values(0);
    for (int i = 1; i < values.size(); i++) {
        if (values(i).real() > maxIndexValue.real()) {
            maxIndex = i;
        }
    }
    //The return quaternion components = eigenvector assocaited with lambda 
    Eigen::Vector4cf maxAmount = vectors.col(maxIndex);
    // IMPORTANT: The matrix K is symmetric -- clearly first row and first column are equal.
    // Furthermore, S is symmetric because s_i,j = b_i,j + b_j,i and s_j,i=b_j,i + b_i,j, so the
    // bottom right 3x3 block of K is symmetric too. Thus all its eigenvalues and eigenvectors
    // are real.
    return Quaternion(maxAmount[0].real(), maxAmount[1].real(), maxAmount[2].real(), maxAmount[3].real());
}

}
