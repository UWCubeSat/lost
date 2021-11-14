#include "attitude-estimators.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace lost {

Attitude DavenportQAlgorithm::Go(const Camera &camera,
                                 const Stars &stars,
                                 const Catalog &catalog,
                                 const StarIdentifiers &starIdentifiers) {

    //create a vector that'll hold {bi} (Stars in our frame)
    //create a vector that'll hold {ri} (Stars in catalog frame)
    Eigen::Matrix3f B;
    B.setZero();
    for (const StarIdentifier &s: starIdentifiers) {
        Star bStar = stars[s.starIndex];
        Vec3 bStarSpatial = camera.CameraToSpatial(bStar.position);
        Eigen::Vector3f bi;
        bi << bStarSpatial.y, bStarSpatial.x, bStarSpatial.z;

        CatalogStar rStar = catalog[s.catalogIndex];
        Eigen::Vector3f ri;
        ri << rStar.spatial.x, rStar.spatial.y, rStar.spatial.z;

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
    float maxEigenvalue = values(0).real();
    for (int i = 1; i < values.size(); i++) {
        if (values(i).real() > maxEigenvalue) {
            maxIndex = i;
            maxEigenvalue = values(i).real();
        }
    }
    //The return quaternion components = eigenvector assocaited with lambda 
    Eigen::Vector4cf maxEigenvector = vectors.col(maxIndex);
    // IMPORTANT: The matrix K is symmetric -- clearly first row and first column are equal.
    // Furthermore, S is symmetric because s_i,j = b_i,j + b_j,i and s_j,i=b_j,i + b_i,j, so the
    // bottom right 3x3 block of K is symmetric too. Thus all its eigenvalues and eigenvectors
    // are real.
    return Attitude(Quaternion(maxEigenvector[0].real(),
                               maxEigenvector[1].real(),
                               maxEigenvector[2].real(),
                               maxEigenvector[3].real()));
}

// return a matrix whose columns are the axes of the frame
static Mat3 TriadCoordinateFrame(Vec3 v1, Vec3 v2) {
    Vec3 d1 = v1.Normalize();
    Vec3 d2 = v1.crossProduct(v2).Normalize();
    Vec3 d3 = d1.crossProduct(d2).Normalize();
    return {
        d1.x, d2.x, d3.x,
        d1.y, d2.y, d3.y,
        d1.z, d2.z, d3.z,
    };
}

Attitude TriadAlgorithm::Go(const Camera &camera,
                            const Stars &stars,
                            const Catalog &catalog,
                            const StarIdentifiers &starIds) {
    if ((int)stars.size() < 2) {
        return Quaternion(1,0,0,0);
    }

    // TODO: Better way of picking the two stars
    StarIdentifier
        a = starIds[0],
        b = starIds[starIds.size()/2];

    Mat3 photoFrame = TriadCoordinateFrame(camera.CameraToSpatial(stars[a.starIndex].position),
                                           camera.CameraToSpatial(stars[b.starIndex].position));
    Mat3 catalogFrame = TriadCoordinateFrame(catalog[a.catalogIndex].spatial,
                                             catalog[b.catalogIndex].spatial);

    // Now, we calculate the rotation that would bring the catalog vector to the photo vector. We
    // have the relationship photoFrame = attitude*catalogFrame, and since catalogFrame is an
    // orthogonal matrix, we have attitude = photoFrame*transpase(catalogFrame)
    return Attitude(photoFrame*catalogFrame.Transpose());
}

}
