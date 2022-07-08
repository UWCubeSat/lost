#include "attitude-estimators.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace lost {

#define EPSILON 0.0001       // threshold for 0 for Newton-Raphson method

DqmQuestHelperMatrices DqmQuestHelperMatricesConstructor(const Camera &camera,
                                 const Stars &stars,
                                 const Catalog &catalog,
                                 const StarIdentifiers &starIdentifiers) {

    // attitude profile matrix
    Eigen::Matrix3f B;
    B.setZero();
    for (const StarIdentifier &s : starIdentifiers) {
        Star bStar = stars[s.starIndex];
        Vec3 bStarSpatial = camera.CameraToSpatial(bStar.position);
        Eigen::Vector3f bi;
        bi << bStarSpatial.x, bStarSpatial.y, bStarSpatial.z;

        CatalogStar rStar = catalog[s.catalogIndex];
        Eigen::Vector3f ri;
        ri << rStar.spatial.x, rStar.spatial.y, rStar.spatial.z;

        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        // Calculate matrix B = sum({ai}{bi}{ri}T)
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

    // only used for Davenport matrix
    // K =  [[[sigma], [Z[0]], [Z[1]], [Z[2]]], [[Z[0]], [S[0][0] - sigma], [S[0][1]], [S[0][2]]], [[Z[1]], [S[1][0]], [S[1][1] - sigma], [S[1][2]]], [[Z[2]], [S[2][0]], [S[2][1]], [S[2][2] - sigma]]]
    Eigen::Matrix4f K;
    K << sigma, Z(0), Z(1), Z(2),
        Z(0), S(0,0) - sigma, S(0,1), S(0,2),
        Z(1), S(1,0), S(1,1) - sigma, S(1,2),
        Z(2), S(2,0), S(2,1), S(2,2) - sigma;

    DqmQuestHelperMatrices m;
    m.K = K;
    m.S = S;
    m.Z = Z;
    m.sigma = sigma;
    return m;
}

Attitude DavenportQAlgorithm::Go(const Camera &camera,
                                 const Stars &stars,
                                 const Catalog &catalog,
                                 const StarIdentifiers &starIdentifiers) {

    // Create matrix
    Eigen::Matrix4f K = DqmQuestHelperMatricesConstructor(camera, stars, catalog, starIdentifiers).K;

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

/// return a matrix whose columns are the axes of the frame
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
    if ((int)stars.size() < 2 || (int)starIds.size() < 2) {
        return Attitude(Quaternion(1,0,0,0));
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

/**
 * Characteristic polynomial of the quest K-matrix
 * @see equation 19b of https://arc.aiaa.org/doi/pdf/10.2514/1.62549
 */
float QuestCharPoly(float x, float a, float b, float c, float d, float s) {return (pow(x,2)-a) * (pow(x,2)-b) - (c*x) + (c*s) - d;}

/**
 * Derivitive of the characteristic polynomial of the quest K-matrix
 */
float QuestCharPolyPrime(float x, float a, float b, float c) {return 4*pow(x,3) - 2*(a+b)*x - c;}

/**
 * Approximates roots of a real function using the Newton-Raphson algorithm 
 * @see https://www.geeksforgeeks.org/program-for-newton-raphson-method/
 */
float QuestEigenvalueEstimator(float guess, float a, float b, float c, float d, float s) {
    float height;
    do {
        height = QuestCharPoly(guess, a, b, c, d, s) / QuestCharPolyPrime(guess, a, b, c);
        guess -= height;
    } while (abs(height) >= EPSILON);

    return guess;
}

Attitude QuestAlgorithm::Go(const Camera &camera,
                            const Stars &stars,
                            const Catalog &catalog,
                            const StarIdentifiers &starIdentifiers) {

    // create matrix
    DqmQuestHelperMatrices m = DqmQuestHelperMatricesConstructor(camera, stars, catalog, starIdentifiers);

    // calculate coefficients for characteristic polynomial
    float delta = m.S.determinant();
    float kappa = (delta * m.S.inverse()).trace();
    float a = pow(m.sigma,2) - kappa;
    float b = pow(m.sigma,2) + (m.Z.transpose() * m.Z);
    float c = delta + (m.Z.transpose() * m.S * m.Z);
    float d = m.Z.transpose() * (m.S * m.S) * m.Z;

    // sum up the weights for initial guess of eigenvalue
    float guess = 0;
    for (const StarIdentifier &s : starIdentifiers) {
        guess += s.weight;
    }

    // Newton-Raphson method for estimating the largest eigenvalue
    float eig = QuestEigenvalueEstimator(guess, a, b, c, d, m.sigma);

    // solve for the optimal quaternion: from https://ahrs.readthedocs.io/en/latest/filters/quest.html
    float alpha = pow(eig,2) - pow(m.sigma, 2) + kappa;
    float beta = eig - m.sigma;
    float gamma = (eig + m.sigma) * alpha - delta;

    Eigen::Vector3f X = ((alpha * Eigen::Matrix3f::Identity()) + (beta * m.S) + (m.S * m.S)) * m.Z;

    Eigen::Vector4f optQuat;
    optQuat << gamma, X[0], X[1], X[2];
    optQuat /= optQuat.norm();

    return Attitude(Quaternion(optQuat[0], optQuat[1], optQuat[2], optQuat[3]));
}

}
