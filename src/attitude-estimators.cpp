#include "attitude-estimators.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include "decimal.hpp"

namespace lost {

#define EPSILON (decimal) 0.0001       // threshold for 0 for Newton-Raphson method

Attitude DavenportQAlgorithm::Go(const Camera &camera,
                                 const Stars &stars,
                                 const Catalog &catalog,
                                 const StarIdentifiers &starIdentifiers) {
    if (starIdentifiers.size() < 2) {
        return Attitude();
    }
    assert(stars.size() >= 2);

    // attitude profile matrix
    #ifdef LOST_FLOAT_MODE
        Eigen::Matrix3f B;
    #else
        Eigen::Matrix3d B;
    #endif


    B.setZero();
    for (const StarIdentifier &s : starIdentifiers) {
        Star bStar = stars[s.starIndex];
        Vec3 bStarSpatial = camera.CameraToSpatial(bStar.position);

        #ifdef LOST_FLOAT_MODE
            Eigen::Vector3f bi;
            Eigen::Vector3f ri;
        #else
            Eigen::Vector3d bi;
            Eigen::Vector3d ri;
        #endif

        bi << bStarSpatial.x, bStarSpatial.y, bStarSpatial.z;

        CatalogStar rStar = catalog[s.catalogIndex];
        ri << rStar.spatial.x, rStar.spatial.y, rStar.spatial.z;
        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        // Calculate matrix B = sum({ai}{bi}{ri}T)
        B += ri * bi.transpose() * s.weight;
    }

    // S = B + Transpose(B)
    #ifdef LOST_FLOAT_MODE
        Eigen::Matrix3f S = B + B.transpose();
    #else
        Eigen::Matrix3d S = B + B.transpose();
    #endif

    //sigma = B[0][0] + B[1][1] + B[2][2]
    decimal sigma = B.trace();

    //Z = [[B[1][2] - B[2][1]], [B[2][0] - B[0][2]], [B[0][1] - B[1][0]]]
    #ifdef LOST_FLOAT_MODE
        Eigen::Vector3f Z;
    #else
        Eigen::Vector3d Z;
    #endif

    Z << B(1,2) - B(2,1),
        B(2,0) - B(0,2),
        B(0,1) - B(1,0);

    // K =  [[[sigma], [Z[0]], [Z[1]], [Z[2]]], [[Z[0]], [S[0][0] - sigma], [S[0][1]], [S[0][2]]], [[Z[1]], [S[1][0]], [S[1][1] - sigma], [S[1][2]]], [[Z[2]], [S[2][0]], [S[2][1]], [S[2][2] - sigma]]]
    #ifdef LOST_FLOAT_MODE
        Eigen::Matrix4f K;
    #else
        Eigen::Matrix4d K;
    #endif

    K << sigma, Z(0), Z(1), Z(2),
        Z(0), S(0,0) - sigma, S(0,1), S(0,2),
        Z(1), S(1,0), S(1,1) - sigma, S(1,2),
        Z(2), S(2,0), S(2,1), S(2,2) - sigma;

    //Find eigenvalues of K, store the largest one as lambda
    //find the maximum index
    #ifdef LOST_FLOAT_MODE
        Eigen::EigenSolver<Eigen::Matrix4f> solver(K);
        Eigen::Vector4cf values = solver.eigenvalues();
        Eigen::Matrix4cf vectors = solver.eigenvectors();
    #else
        Eigen::EigenSolver<Eigen::Matrix4d> solver(K);
        Eigen::Vector4cd values = solver.eigenvalues();
        Eigen::Matrix4cd vectors = solver.eigenvectors();
    #endif

    int maxIndex = 0;
    decimal maxEigenvalue = values(0).real();
    for (int i = 1; i < values.size(); i++) {
        if (values(i).real() > maxEigenvalue) {
            maxIndex = i;
            maxEigenvalue = values(i).real();
        }
    }

    //The return quaternion components = eigenvector assocaited with lambda
    #ifdef LOST_FLOAT_MODE
        Eigen::Vector4cf maxEigenvector = vectors.col(maxIndex);
    #else
        Eigen::Vector4cd maxEigenvector = vectors.col(maxIndex);
    #endif

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
    Vec3 d2 = v1.CrossProduct(v2).Normalize();
    Vec3 d3 = d1.CrossProduct(d2).Normalize();
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
    if (starIds.size() < 2) {
        return Attitude();
    }
    assert(stars.size() >= 2);

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
decimal QuestCharPoly(decimal x, decimal a, decimal b, decimal c, decimal d, decimal s) {return (DECIMAL_POW(x,2)-a) * (DECIMAL_POW(x,2)-b) - (c*x) + (c*s) - d;}

/**
 * Derivitive of the characteristic polynomial of the quest K-matrix
 */
decimal QuestCharPolyPrime(decimal x, decimal a, decimal b, decimal c) {return 4*DECIMAL_POW(x,3) - 2*(a+b)*x - c;}

/**
 * Approximates roots of a real function using the Newton-Raphson algorithm 
 * @see https://www.geeksforgeeks.org/program-for-newton-raphson-method/
 */
decimal QuestEigenvalueEstimator(decimal guess, decimal a, decimal b, decimal c, decimal d, decimal s) {
    decimal height;
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

    if (starIdentifiers.size() < 2) {
        return Attitude();
    }
    assert(stars.size() >= 2);

    // initial guess for eigenvalue (sum of the weights)
    decimal guess = 0;

    // attitude profile matrix
    Mat3 B = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (const StarIdentifier &s : starIdentifiers) {
        Star bStar = stars[s.starIndex];
        Vec3 bStarSpatial = camera.CameraToSpatial(bStar.position);

        CatalogStar rStar = catalog[s.catalogIndex];
        Vec3 rStarSpatial = {rStar.spatial.x, rStar.spatial.y, rStar.spatial.z};

        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        // Calculate matrix B = sum({ai}{bi}{ri}T)
        B = B + (rStarSpatial.OuterProduct(bStarSpatial) * s.weight);

        // sum up the weights for initial guess of eigenvalue
        guess += s.weight;
    }

    // S = B + Transpose(B)
    Mat3 S = B + B.Transpose();
    //sigma = B[0][0] + B[1][1] + B[2][2]
    decimal sigma = B.Trace();
    //Z = [[B[1][2] - B[2][1]], [B[2][0] - B[0][2]], [B[0][1] - B[1][0]]]
    Vec3 Z = {
        B.At(1,2) - B.At(2,1),
        B.At(2,0) - B.At(0,2),
        B.At(0,1) - B.At(1,0)
    };

    // calculate coefficients for characteristic polynomial
    decimal delta = S.Det();
    decimal kappa = (S.Inverse() * delta).Trace();
    decimal a = DECIMAL_POW(sigma,2) - kappa;
    decimal b = DECIMAL_POW(sigma,2) + (Z * Z);
    decimal c = delta + (Z * S * Z);
    decimal d = Z * (S * S) * Z;

    // Newton-Raphson method for estimating the largest eigenvalue
    decimal eig = QuestEigenvalueEstimator(guess, a, b, c, d, sigma);

    // solve for the optimal quaternion: from https://ahrs.readthedocs.io/en/latest/filters/quest.html
    decimal alpha = DECIMAL_POW(eig,2) - DECIMAL_POW(sigma, 2) + kappa;
    decimal beta = eig - sigma;
    decimal gamma = (eig + sigma) * alpha - delta;

    Vec3 X = ((kIdentityMat3 * alpha) + (S * beta) + (S * S)) * Z;
    decimal scalar = 1 / DECIMAL_SQRT(DECIMAL_POW(gamma,2) + X.MagnitudeSq());
    X = X * scalar;
    gamma *= scalar;

    return Attitude(Quaternion(gamma, X.x, X.y, X.z));
}

}
