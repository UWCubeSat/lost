#include "attitude-estimators.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace lost {
    Quaternion DavenportQAlgorithm::Go(const Camera &cameraBoy, const Stars &starBoy, const Catalog &catalogBoy, const StarIdentifiers &StarIdentifiersBoy) {
        //create a vector that'll hold {bi} (Stars in our frame)
        //create a vector that'll hold {ri} (Stars in catalog frame)
        std::vector<Eigen::Vector3f> b;
        std::vector<Eigen::Vector3f> r;
        for (const StarIdentifier &s: StarIdentifiersBoy) {
            Star bStar = starBoy[s.starIndex];
            float ra;
            float dej;
            cameraBoy.CoordinateAngles({bStar.x, bStar.y}, &ra, &dej);
            Eigen::Vector3f bi;
            bi << cos(ra)*cos(dej),
                sin(ra)*cos(dej),
                sin(dej);
            r.push_back(bi);

            CatalogStar rStar = catalogBoy[s.catalogIndex];
            Eigen::Vector3f ri;
            ri << cos(rStar.raj2000)*cos(rStar.dej2000),
                sin(rStar.raj2000)*cos(rStar.dej2000),
                sin(rStar.dej2000);
            r.push_back(ri);
        }
        
        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        //Calculate matrix B = sum({ai}{bi}{ri}T)
        Eigen::Matrix3f B;
        for (int i = 0; i < b.size(); i++) {
            B += b[i] * r[i].transpose();
        }
        // S = B + Transpose(B)
        Eigen::Matrix3f S = B + B.transpose();
        //sigma = B[0][0] + B[1][1] + B[2][2]
        float sigma = B(0,0) + B(1,1) + B(2,2);
        //Z = [[B[1][2] - B[2][1]], [B[2][0] - B[0][2]], [B[0][1] - B[1][0]]]
        Eigen::Matrix3f Z;
        Z << B(1,2) - B(2,1), 
            B(2,0) - B(0,2), 
            B(0,1) - B(1,0);
        //K =  [[[sigma], [Z[0]], [Z[1]], [Z[2]]], [[Z[0]], [S[0][0] - sigma], [S[0][1]], [S[0][2]]], [[Z[1]], [S[1][0]], [S[1][1] - sigma], [S[1][2]]], [[Z[2]], [S[2][0]], [S[2][1]], [S[2][2] - sigma]]]
        Eigen::Matrix3f K;
        K << sigma, Z(0), Z(1), Z(2), 
            Z(0), S(0,0) - sigma, S(0,1), S(0,2),
            Z(1), S(1,0), S(1,1) - sigma, S(1,2),
            Z(2), S(2,0), S(2,1), S(2,2) - sigma;
        //Find eigenvalues of K, store the largest one as lambda 
        Eigen::Vector3cf eigens = K.eigenvalues();
        //find the maximum index
        Eigen::EigenSolver<Eigen::Matrix3f> solver(K);
        Eigen::Vector3cf values = solver.eigenvalues();
        Eigen::Matrix3cf vectors = solver.eigenvectors();
        int maxIndex = 0;
        std::complex<float> maxIndexValue = values(0);
        for (int i = 1; i < values.size(); i++) {
            if (abs(values(i)) > abs(maxIndexValue)) {
                maxIndex = i;
            }
        }
        //The return quaternion components = eigenvector assocaited with lambda 
        auto maxAmount = vectors.col(maxIndex);
        return Quaternion(abs(maxAmount(0)), abs(maxAmount(1)), abs(maxAmount(2)), abs(maxAmount(3)));
    }
}