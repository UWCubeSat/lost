#include "attitude-estimators.hpp"

namespace lost {
    Quaternion AttitudeEstimationAlgorithm::Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &) {
        //create a vector that'll hold {bi} (Stars in our frame)
        //create a vector that'll hold {ri} (Stars in catalog frame)
        //Weight = 1 (can be changed later, in which case we want to make a vector to hold all weights {ai})
        //Calculate matrix B = sum({ai}{bi}{ri}T)
        // S = B + Transpose(B)
        //sigma = B[0][0] + B[1][1] + B[2][2]
        //Z = [[B[1][2] - B[2][1]], [B[2][0] - B[0][2]], [B[0][1] - B[1][0]]]
        //K =  [[[sigma], [Z[0]], [Z[1]], [Z[2]]], [[Z[0]], [S[0][0] - sigma], [S[0][1]], [S[0][2]]], [[Z[1]], [S[1][0]], [S[1][1] - sigma], [S[1][2]]], [[Z[2]], [S[2][0]], [S[2][1]], [S[2][2] - sigma]]]
        //Find eigenvalues of K, store the largest one as lambda 
        //The return quaternion components = eigenvector assocaited with lambda 
        return Quaternion();
    }
}