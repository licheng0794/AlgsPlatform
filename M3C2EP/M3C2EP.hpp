#pragma once

#include "M3C2.hpp"
#include <cmath>

struct scannerparams
{
    double range;
    double scan;
    double yaw;
    Eigen::Vector3d pos;
};


void M3C2EP(const char* incloud1, const char* incloud2, const char* cpcloud, const char* outcloud);

void computeM3C2EP(int index);

void GetMeanCov(Eigen::MatrixXd Cxx, Eigen::MatrixXd trfMat, Eigen::Vector3d corepoint,
    Eigen::Vector3d referPoint, vector<Eigen::Vector3d> candidate,
    vector <scannerparams> candidate_scanner, Eigen::Vector3d& mean, Eigen::Matrix3d& covariance);