#pragma once


# include "FileLAS.hpp"
# include "PointCloudAI.hpp"

void M3C2(const char* incloud1, const char* incloud2, const char* cpcloud, const char* outcloud);
void ParamEstimate(PointCloudAI &pCloud1, PointCloudAI &pCloud2, double& bestnormalScale, double& bestprojScale, double& bestprojDepth);
bool planefit(vector<Eigen::Vector3d> points, double n[4]);
void computeM3C2(int index);
void computeAlongAcross(Eigen::Vector3d invector, Eigen::Vector3d N, double& along, double& across);
int GetRealNeighbours(vector <Eigen::Vector3d>neigh_candidates,
	                                       Eigen::Vector3d corepoint, Eigen::Vector3d N, 
	                                       double radius, double projDepth, vector<double>& dist);
void GetCoreNormals(vector <Eigen::Vector3d>originalpoints, vector <Eigen::Vector3d>originalnormals, vector <Eigen::Vector3d>corepoints);
int getIndex(vector<Eigen::Vector3d> v, Eigen::Vector3d K);

template<typename T>
T stdvariance(const std::vector<T>& vec);

template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random);