#pragma once

# include "open3d/Open3D.h"

using namespace std;


class PointCloudAI{
public:
	PointCloudAI();
	PointCloudAI(const vector<Eigen::Vector3d> Inpoints, vector<Eigen::Vector3d> bbox);
	~PointCloudAI();

	PointCloudAI(const PointCloudAI&) = default;
	PointCloudAI(PointCloudAI&&) = default;
	PointCloudAI& operator= (const PointCloudAI&) = default;
	PointCloudAI& operator= (PointCloudAI&&) = default;
	//PointCloudAI& operator=(PointCloudAI&) = default;

	vector<Eigen::Vector3d> NormalEstimationHybrid(double searchradius, int max_nn=30);
	vector<Eigen::Vector3d> NormalEstimationKNN(int max_nn=30);
	vector<Eigen::Vector3d> NormalEstimationRadius(double searchradius);
	Eigen::Vector3d averagePointCloud();

	int NeighbourSearchHybrid(Eigen::Vector3d querypoint, double searchradius,
		std::vector<int>& indices, std::vector<double>& distance, int max_nn=30);

	int NeighbourSearchKNN(Eigen::Vector3d querypoint,
		std::vector<int>& indices, std::vector<double>& distance, int max_nn = 30);

	int NeighbourSearchRadius(Eigen::Vector3d querypoint, double searchradius,
		std::vector<int>& indices, std::vector<double>& distance);
	
	//void BuildSearchTree();

public:
	vector<Eigen::Vector3d> m_point3D;
	Eigen::Vector3d m_minbound;
	Eigen::Vector3d m_maxbound;
	open3d::geometry::KDTreeFlann m_kdtree;
};

void print2Dvectors(vector<Eigen::Vector3d> inVectors, const char* savedfile);



