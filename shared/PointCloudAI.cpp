// @ author: Cheng Li

#include <iostream>
#include "PointCloudAI.hpp"
#include <fstream>

// @ author: Cheng Li

PointCloudAI::PointCloudAI() {}

PointCloudAI::PointCloudAI(const vector<Eigen::Vector3d> Inpoints, vector<Eigen::Vector3d> bbox /*size 2*/)
{
    m_point3D = Inpoints;
    m_minbound = bbox.at(0);
    m_maxbound = bbox.at(1);
}

PointCloudAI::~PointCloudAI() {}

Eigen::Vector3d PointCloudAI::averagePointCloud()
{
    if (m_point3D.empty()) {
        return Eigen::Vector3d(0, 0, 0);
    }

    auto const count = static_cast<float>(m_point3D.size());

    return std::accumulate(m_point3D.begin(), m_point3D.end(), Eigen::Vector3d(0, 0, 0)) / count;
}

vector<Eigen::Vector3d> PointCloudAI::NormalEstimationHybrid(double searchradius, int max_nn)
{
    //centralized the data
    using namespace open3d::geometry;
    PointCloud pointcloud;

    Eigen::Vector3d avgVal = averagePointCloud();
    for (vector<Eigen::Vector3d>::iterator it = m_point3D.begin(); it != m_point3D.end(); it++) {
        pointcloud.points_.push_back(*it - avgVal);
    }

    const KDTreeSearchParam& search_param = KDTreeSearchParamHybrid(searchradius, max_nn);
    pointcloud.EstimateNormals(search_param);

    return pointcloud.normals_;
}

vector<Eigen::Vector3d> PointCloudAI::NormalEstimationKNN(int max_nn)
{
    // centralized the data
    using namespace open3d::geometry;
    PointCloud pointcloud;
    Eigen::Vector3d avgVal = averagePointCloud();
    for (vector<Eigen::Vector3d>::iterator it = m_point3D.begin(); it != m_point3D.end(); it++) {
        pointcloud.points_.push_back(*it - avgVal);
    }

    const KDTreeSearchParam& search_param = KDTreeSearchParamKNN(max_nn);
    pointcloud.EstimateNormals(search_param);

    return pointcloud.normals_;
}

vector<Eigen::Vector3d> PointCloudAI::NormalEstimationRadius(double searchradius)
{
    //centralized the data
    using namespace open3d::geometry;
    PointCloud pointcloud;
    Eigen::Vector3d avgVal = averagePointCloud();
    for (vector<Eigen::Vector3d>::iterator it = m_point3D.begin(); it != m_point3D.end(); it++) {
        pointcloud.points_.push_back(*it - avgVal);
    }

    const KDTreeSearchParam& search_param = KDTreeSearchParamRadius(searchradius);
    pointcloud.EstimateNormals(search_param);
    pointcloud.OrientNormalsToAlignWithDirection(); // align with Z oritention
    // actually it is simple
    /*const Eigen::Vector3d& orientation_reference = Eigen::Vector3d(0.0, 0.0, 1.0)
    if (normal.norm() == 0.0) {
        normal = orientation_reference;
    }
    else if (normal.dot(orientation_reference) < 0.0) {
        normal *= -1.0;
    }*/
    return pointcloud.normals_;
}

//void PointCloudAI::BuildSearchTree()
//{
//    using namespace open3d::geometry;
//    PointCloud pointcloud;
//    pointcloud.points_ = point3D;
//    kdtree.SetGeometry(pointcloud);
//}

// before calling it, call BuildSearchTree
int PointCloudAI::NeighbourSearchHybrid(Eigen::Vector3d querypoint, double searchradius, 
    std::vector<int>& indices, std::vector<double>& distance, int max_nn)
{
    return m_kdtree.SearchHybrid(querypoint, searchradius, max_nn, indices, distance);
}

int PointCloudAI::NeighbourSearchKNN(Eigen::Vector3d querypoint, 
    std::vector<int>& indices, std::vector<double>& distance, int max_nn)
{
    return m_kdtree.SearchKNN(querypoint, max_nn, indices, distance);
}

int PointCloudAI::NeighbourSearchRadius(Eigen::Vector3d querypoint, double searchradius,
    std::vector<int>& indices, std::vector<double>& distance)
{
    //using namespace open3d::geometry;
    //PointCloud pointcloud;
    //pointcloud.points_ = point3D;
    //kdtree.SetGeometry(pointcloud);

    int k= m_kdtree.SearchRadius(querypoint, searchradius, indices, distance);
    return k;
}

void print2Dvectors(vector<Eigen::Vector3d> inVectors, const char *savedfile)
{
    if (inVectors.empty()) {

        cout << "the input data is empty" << endl;
    }
    else
    {
        ofstream output_file(savedfile);
        ostream_iterator<double> output_iterator(output_file, "\t");
        for (int i = 0; i < inVectors.size(); i++)
        {
            copy(inVectors.at(i).begin(), inVectors.at(i).end(), output_iterator);
            output_file << '\n';
        }
        cout << "the data is saved successfully" << endl;
    }

    return;
}
















