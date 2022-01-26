#include <iostream>
#include "CSF.h"
#include "../shared/FileLAS.hpp"
#include "../shared/PointCloudAI.hpp"
#include <ctime>

// note: I did not use dll and lib since they did not work may because of the version compatibility
// I directly add source code of http://ramm.bnu.edu.cn/projects/CSF/download/#source-code

void main(int argc, char* argv[])
{

    //step 1 read point cloud

    if (argc != 2)
    {
        cout << "the number of the input parameters is not 1!" << endl;
        cout << "the correct form is 'CSF inputcloud' (other parameters are fixed)" << endl;
        return;
    }

    const char* inputfile1 = argv[1];

    if (!is_file_exist(inputfile1))
    {
        cout << "the input file not found!" << endl;
        return;
    }

    CSF csf;
    cout << "The CSF for ground filter is running!" << endl;
    vector<Eigen::Vector3d> bbox(2);
    vector<Eigen::Vector3d> point3D = ReadLas(inputfile1, bbox);

    csf::PointCloud pointcloud;

    for (int i = 0; i < point3D.size(); ++i)
    {
        csf::Point point;
        point.x = point3D[i][0];
        point.y = point3D[i][1];
        point.z = point3D[i][2];

        pointcloud.push_back(point);
    }

    csf.setPointCloud(pointcloud);

    //step 2 parameter settings
    //Among these paramters:  
    //time_step  interations class_threshold can remain as defualt in most cases.
    csf.params.bSloopSmooth = true;
    csf.params.cloth_resolution = 2;
    csf.params.rigidness = 3;

    csf.params.time_step = 0.65;
    csf.params.class_threshold = 0.5;
    csf.params.interations = 500;

    //step 3 do filtering
    //result stores the index of ground points or non-ground points in the original point cloud
    std::clock_t start;
    double duration;
    start = std::clock(); // get current time
    
    std::vector<int> GroundIndexes, offGroundIndexes;
    csf.do_filtering(GroundIndexes, offGroundIndexes);

    vector<Eigen::Vector3d>  groundlas;
    for (int i = 0; i < GroundIndexes.size(); i++)
    {
        groundlas.push_back(point3D[GroundIndexes[i]]);
    }

    SaveLas("ground.las", groundlas);

    vector<Eigen::Vector3d>  nongroundlas;
    for (int i = 0; i < offGroundIndexes.size(); i++)
    {
        nongroundlas.push_back(point3D[offGroundIndexes[i]]);
    }
    SaveLas("non-ground.las", nongroundlas);


    cout << "the ground precentage is " << GroundIndexes.size()*100.0 / (GroundIndexes.size() + offGroundIndexes.size()) << "%" << endl;

    cout << "successfully save results in 'ground.las' and 'non-ground.las'" << endl;

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "Operation took " << duration << "seconds" << endl;

}

