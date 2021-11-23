#include <iostream>
#include "FileLAS.hpp"
#include "PointCloudAI.hpp"
#include "M3C2.hpp"

// @ author: Cheng Li 

#include <thread>

#include <ctime>

#include <ppl.h>

using namespace concurrency;

int main()
{
    //const char* inputfile = "D:/work/Canlvl325_01_Oct2020/test.las";
    const char* inputfile1 = //"D:/work/Canlvl325_01_Oct2020/Canlvl325_01_13Jan2021.las";
        //"D:/work/python_lib_code/NZ19_Wellington.las";
        "D:/work/python_lib_code/F5_March_2019_1cm_crop_crop_small.las";

    //vector<Eigen::Vector3d> bbox(2);
    //vector<Eigen::Vector3d> point3D = ReadLas(inputfile1, bbox);
    //PointCloudAI pCloud1(point3D, bbox);
    /*cout << "the point size of the cloud: " << point3D.size() << endl;
    cout << "the x range  of box is: " << (bbox.at(1)[0] - bbox.at(0)[0]) << endl;
    cout << "the y range  of box is: " << (bbox.at(1)[1] - bbox.at(0)[1]) << endl;
    cout << "the z range  of box is: " << (bbox.at(1)[2] - bbox.at(0)[2]) << endl;*/

    const char* inputfile2 = 
        "D:/work/python_lib_code/F5_October_2019_1cm_crop_crop_register_smalla.las";

    //point3D = ReadLas(inputfile2, bbox);
    //PointCloudAI pCloud2(point3D, bbox);
   /* cout << "the point size of the cloud: " << point3D.size() << endl;
    cout << "the x range  of box is: " << (bbox.at(1)[0] - bbox.at(0)[0]) << endl;
    cout << "the y range  of box is: " << (bbox.at(1)[1] - bbox.at(0)[1]) << endl;
    cout << "the z range  of box is: " << (bbox.at(1)[2] - bbox.at(0)[2]) << endl;

    cout << "running here!" << endl;*/

    //SaveLas("result2.las", pCloud1.m_point3D);

    //vector<double> val1(pCloud1.m_point3D.size(), 1);
    //vector<int> val2(pCloud1.m_point3D.size(), 1);
    //vector<double> val3(pCloud1.m_point3D.size(), 1);

    //SaveLasExtraDims("result2.las", pCloud1.m_point3D, val1, val2, val3);

    const char* resultfile = "M3C2results2.las";
    const char* corepointfile = 
        "D:/work/python_lib_code/F5_March_2019_1cm_crop_crop_small_subsampling.las";

    // test M3C2 function
    M3C2(inputfile1, inputfile2, corepointfile, resultfile);


    /********test Las saved function and print point cloud  ***********************************************************
    // matrix.row(i) // matrix.col(j)
    //print2Dvectors(point3D, "savedpoints.txt");
    
    const char* savedfile = "savedfile.txt";
    ofstream output_file(savedfile);
    ostream_iterator<double> output_iterator(output_file, "\t");
    for (int i = 0; i < vector_2d.size(); i++)
    {
        copy(vector_2d.at(i).begin(), vector_2d.at(i).end(), output_iterator);
        output_file << '\n';
    }
    cout << "the data is saved successfully" << endl;*/

    //test normal estimation
    //vector<Eigen::Vector3d> normals = pCloud.NormalEstimationRadius(0.136153);
    //print2Dvectors(normals, R"(savednormal3.txt)");
    /************************************************************************************************/

    
    // test NeighbourSearchHybrid
    
    //int nn = 2;
    //std::vector<int> new_indices_vec(nn);
    //std::vector<double> new_dists_vec(nn);

    //// necessary steps to call NeighbourSearchHybrid
    //using namespace open3d::geometry;
    //PointCloud pointcloud;
    //pointcloud.points_ = pCloud1.m_point3D;
    //pCloud1.m_kdtree.SetGeometry(pointcloud);

    //cout << "running here" << endl;

    ////Eigen::Vector3d testpoint(pCloud1.m_point3D.at(1)[0]-0.01, pCloud1.m_point3D.at(1)[1], pCloud1.m_point3D.at(1)[2]);

    //double EffectiveSearchRadius;
    //EffectiveSearchRadius = std::sqrt(5.759067 * 5.759067 + 0.394492 / 2 * 0.394492 / 2);
    ////int k = pCloud.NeighbourSearchHybrid(pCloud.point3D.at(0), nn, new_indices_vec, new_dists_vec);
    //int k = pCloud1.NeighbourSearchRadius(pCloud1.m_point3D.at(99910), EffectiveSearchRadius, new_indices_vec, new_dists_vec);
    //
    //cout << "the number of neighbours is " << k << endl;
    //for (int i = 0; i < new_indices_vec.size(); i++) {
    //    cout << new_indices_vec[i] << endl;
    //}
    

    

    /*
    int nodata = 10;
    vector <Eigen::Vector3d> data(nodata, Eigen::Vector3d(0, 0, 0));

    vector <int> pointindex(nodata, 0);
    for (int i = 0; i < nodata; i++)
    {
        pointindex[i] = i;
    }

    std::clock_t start;
    double duration;

    start = std::clock(); // get current time

    parallel_for_each(pointindex.begin(), pointindex.end(), [](int value) {
         cout<<sumvale(value)<<endl;
        });

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

    std::cout << "Operation took " << duration << "seconds" << std::endl;
    */

   

}