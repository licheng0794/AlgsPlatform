# include "M3C2EP.hpp"
# define IF_DEBUG false
# define PI 3.141592653589793238
# define Chi 7.814727903251179

struct M3C2EPparams
{
    vector<double> M3C2EPdists; // distances
    vector<double> M3C2EPsigchg; // significant change (0,1)
    vector<double> M3C2EPlod; // level of detection

    PointCloudAI* pCloud1;
    PointCloudAI* pCloud2;

    vector<Eigen::Vector3d> corepointsnomrals;
    vector<Eigen::Vector3d> corepoints;

    int maxneighbours = 100000;
    int minneighbours = 4; // Lague et.al 2013 preferably 30 

    double projDepth = 5.759067;
    double normalScale = 0.394492;
    double projScale = 0.394492;

    vector <scannerparams> Ref_scanners; // scanner parameters for the reference point cloud
    vector <scannerparams> Chg_scanners; // scanner parameters for the changed point cloud

    Eigen::MatrixXd Cxx; // using resize to set 12*12 matrix
    Eigen::MatrixXd TrfMat; // using resize to set 3*4 matrix without the last row (0,0,0,1)

    vector<int> RPtSrcId; // the point source id of reference point cloud
    vector<int> CPtSrcId; // the point source id of changed point cloud'

    Eigen::Vector3d referPoint; // the reference point obtained from ICP or default as the mean of the reference point cloud
};

static M3C2EPparams g_M3C2EPparams;

//set scanner parameters. 
//If productization, those parameters may be the inputs from UI

double VZ_2000_sigmaRange = 0.005;
double VZ_2000_sigmaScan = 0.00027 / 4;
double VZ_2000_sigmaYaw = 0.00027 / 4;

double VZ_4000_sigmaRange = 0.005;
double VZ_4000_sigmaScan = 0.0003 / 4;
double VZ_4000_sigmaYaw = 0.0003 / 4;

// the scanners of the reference point cloud
scannerparams RS1 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652831.6805, 5189073.5765, 2523.7454)};
scannerparams RS2 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652812.0904, 5189246.1069, 2433.7296) };
scannerparams RS3 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652992.6490, 5189116.7022, 2526.2710) };
scannerparams RS4 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652917.8332, 5189284.5585, 2423.7433) };
scannerparams RS5 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652804.5869, 5189190.5423, 2456.0348) };
scannerparams RS6 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652862.9167, 5189292.7994, 2403.6955) };
scannerparams RS7 = { VZ_4000_sigmaRange, VZ_4000_sigmaScan, VZ_4000_sigmaYaw,
                           Eigen::Vector3d(652949.7085, 5189182.2139, 2473.9986) };

// the scanners of the changed point cloud
scannerparams CS1 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652831.6805, 5189073.5765, 2523.7454) };
scannerparams CS2 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652812.0904, 5189246.1069, 2433.7296) };
scannerparams CS3 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652992.6490, 5189116.7022, 2526.2710) };
scannerparams CS4 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652917.8332, 5189284.5585, 2423.7433) };
scannerparams CS5 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652804.5869, 5189190.5423, 2456.0348) };
scannerparams CS6 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652862.9167, 5189292.7994, 2403.6955) };
scannerparams CS7 = { VZ_2000_sigmaRange, VZ_2000_sigmaScan, VZ_2000_sigmaYaw,
                           Eigen::Vector3d(652949.7085, 5189182.2139, 2473.9986) };


void M3C2EP(const char* incloud1, const char* incloud2, const char* cpcloud, const char* outcloud)
{
    
    cout << "M3C2 Error propagation (M3C2-EP) running..." << endl;

    cout << "Configure scanner parameters and covariance matrix ... " << endl;

    g_M3C2EPparams.Ref_scanners.push_back(RS1);
    g_M3C2EPparams.Ref_scanners.push_back(RS2);
    g_M3C2EPparams.Ref_scanners.push_back(RS3);
    g_M3C2EPparams.Ref_scanners.push_back(RS4);
    g_M3C2EPparams.Ref_scanners.push_back(RS5);
    g_M3C2EPparams.Ref_scanners.push_back(RS6);
    g_M3C2EPparams.Ref_scanners.push_back(RS7);

    g_M3C2EPparams.Chg_scanners.push_back(CS1);
    g_M3C2EPparams.Chg_scanners.push_back(CS2);
    g_M3C2EPparams.Chg_scanners.push_back(CS3);
    g_M3C2EPparams.Chg_scanners.push_back(CS4);
    g_M3C2EPparams.Chg_scanners.push_back(CS5);
    g_M3C2EPparams.Chg_scanners.push_back(CS6);
    g_M3C2EPparams.Chg_scanners.push_back(CS7);

    g_M3C2EPparams.Cxx.resize(12, 12);
    g_M3C2EPparams.TrfMat.resize(3, 4);

    // TrfMat and Cxx returned by ICP with covariance
    g_M3C2EPparams.TrfMat << 9.99848e-01, 1.48736e-04, 4.47580e-04, 1.18134e-01,
        -4.22126e-07, 9.99854e-01, 1.00573e-04, -6.15470e-02,
        5.44780e-05, 2.76962e-05, 9.99969e-01, 2.09701e-02;
    //Cxx << ;

    g_M3C2EPparams.Cxx << 8.27177e-10, -7.89891e-10, -2.08458e-09, -9.92585e-11,
        1.66127e-10, 3.88368e-10, 2.22160e-10, -3.11890e-10,
        -7.02886e-10, -5.46437e-07, 9.54442e-08, -1.83869e-07,
        -7.89891e-10, 9.50467e-10, 2.31625e-09, 1.26183e-10,
        -2.25226e-10, -5.60243e-10, -3.06693e-10, 4.80385e-10,
        1.07621e-09, 6.14116e-07, -1.28964e-07, 2.85140e-07,
        -2.08458e-09, 2.31625e-09, 6.01543e-09, 3.81831e-10,
        -6.01301e-10, -1.51542e-09, -7.11688e-10, 1.09350e-09,
        2.53973e-09, 1.55156e-06, -3.52932e-07, 6.55429e-07,
        -9.92585e-11, 1.26183e-10, 3.81831e-10, 1.21197e-09,
        -1.27871e-09, -3.29420e-09, -3.86765e-10, 3.80609e-10,
        1.03348e-09, 8.62728e-08, -8.49571e-07, 2.61477e-07,
        1.66127e-10, -2.25226e-10, -6.01301e-10, -1.27871e-09,
        1.59973e-09, 3.94055e-09, 3.84302e-10, -4.56767e-10,
        -1.15722e-09, -1.45875e-07, 1.01226e-06, -2.93230e-07,
        3.88368e-10, -5.60243e-10, -1.51542e-09, -3.29420e-09,
        3.94055e-09, 1.02388e-08, 1.04914e-09, -1.18262e-09,
        -3.15855e-09, -3.56655e-07, 2.54210e-06, -7.79504e-07,
        2.22160e-10, -3.06693e-10, -7.11688e-10, -3.86765e-10,
        3.84302e-10, 1.04914e-09, 3.18742e-10, -3.38928e-10,
        -8.69457e-10, -1.84707e-07, 2.61646e-07, -2.25126e-07,
        -3.11890e-10, 4.80385e-10, 1.09350e-09, 3.80609e-10,
        -4.56767e-10, -1.18262e-09, -3.38928e-10, 4.49591e-10,
        1.07453e-09, 2.90409e-07, -2.89552e-07, 2.81245e-07,
        -7.02886e-10, 1.07621e-09, 2.53973e-09, 1.03348e-09,
        -1.15722e-09, -3.15855e-09, -8.69457e-10, 1.07453e-09,
        2.71223e-09, 6.55167e-07, -7.58913e-07, 6.92575e-07,
        -5.46437e-07, 6.14116e-07, 1.55156e-06, 8.62728e-08,
        -1.45875e-07, -3.56655e-07, -1.84707e-07, 2.90409e-07,
        6.55167e-07, 4.08394e-04, -8.40420e-05, 1.72451e-04,
        9.54442e-08, -1.28964e-07, -3.52932e-07, -8.49571e-07,
        1.01226e-06, 2.54210e-06, 2.61646e-07, -2.89552e-07,
        -7.58913e-07, -8.40420e-05, 6.53577e-04, -1.91743e-04,
        -1.83869e-07, 2.85140e-07, 6.55429e-07, 2.61477e-07,
        -2.93230e-07, -7.79504e-07, -2.25126e-07, 2.81245e-07,
        6.92575e-07, 1.72451e-04, -1.91743e-04, 1.80198e-04;


    vector<Eigen::Vector3d> bbox(2);
    vector<int> Rsrcid;
    cout << "1. Reading the reference point cloud!" << endl;
    vector<Eigen::VectorXd> MetaData;
    vector<ExtraDim> extraDims;
    vector<Eigen::Vector3d> point3D = ReadLas(incloud1, bbox, MetaData, extraDims);
    // save PointSourceId 
    for (int i = 0; i < MetaData.size(); i++)
    {
        Rsrcid.push_back(MetaData[i][8]);
    }
    

    //vector<Eigen::Vector3d> point3D = ReadLasPtSrc(incloud1, bbox, Rsrcid);

    size_t sz = Rsrcid.size();
    // Calculate the mean
    double mean = std::accumulate(Rsrcid.begin(), Rsrcid.end(), 0.0) / sz;
    if (mean == 0)
    {
        cout << "The reference point cloud in M3C2-EP must include the dimension 'PointSourceId' !" << endl;
        cout << "Stop here !" << endl;
        return;
    }

    g_M3C2EPparams.RPtSrcId = Rsrcid;
    g_M3C2EPparams.pCloud1 = new PointCloudAI(point3D, bbox);
    g_M3C2EPparams.referPoint = g_M3C2EPparams.pCloud1->averagePointCloud(); // default
    cout << "   the reference point cloud has " << point3D.size() << " points!" << endl;

    cout << "2. Reading the target point cloud!" << endl;

    vector<int> Csrcid;
    vector<Eigen::VectorXd> MetaDataTarget;
    vector<ExtraDim> extraDimsTarget;
    point3D = ReadLas(incloud1, bbox, MetaDataTarget, extraDimsTarget);
    // save PointSourceId 
    for (int i = 0; i < MetaDataTarget.size(); i++)
    {
        Csrcid.push_back(MetaDataTarget[i][8]);
    }

    sz = Csrcid.size();
    // Calculate the mean
    mean = std::accumulate(Csrcid.begin(), Csrcid.end(), 0.0) / sz;
    if (mean == 0)
    {
        cout << "The target point cloud in M3C2-EP must include the dimension 'PointSourceId' !" << endl;
        cout << "Stop here !" << endl;
        return;
    }

    g_M3C2EPparams.CPtSrcId = Csrcid;

    for (auto it = point3D.begin(); it != point3D.end(); ++it)
    {
        
        Eigen::Vector3d diff =  (*it) - g_M3C2EPparams.referPoint;
        Eigen::MatrixXd diffconvert;
        // covert vector3d to matrixXd
        diffconvert = Eigen::MatrixXd::Map(diff.data(),1, diff.size());      

        Eigen::MatrixXd multipval = diffconvert * g_M3C2EPparams.TrfMat.block<3, 3>(0, 0);
        // convert Matrix to vector3d
        Vector3d B(Map<Vector3d>(multipval.data(), multipval.rows() * multipval.cols()));

        Eigen::MatrixXd trfMat3 = g_M3C2EPparams.TrfMat.block<3, 1>(0, 3);
        Vector3d C(Map<Vector3d>(trfMat3.data(), trfMat3.rows() * trfMat3.cols())); 

        *it = B + C + g_M3C2EPparams.referPoint;
    }

    g_M3C2EPparams.pCloud2 = new PointCloudAI(point3D, bbox);

    cout << "   the target point cloud has " << point3D.size() << " points!" << endl;

    using namespace open3d::geometry;
    PointCloud pointcloud1;
    pointcloud1.points_ = g_M3C2EPparams.pCloud1->m_point3D;
    cout << "3. Building the kdtree of the reference point cloud!" << endl;
    g_M3C2EPparams.pCloud1->m_kdtree.SetGeometry(pointcloud1);
    PointCloud pointcloud2;
    pointcloud2.points_ = g_M3C2EPparams.pCloud2->m_point3D;

    cout << "4. Building the kdtree of the target point cloud!" << endl;
    g_M3C2EPparams.pCloud2->m_kdtree.SetGeometry(pointcloud2);

    // parameter estimation, return projectionscale, normalscale, projDepth
    double normalScale = 0;
    double projScale = 0;
    double projDepth = 0;

    cout << "5. Estimating proper parameters ..." << endl;
    ParamEstimate(*g_M3C2EPparams.pCloud1, *g_M3C2EPparams.pCloud2, normalScale, projScale, projDepth);

    //g_M3C2EPparams.corepoints = g_M3C2EPparams.pCloud1->m_point3D; // temporary
    g_M3C2EPparams.normalScale = normalScale;
    g_M3C2EPparams.projScale = projScale;
    g_M3C2EPparams.projDepth = projDepth;

    /*g_M3C2EPparams.normalScale = 0.828492;
    g_M3C2EPparams.projScale = 0.207123;
    g_M3C2EPparams.projDepth = 7.45075;*/

    cout << "6. Reading the core points!" << endl;

    vector<Eigen::VectorXd> MetaDataCore;
    vector<ExtraDim> extraDimsCore;

    if (strcmp(incloud1, cpcloud) == 0)
    {
        g_M3C2EPparams.corepoints = g_M3C2EPparams.pCloud1->m_point3D;
        MetaDataCore = MetaData;
        extraDimsCore = extraDims;
        cout << "   the core point is the same with the reference point cloud!" << endl;
    }
    else
    {
        point3D = ReadLas(cpcloud, bbox, MetaDataCore, extraDimsCore);
        g_M3C2EPparams.corepoints = point3D;
    }

    int cplen = g_M3C2EPparams.corepoints.size();
    cout << "   the core point cloud has " << cplen << " points!" << endl;


    cout << "7. Computing normal of core points!" << endl;
    // here we get normals for the whole point cloud !!!
    vector<Eigen::Vector3d> normals = g_M3C2EPparams.pCloud1->NormalEstimationRadius(g_M3C2EPparams.normalScale / 2);


    g_M3C2EPparams.corepointsnomrals = normals;

    vector <int> corepointindex(cplen, 0);
    for (int i = 0; i < cplen; i++)
    {
        corepointindex[i] = i;
    }

    g_M3C2EPparams.M3C2EPdists = vector<double>(cplen, NAN);
    g_M3C2EPparams.M3C2EPsigchg = vector<double>(cplen, NAN);
    g_M3C2EPparams.M3C2EPlod = vector<double>(cplen, NAN);

    cout << "8. Computing M3C2-EP distances ... wait..." << endl;
    cout << "   the estimated waiting time is about " << cplen / 10000 * 25 << "s (not accurate)" << endl;

    if (IF_DEBUG)
    {
        std::clock_t start;
        double duration;
        start = std::clock(); // get current time
        // 10000 core points take 18.887s
        // 20000 core points take 48.211s
        //computeM3C2EP(419);
        //computeM3C2EP(419);
        for (int index = 0; index < corepointindex.size(); index++)
        {
            computeM3C2EP(index);
            // duration 0.001s for each
        }
        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "Operation took " << duration << "seconds" << std::endl;

    }
    else
    {
        std::clock_t start;
        double duration;
        start = std::clock(); // get current time
        // 10000 core points take 4.658s
        // 20000 core points takes 12.874s
////pragma omp parallel for schedule(static) \
//        num_threads(8)
//        for (int index = 0; index < corepointindex.size(); index++) {
//            computeM3C2(index);
//        }
        parallel_for_each(begin(corepointindex), end(corepointindex), [](int index) {
            computeM3C2EP(index);
            });

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "Operation took " << duration << "seconds" << std::endl;
    }

    // save results to a las file;
    if (SaveLasExtraDims(outcloud, g_M3C2EPparams.corepoints, MetaDataCore, extraDimsCore, g_M3C2EPparams.M3C2EPdists, g_M3C2EPparams.M3C2EPsigchg, g_M3C2EPparams.M3C2EPlod))
    {
        cout << "Saving running results successfully!" << endl;
    }

}

void computeM3C2EP(int index)
{
    //index = 363;
    //cout << "the current point index: " << index << endl;
   
    // get normal
    //cout << "the current core point index is " << index << endl;
    Eigen::Vector3d N = g_M3C2EPparams.corepointsnomrals[index];
    Eigen::Vector3d corepoint = g_M3C2EPparams.corepoints[index];

    std::vector<int> neigh_candidate_indices(g_M3C2EPparams.maxneighbours);
    std::vector<double> neigh_candidate_dists(g_M3C2EPparams.maxneighbours);

    double EffectiveSearchRadius;
    EffectiveSearchRadius = std::sqrt(g_M3C2EPparams.projDepth * g_M3C2EPparams.projDepth + g_M3C2EPparams.projScale / 2 * g_M3C2EPparams.projScale / 2);

    // for the refence point cloud pCloud1
    // get neighbour candidates
    int k = g_M3C2EPparams.pCloud1->NeighbourSearchRadius(g_M3C2EPparams.corepoints[index], EffectiveSearchRadius,
        neigh_candidate_indices, neigh_candidate_dists);

    if (k < g_M3C2EPparams.minneighbours)
    {
        return; // not able to compute statistcs
    }

    vector<Eigen::Vector3d> candidate;
    for (int i = 0; i < k; i++)
    {
        candidate.push_back(g_M3C2EPparams.pCloud1->m_point3D[neigh_candidate_indices[i]]);
    }

    vector<double> projDist1; // the distances of points projected to the normal direction
    // the function will output neigh_candidate_indices and projDist
    int n1 = GetRealNeighbours(candidate, neigh_candidate_indices, corepoint, N, g_M3C2EPparams.projScale / 2, g_M3C2EPparams.projDepth, projDist1);
    //cout << "n1 is: " << n1 << endl;
    if (n1 < g_M3C2EPparams.minneighbours)
    {
        return; // not able to compute statistcs
    }

    // compute mean and covariance for the reference point cloud
    Eigen::MatrixXd Cxx = g_M3C2EPparams.Cxx;
    Eigen::MatrixXd trfMat = g_M3C2EPparams.TrfMat;
    vector <scannerparams> scanners = g_M3C2EPparams.Ref_scanners;
    Eigen::Vector3d referPoint = g_M3C2EPparams.referPoint;

    vector <scannerparams> candidate_scanners;

    for (int i = 0; i < neigh_candidate_indices.size(); i++)
        candidate_scanners.push_back(scanners[g_M3C2EPparams.RPtSrcId[neigh_candidate_indices[i]]-1]);

    Eigen::Vector3d local_mean_R;
    Eigen::Matrix3d local_covariance_R;
    GetMeanCov(Cxx, trfMat, corepoint, referPoint, candidate, candidate_scanners, local_mean_R, local_covariance_R);
    //if (index == 419)
    //{
    //    cout << "local_mean_R: " << local_mean_R << endl;
    //    cout << "local_covariance_R: " << local_covariance_R << endl;
    //}
    

    // for the target point cloud pCloud2
    // get neighbour candidates
    k = g_M3C2EPparams.pCloud2->NeighbourSearchRadius(g_M3C2EPparams.corepoints[index], EffectiveSearchRadius,
        neigh_candidate_indices, neigh_candidate_dists);

    if (k < g_M3C2EPparams.minneighbours)
    {
        return; // not able to compute statistcs
    }

    candidate.clear();
    for (int i = 0; i < k; i++)
    {
        candidate.push_back(g_M3C2EPparams.pCloud2->m_point3D[neigh_candidate_indices[i]]);
    }

    vector<double> projDist2; // the distances of points projected to the normal direction
    // the function will output neigh_candidate_indices and projDist
    int n2 = GetRealNeighbours(candidate, neigh_candidate_indices, corepoint, N, g_M3C2EPparams.projScale / 2, g_M3C2EPparams.projDepth, projDist2);
    //cout << "n2 is: " << n2 << endl;
    
    if (n2 < g_M3C2EPparams.minneighbours)
    {
        return; // not able to compute statistcs
    }

    // compute mean and covariance for the target point cloud
    scanners = g_M3C2EPparams.Chg_scanners;
    
    candidate_scanners.clear();
    for (int i = 0; i < neigh_candidate_indices.size(); i++)
    {
        //cout << "the neighbour index is " << neigh_candidate_indices[i] << endl;
        candidate_scanners.push_back(scanners[g_M3C2EPparams.CPtSrcId[neigh_candidate_indices[i]]-1]);
    }

    Eigen::Vector3d local_mean_C;
    Eigen::Matrix3d local_covariance_C;

    GetMeanCov(Cxx, trfMat, corepoint, referPoint, candidate, candidate_scanners, local_mean_C, local_covariance_C);

    //if (index == 419)
    //{
    //    cout << "local_mean_C: " << local_mean_C << endl;
    //    cout << "local_covariance_C: " << local_covariance_C << endl;
    //}

    /*Eigen::MatrixXd local_p1_p2_covariance;
    local_p1_p2_covariance.resize(6, 6);
    local_p1_p2_covariance = Eigen::Matrix<double, 6, 6>::Zero();
    local_p1_p2_covariance.block<3, 3>(0, 0) = local_covariance_R;
    local_p1_p2_covariance.block<3, 3>(3, 3) = local_covariance_C;*/

   
    Eigen::Matrix3d SigmaD = local_covariance_R + local_covariance_C;

    Eigen::MatrixXd Normal;
    Normal = Eigen::MatrixXd::Map(N.data(), 3, 1);

    Eigen::MatrixXd Tsqalt = Normal.transpose()* SigmaD.inverse()* Normal;
    double T = Tsqalt(0, 0);
    //cout << T << endl;

   /* Eigen::MatrixXd normalstack;
    normalstack.resize(6, 1);
    normalstack(0, 0) = -1 * N[0];
    normalstack(0, 1) = -1 * N[1];
    normalstack(0, 2) = -1 * N[2];
    normalstack(0, 3) = N[0];
    normalstack(0, 4) = N[1];
    normalstack(0, 5) = N[2];*/


    double dist = N.dot(local_mean_R - local_mean_C);
    double lod = std::sqrt(Chi / T);

    g_M3C2EPparams.M3C2EPlod[index] = lod;
    g_M3C2EPparams.M3C2EPdists[index] = dist;

    //cout << lod << endl;
    //cout << dist << endl;

    if (dist<-lod || dist > lod) // significant change
    {
        g_M3C2EPparams.M3C2EPsigchg[index] = 1;
    }
}

 // compute covariance and mean for a given core point
void GetMeanCov(Eigen::MatrixXd Cxx, Eigen::MatrixXd trfMat, Eigen::Vector3d corepoint, 
    Eigen::Vector3d referPoint, vector<Eigen::Vector3d> candidate, 
    vector <scannerparams> candidate_scanner, Eigen::Vector3d& mean, Eigen::Matrix3d& covariance)
{
    //scanners
    const int npt = candidate.size(); // the number of candidate neighbours

    double C11=0, C12=0, C13=0, C22=0, C23=0, C33 = 0;

    Eigen::Matrix3d local_Cxx;

    Eigen::MatrixXd ATP;
    ATP.resize(3, 3 * npt); // 3* (3*npt)
    ATP = Eigen::ArrayXXd::Zero(3, 3 * npt); // zero matrix


    Eigen::MatrixXd A = Eigen::Matrix<double, 3, 3>::Identity();
    A = A.replicate(npt, 1);

    for (int i = 0; i < npt; i++)
    {

        double dx, dy, dz, rrange, sinscan, cosscan, sinyaw, cosyaw, sigmaRange, sigmaYaw, sigmaScan;

        Eigen::Vector3d dist = candidate[i] - candidate_scanner[i].pos;
        double dlx = dist[0];
        double dly = dist[1];
        double dlz = dist[2];
        double yaw = atan2(dly, dlx);
        double planar_dist = std::sqrt(dlx*dlx + dly*dly);
        double scan = PI / 2 - atan(dlz / planar_dist);

        rrange = std::sqrt(planar_dist * planar_dist + dlz * dlz);
        sinscan = sin(scan);
        cosscan = cos(scan);
        sinyaw = sin(yaw);
        cosyaw = cos(yaw);

        Eigen::Vector3d dr = candidate[i] - referPoint;
        dx = dr[0];
        dy = dr[1];
        dz = dr[2];

        sigmaRange = std::sqrt(candidate_scanner[i].range * candidate_scanner[i].range 
                        + candidate_scanner[i].range * 1e-6 * rrange * rrange);
        sigmaYaw = candidate_scanner[i].scan;
        sigmaScan = candidate_scanner[i].yaw;

        double SigmaXiXj = dx * dx * Cxx(0,0)
            + 2 * dx * dy * Cxx(0, 1)
            + dy * dy * Cxx(1, 1)
            + 2 * dy * dz * Cxx(1, 2)
            + dz * dz * Cxx(2, 2)
            + 2 * dz * dx * Cxx(0, 2)
            + 2 * (dx * Cxx(0, 9) + dy * Cxx(1, 9) + dz * Cxx(2, 9)) + Cxx(9, 9);

        double SigmaYiYj = dx * dx * Cxx(3, 3) +
            2 * dx * dy * Cxx(3, 4) +
            dy * dy * Cxx(4, 4) +
            2 * dy * dz * Cxx(4, 5) +
            dz * dz * Cxx(5, 5) +
            2 * dz * dx * Cxx(3, 5) +
            2 * (dx * Cxx(3, 10) +
                dy * Cxx(4, 10) +
                dz * Cxx(5, 10)) +
            Cxx(10, 10);

        double SigmaZiZj = dx * dx * Cxx(6, 6) +
            2 * dx * dy * Cxx(6, 7) +
            dy * dy * Cxx(7, 7) +
            2 * dy * dz * Cxx(7, 8) +
            dz * dz * Cxx(8, 8) +
            2 * dz * dx * Cxx(6, 8) +
            2 * (dx * Cxx(6, 11) +
                dy * Cxx(7, 11) +
                dz * Cxx(8, 11)) +
            Cxx(11, 11);

        double SigmaXiYj = Cxx(9, 10) +
            dx * Cxx(0, 10) +
            dy * Cxx(1, 10) +
            dz * Cxx(2, 10) +
            dx * (Cxx(3, 9) +
                Cxx(0, 3)* dx +
                Cxx(1, 3)* dy +
                Cxx(2, 3)* dz) +
            dy * (Cxx(4, 9) +
                Cxx(0, 4)* dx +
                Cxx(1, 4)* dy +
                Cxx(2, 4)* dz) +
            dz * (Cxx(5, 9) +
                Cxx(0, 5)* dx +
                Cxx(1, 5)* dy +
                Cxx(2, 5)* dz);

        double SigmaXiZj = Cxx(9, 11) +
            dx * Cxx(0, 11) +
            dy * Cxx(1, 11) +
            dz * Cxx(2, 11) +
            dx * (Cxx(6, 9) +
                Cxx(0, 6) * dx +
                Cxx(1, 6) * dy +
                Cxx(2, 6) * dz) +
            dy * (Cxx(7, 9) +
                Cxx(0, 7) * dx +
                Cxx(1, 7) * dy +
                Cxx(2, 7) * dz) +
            dz * (Cxx(8, 9) +
                Cxx(0, 8) * dx +
                Cxx(1, 8) * dy +
                Cxx(2, 8) * dz);

        double SigmaYiZj = Cxx(10, 11) +
            dx * Cxx(6, 10) +
            dy * Cxx(7, 10) +
            dz * Cxx(8, 10) +
            dx * (Cxx(3, 11) +
                Cxx(3, 6) * dx +
                Cxx(3, 7) * dy +
                Cxx(3, 8) * dz) +
            dy * (Cxx(4, 11) +
                Cxx(4, 6) * dx +
                Cxx(4, 7) * dy +
                Cxx(4, 8) * dz) +
            dz * (Cxx(5, 11) +
                Cxx(5, 6) * dx +
                Cxx(5, 7) * dy +
                Cxx(5, 8) * dz);

        C11 += SigmaXiXj;
        C12 += SigmaXiYj;
        C13 += SigmaXiZj;
        C22 += SigmaYiYj;
        C23 += SigmaYiZj;
        C33 += SigmaZiZj;

        double C11p = (std::pow((trfMat(0, 0) * cosyaw * sinscan + // dX / dRange - measurements
            trfMat(0, 1) * sinyaw * sinscan +
            trfMat(0, 2) * cosscan), 2) * sigmaRange * sigmaRange +
            std::pow((-1 * trfMat(0, 0) * rrange * sinyaw * sinscan + // dX / dYaw
                trfMat(0, 1) * rrange * cosyaw * sinscan), 2) * sigmaYaw * sigmaYaw +
            std::pow((trfMat(0, 0) * rrange * cosyaw * cosscan + // dX / dScan
                trfMat(0, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(0, 2) * rrange * sinscan), 2) * sigmaScan * sigmaScan);

        double C12p = ((trfMat(1, 0) * cosyaw * sinscan + // dY / dRange - measurements
            trfMat(1, 1) * sinyaw * sinscan +
            trfMat(1, 2) * cosscan) *
            (trfMat(0, 0) * cosyaw * sinscan + // dX / dRange - measurements
                trfMat(0, 1) * sinyaw * sinscan +
                trfMat(0, 2) * cosscan) * pow(sigmaRange, 2) +
            (-1 * trfMat(1, 0) * rrange * sinyaw * sinscan + // dY / dYaw
                trfMat(1, 1) * rrange * cosyaw * sinscan) *
            (-1 * trfMat(0, 0) * rrange * sinyaw * sinscan + // dX / dYaw
                trfMat(0, 1) * rrange * cosyaw * sinscan) * pow(sigmaYaw, 2) +
            (trfMat(0, 0) * rrange * cosyaw * cosscan + // dX / dScan
                trfMat(0, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(0, 2) * rrange * sinscan) *
            (trfMat(1, 0) * rrange * cosyaw * cosscan + // dY / dScan
                trfMat(1, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(1, 2) * rrange * sinscan) * pow(sigmaScan, 2));

        double C22p = (pow(trfMat(1, 0) * cosyaw * sinscan + // dY / dRange - measurements
            trfMat(1, 1) * sinyaw * sinscan +
            trfMat(1, 2) * cosscan, 2) * pow(sigmaRange, 2) +
            pow(-1 * trfMat(1, 0) * rrange * sinyaw * sinscan + // dY / dYaw
                trfMat(1, 1) * rrange * cosyaw * sinscan, 2) * pow(sigmaYaw, 2) +
            pow(trfMat(1, 0) * rrange * cosyaw * cosscan + // dY / dScan
                trfMat(1, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(1, 2) * rrange * sinscan, 2) * pow(sigmaScan, 2));

        double C23p = ((trfMat(1, 0) * cosyaw * sinscan + // dY / dRange - measurements
            trfMat(1, 1) * sinyaw * sinscan +
            trfMat(1, 2) * cosscan) *
            (trfMat(2, 0) * cosyaw * sinscan + // dZ / dRange - measurements
                trfMat(2, 1)* sinyaw * sinscan +
                trfMat(2, 2)* cosscan) * pow(sigmaRange, 2) +
            (-1 * trfMat(1, 0) * rrange * sinyaw * sinscan + // dY / dYaw
                trfMat(1, 1)* rrange * cosyaw * sinscan) *
            (-1 * trfMat(2, 0)* rrange * sinyaw * sinscan + // dZ / dYaw
                trfMat(2, 1)* rrange * cosyaw * sinscan) * pow(sigmaYaw, 2) +
            (trfMat(2, 0) * rrange * cosyaw * cosscan + // dZ / dScan
                trfMat(2, 1)* rrange * sinyaw * cosscan +
                -1 * trfMat(2, 2)* rrange * sinscan) *
            (trfMat(1, 0)* rrange * cosyaw * cosscan + // dY / dScan
                trfMat(1, 1)* rrange * sinyaw * cosscan +
                -1 * trfMat(1, 2)* rrange * sinscan) * pow(sigmaScan, 2));

        double C33p = (pow(trfMat(2, 0) * cosyaw * sinscan + // dZ / dRange - measurements
            trfMat(2, 1) * sinyaw * sinscan +
            trfMat(2, 2) * cosscan, 2) * pow(sigmaRange, 2) +
            pow(-1 * trfMat(2, 0) * rrange * sinyaw * sinscan + // dZ / dYaw
                trfMat(2, 1) * rrange * cosyaw * sinscan, 2) * pow(sigmaYaw, 2) +
            pow(trfMat(2, 0) * rrange * cosyaw * cosscan + // dZ / dScan
                trfMat(2, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(2, 2) * rrange * sinscan, 2) * pow(sigmaScan, 2));

        double C13p = ((trfMat(2, 0) * cosyaw * sinscan + // dZ / dRange - measurements
            trfMat(2, 1) * sinyaw * sinscan +
            trfMat(2, 2) * cosscan) *
            (trfMat(0, 0) * cosyaw * sinscan + // dX / dRange - measurements
                trfMat(0, 1) * sinyaw * sinscan +
                trfMat(0, 2) * cosscan) * pow(sigmaRange, 2) +
            (-1 * trfMat(2, 0) * rrange * sinyaw * sinscan + // dZ / dYaw
                trfMat(2, 1) * rrange * cosyaw * sinscan) *
            (-1 * trfMat(0, 0) * rrange * sinyaw * sinscan + // dX / dYaw
                trfMat(0, 1) * rrange * cosyaw * sinscan) * pow(sigmaYaw, 2) +
            (trfMat(2, 0) * rrange * cosyaw * cosscan + // dZ / dScan
                trfMat(2, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(2, 2) * rrange * sinscan) *
            (trfMat(0, 0) * rrange * cosyaw * cosscan + // dX / dScan
                trfMat(0, 1) * rrange * sinyaw * cosscan +
                -1 * trfMat(0, 2) * rrange * sinscan) * pow(sigmaScan, 2));

        C11 += C11p;
        C12 += C12p;
        C22 += C22p;
        C13 += C13p;
        C23 += C23p;
        C33 += C33p;
        Eigen::MatrixXd CxxTmp;
        CxxTmp.resize(3, 3);
        CxxTmp << C11p, C12p, C13p,
            C12p, C22p, C23p,
            C13p, C23p, C33p;
        if (CxxTmp.determinant() == 0)
        {
            CxxTmp = Matrix<double, 3, 3>::Identity();
        }
            
        //CxxTmp = CxxTmp.inverse();
        CxxTmp = CxxTmp.completeOrthogonalDecomposition().pseudoInverse();

        ATP.block<3, 3>(0, i*3) = CxxTmp; // give CxxTmp to ATP

    }

    local_Cxx << C11, C12, C13,
        C12, C22, C23,
        C13, C23, C33;

    Eigen::Matrix3d N = ATP * A;

    Eigen::Matrix3d Qxx = N.completeOrthogonalDecomposition().pseudoInverse(); // 3*3 inverse matrix

    

    //cout << "Qxx is " << Qxx << endl;

    Eigen::Vector3d candidate_mean = std::accumulate(candidate.begin(),
        candidate.end(), Eigen::Vector3d(0, 0, 0)) / candidate.size();

    Eigen::MatrixXd l;
    l.resize(3 * npt, 1);

    for (auto i = 0; i<candidate.size(); i++)
    {
        l.block<3,1>(i * 3, 0) = candidate[i] - candidate_mean;
    }

    Eigen::Matrix3d local_mean = (Qxx * (ATP * l)).transpose() ;

    mean[0] = local_mean(0, 0) + candidate_mean[0];
    mean[1] = local_mean(0, 1) + candidate_mean[1];
    mean[2] = local_mean(0, 2) + candidate_mean[2];

    covariance = local_Cxx / npt;
     
    
}