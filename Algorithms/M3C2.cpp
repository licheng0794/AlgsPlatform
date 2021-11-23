// @ author: Cheng Li

# include "M3C2.hpp"
# include <armadillo>
#include <ppl.h>
#include <math.h>
#include <ctime>
#include <string.h>
#include <assert.h>
#include <fstream>
#include <iostream>

using namespace arma;
using namespace concurrency;

# define MAX_OCTREE_LEVEL 14
# define IF_DEBUG false

using namespace std;

vector <Eigen::Vector3d> pointxx(2, Eigen::Vector3d(0,0,0));
vector <Eigen::Vector3d> bboxx(2, Eigen::Vector3d(0, 0, 0));


struct M3C2Params {
    vector<double> M3C2dists; // distances
    vector<double> M3C2sigchg; // significant change (0,1)
    vector<double> M3C2lod; // level of detection

    PointCloudAI *pCloud1;
    PointCloudAI *pCloud2;

    vector<Eigen::Vector3d> corepointsnomrals;
    vector<Eigen::Vector3d> corepoints;

    int maxneighbours = 100000;
    int minneighbours = 4; // Lague et.al 2013 preferably 30 

    double projDepth = 5.759067;
    double normalScale = 0.394492;
    double projScale = 0.394492;

    double regError = 0.17;

};

static M3C2Params g_M3C2Params;

int getIndex(vector<Eigen::Vector3d> v, Eigen::Vector3d K)
{
    int a = 1;

    auto it = find(v.begin(), v.end(), K);

    assert(it != v.end()); // If the element is not present in the vector

    // calculating the index
    // of K
    int index = it - v.begin();
    return index;
    
}

void GetCoreNormals(vector <Eigen::Vector3d>originalpoints, vector <Eigen::Vector3d>originalnormals, vector <Eigen::Vector3d>corepoints)
{
    vector <Eigen::Vector3d> corepointsnormals;

    for (int i=0; i < corepoints.size(); i++)
    {
        int index = getIndex(originalpoints, corepoints[i]);
        corepointsnormals.push_back(originalnormals[index]);
    }
    g_M3C2Params.corepointsnomrals = corepointsnormals;
}


void M3C2(const char* incloud1, const char* incloud2, const char* cpcloud, const char* outcloud)
{

    vector<Eigen::Vector3d> bbox(2);
    cout << "1. Reading the reference point cloud!" << endl;
    vector<Eigen::Vector3d> point3D = ReadLas(incloud1, bbox);
    g_M3C2Params.pCloud1 = new PointCloudAI(point3D, bbox);
    cout << "   The reference point cloud has " << point3D.size() << " points!" << endl;

    cout << "2. Reading the target point cloud!" << endl;
    point3D = ReadLas(incloud2, bbox);
    g_M3C2Params.pCloud2 = new PointCloudAI(point3D, bbox);
    cout << "   The target point cloud has " << point3D.size() << " points!" << endl;

    using namespace open3d::geometry;
    PointCloud pointcloud1;
    pointcloud1.points_ = g_M3C2Params.pCloud1->m_point3D;
    cout << "3. Building the kdtree of the reference point cloud!" << endl;
    g_M3C2Params.pCloud1->m_kdtree.SetGeometry(pointcloud1);
    PointCloud pointcloud2;
    pointcloud2.points_ = g_M3C2Params.pCloud2->m_point3D;

    cout << "4. Building the kdtree of the target point cloud!" << endl;
    g_M3C2Params.pCloud2->m_kdtree.SetGeometry(pointcloud2);

    // parameter estimation, return projectionscale, normalscale, projDepth
    double normalScale = 0;
    double projScale = 0;
    double projDepth = 0;

    cout << "5. Estimating proper parameters ..." << endl;
    //ParamEstimate(*g_M3C2Params.pCloud1, *g_M3C2Params.pCloud2, normalScale, projScale, projDepth);

    //g_M3C2Params.corepoints = g_M3C2Params.pCloud1->m_point3D; // temporary
    /*g_M3C2Params.normalScale = normalScale;
    g_M3C2Params.projScale = projScale;
    g_M3C2Params.projDepth = projDepth;*/

    g_M3C2Params.normalScale = 0.315625;
    g_M3C2Params.projScale = 0.315625;
    g_M3C2Params.projDepth = 2.45;
    
    cout << "6. Reading the core points!" << endl;
    if (strcmp(incloud1, cpcloud)==0)
    {
        g_M3C2Params.corepoints = g_M3C2Params.pCloud1->m_point3D;
        cout << "the core point is the same with the reference point cloud!" << endl;
    }
    else
    {
        point3D = ReadLas(cpcloud, bbox);
        g_M3C2Params.corepoints = point3D;
    }
    //g_M3C2Params.corepoints.clear();
    //for (int i = 0; i < 1000; i++)
    //{
    //    g_M3C2Params.corepoints.push_back(g_M3C2Params.pCloud1->m_point3D[i]);
    //}

    int cplen = g_M3C2Params.corepoints.size();
    cout << "   The core point cloud has " << cplen << " points!" << endl;


    cout << "7. Computing normal of core points!" << endl;
    // here we get normals for the whole point cloud !!!
    vector<Eigen::Vector3d> normals = g_M3C2Params.pCloud1->NormalEstimationRadius(g_M3C2Params.normalScale/2); 
    g_M3C2Params.corepointsnomrals = normals; // temporary. normal size should be equal to the corepoint size. if corepoints are changed, 
                                               // normal should be changed also.
    
    /*vector<Eigen::Vector3d> newnormals(normals.size());
    std::ifstream myfile("D:\\work\\AILib\\SDK\\Algorithms\\Algorithms\\M3C2_small_small_cc.asc");
    string line;
    int _size = 9;
    int n = 0;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            char* _sz;
            char* _string = &line[0];

            for (int i = 0; i < _size; i++)
            {
                
                if (i == 6 || i==7 || i==8)
                {
                    newnormals[n][i-6] = strtod(_string, &_sz);
                }
                else
                {
                    strtod(_string, &_sz);
                }

                strcpy(_string, _sz);
            }

            n++;
        }
        myfile.close();
    }*/

    // we should return the normals of core points
    // if core points are part of the original point cloud, we call the function below to get the normals of core points
    //GetCoreNormals(g_M3C2Params.pCloud1->m_point3D, normals, g_M3C2Params.corepoints);
    g_M3C2Params.corepointsnomrals.clear();
    g_M3C2Params.corepointsnomrals = newnormals;

    vector <int> corepointindex(cplen, 0);
    for (int i = 0; i < cplen; i++)
    {
        corepointindex[i] = i;
    }

    // initialize the outputs
    
    g_M3C2Params.M3C2dists = vector<double> (cplen, NAN);
    g_M3C2Params.M3C2sigchg = vector<double>(cplen, NAN);
    g_M3C2Params.M3C2lod = vector<double>(cplen, NAN);
    
    cout << "8. Computing M3C2 distances ... wait..." << endl;
    cout << "   the estimated waiting time is about " << cplen/10000 * 25  << "s (not accurate)" << endl;

    if (IF_DEBUG)
    {
        std::clock_t start;
        double duration;
        start = std::clock(); // get current time
        // 10000 core points take 18.887s
        // 20000 core points take 48.211s
        for (int index = 0; index < corepointindex.size(); index++)
        {
           
            computeM3C2(index);  
            // duration 0.001s for each
        }
        duration = (std::clock()-start) / (double)CLOCKS_PER_SEC;
        std::cout << "Operation took " << duration << "seconds" << std::endl;
        
    }
    else
    {
        std::clock_t start;
        double duration;
        start = std::clock(); // get current time
        // 10000 core points take 4.658s
        // 20000 core points takes 12.874s
        parallel_for_each(begin(corepointindex), end(corepointindex), [](int index) {
            computeM3C2(index);
            });

        duration = (std::clock()-start) / (double)CLOCKS_PER_SEC;
        std::cout << "Operation took " << duration << "seconds" << std::endl;
    }

    // save results to a las file;
    if (SaveLasExtraDims(outcloud, g_M3C2Params.corepoints, g_M3C2Params.M3C2dists, g_M3C2Params.M3C2sigchg, g_M3C2Params.M3C2lod))
    {
        cout << "Saving running results successfully!" << endl;
    }


}

void computeM3C2(int index)
{
    // get normal
    //cout << "the current core point index is " << index << endl;
    Eigen::Vector3d N = g_M3C2Params.corepointsnomrals[index];
    Eigen::Vector3d corepoint = g_M3C2Params.corepoints[index];

    std::vector<int> neigh_candidate_indices(g_M3C2Params.maxneighbours);
    std::vector<double> neigh_candidate_dists(g_M3C2Params.maxneighbours);

    double EffectiveSearchRadius;
    EffectiveSearchRadius = std::sqrt(g_M3C2Params.projDepth * g_M3C2Params.projDepth + g_M3C2Params.projScale / 2 * g_M3C2Params.projScale / 2);
    
    // for the refence point cloud pCloud1
    // get neighbour candidates
    int k = g_M3C2Params.pCloud1->NeighbourSearchRadius(g_M3C2Params.corepoints[index], EffectiveSearchRadius,
        neigh_candidate_indices, neigh_candidate_dists);

 
    //cout << "the current core point have neighbours " << L << endl;
    //cout << "the effectivesearchradius is " << EffectiveSearchRadius << endl;
    //cout << "the current core point in cloud1 have neighbours " << k << endl;

    if (k < g_M3C2Params.minneighbours)
    {
        return; // not able to compute statistcs
    }
    
    vector<Eigen::Vector3d> candidate;
    for (int i = 0; i < k; i++)
    {
        candidate.push_back(g_M3C2Params.pCloud1->m_point3D[neigh_candidate_indices[i]]);
    }

    vector<double> projDist1; // the distances of points projected to the normal direction
    int n1 = GetRealNeighbours(candidate, corepoint, N, g_M3C2Params.projScale / 2, g_M3C2Params.projDepth, projDist1);
    
    if (n1 < g_M3C2Params.minneighbours)
    {
        return; // not able to compute statistcs
    }

    // for the targeted point cloud pCloud2
    // get neighbour candidates
    k = g_M3C2Params.pCloud2->NeighbourSearchRadius(g_M3C2Params.corepoints[index], EffectiveSearchRadius,
        neigh_candidate_indices, neigh_candidate_dists);
    //cout << "the current core point in cloud2 have neighbours " << k << endl;
    if (k < g_M3C2Params.minneighbours)
    {
        return; // not able to compute statistcs
    }

    candidate.clear();
    for (int i = 0; i < k; i++)
    {
        candidate.push_back(g_M3C2Params.pCloud2->m_point3D[neigh_candidate_indices[i]]);
    }

    vector<double> projDist2; // the distances of points projected to the normal direction

    int n2 = GetRealNeighbours(candidate, corepoint, N, g_M3C2Params.projScale / 2, g_M3C2Params.projDepth, projDist2);

    if (n2 < g_M3C2Params.minneighbours)
    {
        return; // not able to compute statistcs
    }

    // using projDist1 and projDist2 to compute M3C2 distance, level of detection, significant change
    
    double mean1 = std::accumulate(projDist1.begin(), projDist1.end(), 0.0) / n1;
    double mean2 = std::accumulate(projDist2.begin(), projDist2.end(), 0.0) / n2;

    double stddev1 = stdvariance(projDist1);
    double stddev2 = stdvariance(projDist2);

    double dist = mean2 - mean1;
    double stddev = (stddev1 * stddev1) / n1 + (stddev2 * stddev2) / n2;
    double lod = 1.96 * (sqrt(stddev) + g_M3C2Params.regError);
    
    g_M3C2Params.M3C2dists[index] = dist;
    g_M3C2Params.M3C2lod[index] = lod;

    if (dist<-lod || dist > lod) // significant change
    {
        g_M3C2Params.M3C2sigchg[index] = 1;
    }


}

// computing dotproduct between two vectors
double dotproduct(Eigen::Vector3d v1, Eigen::Vector3d v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// computing the distance of a vector to the normal
void computeAlongAcross(Eigen::Vector3d invector, Eigen::Vector3d N, double &along, double &across)
{
    along = dotproduct(invector, N);
    Eigen::Vector3d alongvector = invector - N * along;
    across = std::pow(alongvector.norm(), 2);
}

int GetRealNeighbours(vector <Eigen::Vector3d>candidates, Eigen::Vector3d corepoint, Eigen::Vector3d N, 
                                           double radius, double projDepth, vector<double> &dist)
{

    vector <Eigen::Vector3d> RealNeighbours;

    vector<double> along(candidates.size(), 0);
    vector<double> across(candidates.size(), 0);
    vector<int> candindex(candidates.size(), 0);
    for (int i = 0; i < candidates.size(); i++)
    {
        candindex[i] = i;
    }

    double alongval, acrossval;
    if (IF_DEBUG)
    {
        for (int i = 0; i < candidates.size(); i++)
        {
            
            computeAlongAcross(candidates[i] - corepoint, N, alongval, acrossval);
            along[i] = alongval;
            across[i] = acrossval;
        }
    }
    else
    {
     /*   parallel_for_each(candindex.begin(), candindex.end(), [&](int ind) {
            computeAlongAcross(candidates[ind] - corepoint, N, alongval, acrossval);
            along[ind] = alongval;
            across[ind] = acrossval;
            });*/

        for (int i = 0; i < candidates.size(); i++)
        {
            computeAlongAcross(candidates[i] - corepoint, N, alongval, acrossval);
            along[i] = alongval;
            across[i] = acrossval;
        }
    }

    // make sure along <= projDepth and across <= radius*radius
    for (int i = 0; i < candidates.size(); i++)
    {
        if (std::abs(along[i]) <= projDepth && across[i] <= radius * radius)
        {
            RealNeighbours.push_back(candidates[i]);
            dist.push_back(along[i]);
        }
    }

    return RealNeighbours.size();
 
}



template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
    
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        srand(time(NULL));
        std::advance(r, rand() % left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

bool planefit(vector<Eigen::Vector3d> points, double n[4])
{
    // n[0,1,2]      An array of three doubles containing the unit normal of the best fit plane.
    // n[3] 

    // The Aramadillo library was used in this work:
    // Conrad Sanderson and Ryan Curtin
    // "Armadillo: a template - based C++ library for linear algebra".
    //	Journal of Open Source Software, Vol. 1, pp. 26, 2016.
    // See http://arma.sourceforge.net/ for details.
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 1.0;

    Eigen::Vector3d centriod(0, 0, 0);

    for (size_t i = 0; i < points.size(); ++i) {
        centriod += points[i];
    }

    centriod = centriod / (double)points.size();

    // Instantiate a 3 x N armadillo matrix X, where N is the number of input points.
    // The (arma::uword) cuts out warnings from the compiler concerning
    // potential loass of data from conversion.

    mat X((arma::uword)3, (arma::uword)(points.size()));

    // Initialize X with the co-ordinates of the input points relative to the centroid.
    // This, in effect, 'projects' the intput point set to the origin.

    for (int i = 0; i < (int)points.size(); ++i) {
        X(0, i) = points[i][0] - centriod[0];
        X(1, i) = points[i][1] - centriod[1];
        X(2, i) = points[i][2] - centriod[2];
    }

    // The next step is to perform a singular value decomposition of the 3 x N matrix X.
    // This means finding a diagonal 3 x N matrix S such that X = USV^t, where U is 3x3 orthogonal and
    // V is NxN orthogonal.
    // Because the matrix S is diagonal, we only need a vector to house its values.

    mat U;  // Instantiate an armadillo matrix for U.
    vec S;  // Instantiate an armadillo vector for S.
    mat V;  // Instantiate an armadillo matrix for V.

    // Perform a Singular Value Decomposition of the matrix X using armadillo.

    if (!svd(U, S, V, X, "std")) {
        return false;
    }
    // U is now a 3x3 orthogonal matrix
    // S is a vector of dimension 3.
    // V is now an N x N orthogonal matrix.

    // The normal to the plane is the eigenvector in U that corresponds to the
    // minimal value (eigenvector) in S. Note the () overload used by armadillo
    // on the vector S.

    int imin = 0;
    if (abs(S(1)) < abs(S(0))) imin = 1;
    if (abs(S(2)) < abs(S(imin))) imin = 2;

    n[0] = U(0, imin);
    n[1] = U(1, imin);
    n[2] = U(2, imin);

    n[3] = -1 * (n[0] * centriod[0] + n[1] * centriod[1] + n[2] * centriod[2]);

    return true;

}

template<typename T>
T stdvariance(const std::vector<T>& vec) {
    const size_t sz = vec.size();
    if (sz == 1) {
        return 0.0;
    }

    // Calculate the mean
    const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    // Now calculate the variance
    auto variance_func = [&mean, &sz](T accumulator, const T& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz));
    };

    return std::sqrt(std::accumulate(vec.begin(), vec.end(), 0.0, variance_func));
}

static int m_cellnum[MAX_OCTREE_LEVEL] = { 0 };

bool f_traverse(const std::shared_ptr<open3d::geometry::OctreeNode>& node,
    const std::shared_ptr<open3d::geometry::OctreeNodeInfo>& node_info) {
    if (auto internal_node =
        std::dynamic_pointer_cast<open3d::geometry::OctreeInternalNode>(node)) {
        if (auto internal_point_node = std::dynamic_pointer_cast<
            open3d::geometry::OctreeInternalPointNode>(internal_node)) {
            int num_children = 0;
            for (const auto& c : internal_point_node->children_) {
                if (c) num_children++;
            }
            /*open3d::utility::LogInfo(
                "Internal node at depth {} with origin {} has "
                "{} children and {} points",
                node_info->depth_, node_info->origin_, num_children,
                internal_point_node->indices_.size());*/
            m_cellnum[node_info->depth_] += num_children;
        }
    }
    else if (auto leaf_node = std::dynamic_pointer_cast<
        open3d::geometry::OctreePointColorLeafNode>(node)) {
        m_cellnum[node_info->depth_] += 1;
        /*open3d::utility::LogInfo(
            "Node at depth {} with origin {} has"
            "color {} and {} points",
            node_info->depth_, node_info->origin_, leaf_node->color_,
            leaf_node->indices_.size());*/
        // utility::LogInfo("Indices: {}", leaf_node->indices_);
    }
    else {
        open3d::utility::LogError("Unknown node type");
    }

    return false;
}


void ParamEstimate(PointCloudAI &pCloud1, PointCloudAI &pCloud2, double& bestNormalScale, double& bestProjScale, double& bestProjDepth)
{
    // build octree and get its size
    vector <double> cellsize(MAX_OCTREE_LEVEL, 0);
    auto octree = std::make_shared<open3d::geometry::Octree>(MAX_OCTREE_LEVEL);

    using namespace open3d::geometry;
    PointCloud pointcloud;
    pointcloud.points_ = pCloud1.m_point3D;
    //pCloud1.m_kdtree.SetGeometry(pointcloud); // necessary steps for neighbour search. can be considered to move out the function

    //pointcloud.points_ = pCloud2.m_point3D;
    //pCloud2.m_kdtree.SetGeometry(pointcloud);

    octree->ConvertFromPointCloud(pointcloud);

    octree->Traverse(f_traverse);

    cout << "   the octree size is " << octree->size_ << endl;

    double octreesize = octree->size_;
    cellsize[0] = octreesize;
    double startsearchscale = 0;
    for (int i = 1; i < MAX_OCTREE_LEVEL; i++)
    {
        cellsize[i] = cellsize[i-1] / 2;
      /*  if (cellsize[i] > 0.01 && cellsize[i-1] < 0.01)
            startsearchscale = cellsize[i];*/
    }

    for (int i = 1; i < MAX_OCTREE_LEVEL; i++)
    {
        if (pCloud1.m_point3D.size() / m_cellnum[i] < 1.5)
        {
            startsearchscale = cellsize[i-1];
            break;
        }
    }

    cout << "   the searching basescale is " << startsearchscale << endl;

    // let us start to search
    
    int NoPoints = 10000; // the number of random points to verify the selection criteria
    vector<Eigen::Vector3d> selectedPoints = pCloud1.m_point3D;
    // random points have been saved to the first 'NoPoints' in selectedPoints
    random_unique(selectedPoints.begin(), selectedPoints.end(), NoPoints);
    
    bool findNormalScale = false;
    double normalscale = startsearchscale;
    
    while (!findNormalScale)
    {
        
        int planeNum = 0;
        double sumRoughness = 0;
        
        // compute roughness and find the normalscale corresponding to the lowest roughness & \sigma(D)/D>=25 (Lague et.al 2013)
        int i = 0;
        
        vector<double> vecRoughness;

        for (i = 0; i < NoPoints; i++)
        {
            /*cout << "the current point is " << i << endl;*/
            // compute neighbours
            int nn = 1000;
            double roughness = 0;
            vector<int> new_indices_vec(nn);
            vector<double> new_dists_vec(nn);
            int k = pCloud1.NeighbourSearchRadius(selectedPoints[i], normalscale / 2, new_indices_vec, new_dists_vec);
            if (k >= 10) // Lange et al. 2013
            {
                planeNum++;
                // fit plane
                double N[4];

                // store all neighbours
                vector <Eigen::Vector3d> neighbours;
                for (int j = 0; j < k; j++)
                {
                    neighbours.push_back(pCloud1.m_point3D[new_indices_vec[j]]);
                }

                planefit(neighbours, N); // 
                
                Eigen::Vector3d normals(N[0], N[1], N[2]);

                vector <double> distances;

                for (int j = 0; j < k; j++)
                {
                    distances.push_back(std::abs(pCloud1.m_point3D[new_indices_vec[j]].dot(normals) + N[3]));
                }

                roughness = stdvariance(distances);

                /*double distances;
                distances = std::abs(selectedPoints[i].dot(normals) + N[3]);
                roughness = distances;*/

                sumRoughness += roughness/normalscale; // relative roughness. larger the normalscale, larger the absolutive roughness 
                vecRoughness.push_back(roughness);
            }
            
        }

        double meanRoughness = 1000;
        cout << "   the current normalscale is " << normalscale  << " roughness is "  << sumRoughness / planeNum  << endl;
        if (sumRoughness / planeNum <= meanRoughness)
            meanRoughness = sumRoughness / planeNum;
        else
        {
            if (normalscale/ stdvariance(vecRoughness) >=25)
            {
                bestNormalScale = normalscale;
                findNormalScale = true;
            }
        }
        

        if (normalscale > cellsize[7])
        {
            cout << "   the searching size is over level 7!" << endl;
            findNormalScale = true;
            bestNormalScale = cellsize[7];
        }
        normalscale += startsearchscale;
    }

    cout << "   the best normalscale is " << bestNormalScale << endl;

    // working for projection scale
    bool findProjScale = false;
    double projscale = startsearchscale;

   
    vector<Eigen::Vector3d> selectedPoints2 = pCloud2.m_point3D;
    // random points have been saved to the first 'NoPoints' in selectedPoints
    random_unique(selectedPoints2.begin(), selectedPoints2.end(), NoPoints);

    while (!findProjScale)
    {
        // make sure the number of neighbours in reference and target clouds should be as least 30 (Lague et.al 2013)
        vector <int> pop1;
        vector <int> pop2;
        int nn = 1000;
        vector<int> new_indices_vec(nn);
        vector<double> new_dists_vec(nn);

        for (int i = 0; i < NoPoints; i++)
        {
            int k = pCloud1.NeighbourSearchRadius(selectedPoints[i], projscale / 2, new_indices_vec, new_dists_vec);
            pop1.push_back(k);

            k = pCloud2.NeighbourSearchRadius(selectedPoints2[i], projscale / 2, new_indices_vec, new_dists_vec);
            pop2.push_back(k);
        }

        std::sort(pop1.begin(), pop1.end());
        std::sort(pop2.begin(), pop2.end());

        if (pop1[(int)(NoPoints * 0.03)] >= 30 && pop2[(int)(NoPoints * 0.03)] >= 30)
        {
            findProjScale = true;
            bestProjScale = projscale;
        }
        projscale += startsearchscale;

    }

    cout << "   the best projscale is " << bestProjScale << endl;

    // working for projDepth
    Eigen::Vector3d bbox1 = pCloud1.m_maxbound - pCloud1.m_minbound;
    Eigen::Vector3d bbox2 = pCloud2.m_maxbound - pCloud2.m_minbound;

    bestProjDepth = std::min(bbox1.norm(), bbox2.norm()) * 0.05; // guess

    cout << "   the best projdepth is " << bestProjDepth << endl;


}