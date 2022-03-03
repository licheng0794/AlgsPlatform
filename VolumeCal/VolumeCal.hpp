#pragma once

# include "../shared/FileLAS.hpp"
# include "../shared/PointCloudAI.hpp"
# include "open3d/io/TriangleMeshIO.h"

# include <armadillo>

using namespace std;
using namespace open3d::geometry;
using namespace arma;

// the cell info
struct Cell {
	int Npoints; // the number of points in the cell
	bool Binterpolate; // interpolate or not ?
	double Height;
	double avgHeight; // the average height
	double stdHeight; // the standard dev of heights
	double minHeight; // the minimum height
	double maxHeight; // the maximum height
	unsigned int minPointIndex; // the point index with the lowest height
	unsigned int maxPointIndex; // the point index with the highest height
};

struct VolumeResults
{
	double VolumeDiff = 0;
	double Surface = 0;
	double addedVolume = 0;
	double removedVolume = 0;
	float matchPercentage = 0;
	float nonmatchGroundPercentage = 0;
	float nonmatchCeilPercentage = 0;
	
};


enum class ProjectionType { AVERAGE, MAXIMUM, MINIMUM};
enum class FillEmptyStrategy {EMPTY, INTERPOLATE};
void Volcal(const char* volfile);
void Volcal(const char* groundfile, const char* ceilfile);
void Volcal(const char* groundfile, const char* ceilfile, double gridstep);
bool VolumeDiffCal(vector<Eigen::Vector3d> groundpoints, vector<Eigen::Vector3d> ceilpoints, 
	Eigen::Vector3d miniCorner, Eigen::Vector3d maxiCorner,
	double gridstep, VolumeResults& result, double &ratio, 
	FillEmptyStrategy fillstrategy = FillEmptyStrategy::EMPTY, bool fixedgrid=false, int projtype=0);

void VolPointCloud(const char* volfile);

void SaveResult(const char* SavedReportName, double gridstep, VolumeResults result, bool optimalgrid=false);

void PaintMesh(open3d::geometry::TriangleMesh& mesh,
	const vector<Eigen::Vector3d>& color);

template<typename T>
T stdvariance(const std::vector<T>& vec);

bool planefitABC(vector<Eigen::Vector3d> points, double n[4]);

void ProjDisPoint(vector<Eigen::Vector3d> points, double* n, 
	vector<Eigen::Vector3d>& Projpoints, vector<double>& distoplane);

double trianglearea(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);
