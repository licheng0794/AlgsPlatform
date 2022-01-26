#pragma once

# include "../shared/FileLAS.hpp"
# include"../shared/PointCloudAI.hpp"

using namespace std;

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


