#pragma once


#include <vector>
#include <memory>
#include <pdal/Options.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/SpatialReference.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <Eigen/Dense>
#include <math.h>



using namespace pdal;
using namespace pdal::Dimension;
using namespace std;
using namespace Eigen;

struct boundingbox {
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
};

vector<Eigen::Vector3d> ReadLas(const char* inputfile, vector<Eigen::Vector3d>& bbox);

vector<Eigen::Vector3d> ReadLasPtSrc(const char* inputfile, vector<Eigen::Vector3d>& bbox, vector<int>& srcid);

void SaveLas(const char* outputfile, vector<Eigen::Vector3d> point3D);

bool SaveLasExtraDims(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<double> val1, vector<double> val2, vector<double> val3);

bool is_file_exist(const char* fileName);
