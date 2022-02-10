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
#include <pdal/io/LasVLR.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/Filter.hpp>
#include <pdal/filters/StreamCallbackFilter.hpp>
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

vector<Eigen::Vector3d> ReadLas(const char* inputfile, vector<Eigen::Vector3d>& bbox,
	vector<Eigen::VectorXd>& MetaData, vector<ExtraDim>& extraDims);

void SaveLas(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<Eigen::VectorXd> MetaData, vector<ExtraDim> extraDims);

bool SaveLasExtraDims(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<Eigen::VectorXd> MetaData, vector<ExtraDim> extraDims,
	vector<double> val1, vector<double> val2, vector<double> val3);

static bool ReadExtraBytesVlr(LasHeader& header, std::vector<ExtraDim>& extraDims);

bool is_file_exist(const char* fileName);

const char* CombineFileName(const char* inputfile, const char* extendName);