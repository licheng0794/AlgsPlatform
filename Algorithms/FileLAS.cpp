# include "FileLAS.hpp"


// @ author: Cheng Li
vector<Eigen::Vector3d> ReadLas(const char* inputfile, vector<Eigen::Vector3d>& bbox)
{
	pdal::Option las_opt("filename", inputfile);
	pdal::Options las_opts;
	las_opts.add(las_opt);

	pdal::PointTable table;
	pdal::LasReader las_reader;
	las_reader.setOptions(las_opts);
	las_reader.prepare(table);

	pdal::PointViewSet point_view_set = las_reader.execute(table);
	pdal::PointViewPtr point_view = *point_view_set.begin();
	pdal::Dimension::IdList dims = point_view->dims();
	pdal::LasHeader las_header = las_reader.header();

	unsigned int PointCount = las_header.pointCount();
	
	double scale_x = las_header.scaleX();
	double scale_y = las_header.scaleY();
	double scale_z = las_header.scaleZ();
	double offset_x = las_header.offsetX();
	double offset_y = las_header.offsetY();
	double offset_z = las_header.offsetZ();

	bbox.at(0) = Eigen::Vector3d(las_header.minX(), las_header.minY(), las_header.minZ());
	bbox.at(1) = Eigen::Vector3d(las_header.maxX(), las_header.maxY(), las_header.maxZ());

	vector<Eigen::Vector3d> point3D(PointCount);

	unsigned int i = 0;
	for (pdal::PointId idx = 0, i=0; (idx < point_view->size()) & (i< PointCount); ++idx, ++i)
	{

		double x = point_view->getFieldAs<double>(Id::X, idx);
		double y = point_view->getFieldAs<double>(Id::Y, idx);
		double z = point_view->getFieldAs<double>(Id::Z, idx);
		
		point3D[i][0] = x;
		point3D[i][1] = y;
		point3D[i][2] = z;
	}

	//cout << "offset_x = " << offset_x;
	//cout << "offset_y = " << offset_y;
	//cout << "offset_z = " << offset_z << endl;
	return point3D; // should set different return values
}

void SaveLas(const char* outputfile, vector<Eigen::Vector3d> point3D)
{
	Options las_opt;
	las_opt.add("filename", outputfile);
	//las_opt.add("extra_dims", "all");
	double xoffset = 0.0, yoffset = 0.0, zoffset = 0.0;

	las_opt.add("offset_x", xoffset);
	las_opt.add("offset_y", yoffset);
	las_opt.add("offset_z", zoffset);
	las_opt.add("scale_x", 0.01); // default
	las_opt.add("scale_y", 0.01);
	las_opt.add("scale_z", 0.01);

	PointTable table;
	table.layout()->registerDim(Dimension::Id::X);
	table.layout()->registerDim(Dimension::Id::Y);
	table.layout()->registerDim(Dimension::Id::Z);
	/*table.layout()->registerOrAssignDim("RedA",
		Dimension::Type::Double);
	table.layout()->registerOrAssignDim("RedB",
		Dimension::Type::Double);
	table.layout()->registerOrAssignDim("RedC",
		Dimension::Type::Double);*/

	//Dimension::Id foo = table.layout()->findDim("RedA");

	PointViewPtr view(new PointView(table));
	for (int i = 0; i < point3D.size(); ++i)
	{
		//cout << point3D[i][0] << endl;
		view->setField(pdal::Dimension::Id::X, i, point3D[i][0]);
		view->setField(pdal::Dimension::Id::Y, i, point3D[i][1]);
		view->setField(pdal::Dimension::Id::Z, i, point3D[i][2]);
		//view->setField(foo, i, NAN);

		// view->setField(pdal::Dimension::Id::Intensity, i, pointCloud[i].intensity);
	}

	BufferReader xjBufferReader;
	xjBufferReader.addView(view);

	StageFactory factory;
	std::string r_drivername = factory.inferReaderDriver(outputfile);
	std::string w_drivername = factory.inferWriterDriver(outputfile);
	Stage* writer = factory.createStage("writers.las");

	writer->setInput(xjBufferReader);
	writer->setOptions(las_opt);
	writer->prepare(table);
	writer->execute(table);

	return;
}

bool SaveLasExtraDims(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<double> val1, vector<double> val2, vector<double> val3)
{
	Options las_opt;
	las_opt.add("filename", outputfile);
	las_opt.add("extra_dims", "M3C2dist=double, SigChg=double, LOD=double");
	double xoffset = 0.0, yoffset = 0.0, zoffset = 0.0;
	las_opt.add("offset_x", xoffset);
	las_opt.add("offset_y", yoffset);
	las_opt.add("offset_z", zoffset);
	las_opt.add("scale_x", 0.01); // default
	las_opt.add("scale_y", 0.01);
	las_opt.add("scale_z", 0.01);

	PointTable table;
	table.layout()->registerDim(Dimension::Id::X);
	table.layout()->registerDim(Dimension::Id::Y);
	table.layout()->registerDim(Dimension::Id::Z);
	table.layout()->registerOrAssignDim("M3C2dist",
		Dimension::Type::Double);
	table.layout()->registerOrAssignDim("SigChg", Dimension::Type::Double);
	table.layout()->registerOrAssignDim("LOD", Dimension::Type::Double);

	Dimension::Id dist = table.layout()->findDim("M3C2dist");
	Dimension::Id sgchg = table.layout()->findDim("SigChg");
	Dimension::Id lod = table.layout()->findDim("LOD");

	PointViewPtr view(new PointView(table));
	for (int i = 0; i < point3D.size(); ++i)
	{
		//cout << point3D[i][0] << endl;
		view->setField(pdal::Dimension::Id::X, i, point3D[i][0]);
		view->setField(pdal::Dimension::Id::Y, i, point3D[i][1]);
		view->setField(pdal::Dimension::Id::Z, i, point3D[i][2]);
		view->setField(dist, i, val1[i]);
		view->setField(sgchg, i, val2[i]);
		view->setField(lod, i, val3[i]);

		// view->setField(pdal::Dimension::Id::Intensity, i, pointCloud[i].intensity);
	}

	BufferReader xjBufferReader;
	xjBufferReader.addView(view);

	StageFactory factory;
	std::string r_drivername = factory.inferReaderDriver(outputfile);
	std::string w_drivername = factory.inferWriterDriver(outputfile);
	Stage* writer = factory.createStage("writers.las");

	writer->setInput(xjBufferReader);
	writer->setOptions(las_opt);
	writer->prepare(table);
	writer->execute(table);

	return true;
}


