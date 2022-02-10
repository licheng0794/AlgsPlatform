
# include "FileLAS.hpp"


static bool ReadExtraBytesVlr(LasHeader& header, std::vector<ExtraDim>& extraDims)
{
	const LasVLR* vlr = header.findVlr(SPEC_USER_ID, EXTRA_BYTES_RECORD_ID);
	if (!vlr)
	{
		// cout << "no VLR in the las!" << endl;
		return false;
	}

	size_t size = vlr->dataLen();
	if (size % sizeof(ExtraBytesSpec) != 0)
	{
		// cout << "[LAS] Bad size for extra bytes VLR. Ignoring." << endl;
		return false;
	}
	size_t count = size / sizeof(ExtraBytesSpec);
	//cout << "[LAS] VLR count: " << count << endl;;

	try
	{
		const char* pos = vlr->data();
		for (size_t i = 0; i < count; ++i)
		{
			ExtraBytesIf eb;
			eb.readFrom(pos);
			pos += sizeof(ExtraBytesSpec);

			std::vector<ExtraDim> eds = eb.toExtraDims();
			for (const ExtraDim& ed : eds)
			{
				// cout << "VLR " << (i + 1) << " : " << ed.m_name << endl;
				extraDims.push_back(ed);
			}
		}
	}
	catch (const std::bad_alloc&)
	{
		cout << "[LAS] Not enough memory to retrieve the extra bytes fields." << endl;
		return false;
	}

	return true;
}

void ExtraBytesIf::readFrom(const char* buf)
{
	LeExtractor extractor(buf, sizeof(ExtraBytesSpec));
	uint16_t dummy16;
	uint32_t dummy32;
	uint64_t dummy64;
	double dummyd;
	uint8_t options;
	uint8_t type;

	uint8_t SCALE_MASK = 1 << 3;
	uint8_t OFFSET_MASK = 1 << 4;

	extractor >> dummy16 >> type >> options;
	extractor.get(m_name, 32);
	extractor >> dummy32;
	for (size_t i = 0; i < 3; ++i)
		extractor >> dummy64;  // No data field.
	for (size_t i = 0; i < 3; ++i)
		extractor >> dummyd;  // Min.
	for (size_t i = 0; i < 3; ++i)
		extractor >> dummyd;  // Max.
	for (size_t i = 0; i < 3; ++i)
		extractor >> m_scale[i];
	for (size_t i = 0; i < 3; ++i)
		extractor >> m_offset[i];
	extractor.get(m_description, 32);

	setType(type);
	if (m_type == Dimension::Type::None)
		m_size = options;
	if (!(options & SCALE_MASK))
		for (size_t i = 0; i < 3; ++i)
			m_scale[i] = 1.0;
	if (!(options & OFFSET_MASK))
		for (size_t i = 0; i < 3; ++i)
			m_offset[i] = 0.0;
}

void ExtraBytesIf::setType(uint8_t lastype)
{
	m_fieldCnt = 1;
	while (lastype > 10)
	{
		m_fieldCnt++;
		lastype -= 10;
	}
	using DT = Dimension::Type;
	const Dimension::Type lastypes[] = {
		DT::None, DT::Unsigned8, DT::Signed8, DT::Unsigned16, DT::Signed16,
		DT::Unsigned32, DT::Signed32, DT::Unsigned64, DT::Signed64,
		DT::Float, DT::Double
	};
	m_type = lastypes[lastype];
	if (m_type == Dimension::Type::None)
		m_fieldCnt = 0;
}

std::vector<ExtraDim> ExtraBytesIf::toExtraDims()
{
	std::vector<ExtraDim> eds;

	if (m_type == Dimension::Type::None)
	{
		ExtraDim ed(m_name, Dimension::Type::None);
		ed.m_size = m_size;
		eds.push_back(ed);
	}
	else if (m_fieldCnt == 1)
	{
		ExtraDim ed(m_name, m_type, m_scale[0], m_offset[0]);
		eds.push_back(ed);
	}
	else
	{
		for (size_t i = 0; i < m_fieldCnt; ++i)
		{
			ExtraDim ed(m_name + std::to_string(i), m_type,
				m_scale[i], m_offset[i]);
			eds.push_back(ed);
		}
	}
	return eds;
}


// @ author: Cheng Li
vector<Eigen::Vector3d> ReadLas(const char* inputfile, vector<Eigen::Vector3d>& bbox,
	vector<Eigen::VectorXd>& MetaData, vector<ExtraDim>& extraDims)
{
	pdal::Options las_opts;
	las_opts.add("filename", inputfile);

	LasReader las_reader;
	las_reader.setOptions(las_opts);
	FixedPointTable fields(100);
	PointLayoutPtr layout(fields.layout());
	las_reader.prepare(fields);

	pdal::LasHeader las_header = las_reader.header();
	unsigned int PointCount = las_header.pointCount();
	vector<Eigen::Vector3d> point3D;
	if (PointCount == 0)
	{
		cout << "the input las file is empty!" << endl;
		return point3D;
	}

	std::vector<Id> extraDimensionsIds;
	StringList extraNames;
	string extraDimsArg;

	if (ReadExtraBytesVlr(las_header, extraDims)) // has extra dims in the las
	{
		for (unsigned i = 0; i < extraDims.size(); ++i)
		{
			extraDimsArg += extraDims[i].m_name + "=" + interpretationName(extraDims[i].m_dimType.m_type) + ",";

			extraNames.push_back(extraDims[i].m_name);
		}
	}

	if (!extraNames.empty())
	{
		Options las_opts2;
		las_opts2.add("extra_dims", extraDimsArg);

		las_reader.addOptions(las_opts2);
		las_reader.prepare(fields);

		for (std::string& dim : extraNames)
		{
			extraDimensionsIds.push_back(layout->findDim(dim));
		}

	}

	int nbPointsRead = 0;
	StreamCallbackFilter f;
	f.setInput(las_reader);


	/*
	*   X,
		Y,
		Z,
		Intensity,
		ReturnNumber,
		NumberOfReturns,
		ScanDirectionFlag,
		EdgeOfFlightLine,
		Classification,
		ScanAngleRank,
		UserData,
		PointSourceId,
		Red,
		Green,
		Blue,
		GpsTime,
	*/

	int metalen = 13 + extraDimensionsIds.size();

	auto ccProcessOne = [&](PointRef& point)
	{
		Eigen::Vector3d P(point.getFieldAs<double>(Id::X),
			point.getFieldAs<double>(Id::Y),
			point.getFieldAs<double>(Id::Z));

		Eigen::VectorXd metadata(metalen);

		metadata(0) = point.getFieldAs<unsigned int>(Id::Intensity);
		metadata(1) = point.getFieldAs<unsigned short int>(Id::ReturnNumber);
		metadata(2) = point.getFieldAs<unsigned short int>(Id::NumberOfReturns);
		metadata(3) = point.getFieldAs<unsigned short int>(Id::ScanDirectionFlag);
		metadata(4) = point.getFieldAs<unsigned short int>(Id::EdgeOfFlightLine);
		metadata(5) = point.getFieldAs<unsigned short int>(Id::Classification);
		metadata(6) = point.getFieldAs<float>(Id::ScanAngleRank);
		metadata(7) = point.getFieldAs<unsigned short int>(Id::UserData);
		metadata(8) = point.getFieldAs<unsigned int>(Id::PointSourceId);
		metadata(9) = point.getFieldAs<unsigned int>(Id::Red);
		metadata(10) = point.getFieldAs<unsigned int>(Id::Green);
		metadata(11) = point.getFieldAs<unsigned int>(Id::Blue);
		metadata(12) = point.getFieldAs<double>(Id::GpsTime);

		if (!extraDimensionsIds.empty())
		{
			for (int l = 0; l < extraDimensionsIds.size(); l++)
			{
				metadata(int(13 + l)) = point.getFieldAs<double>(extraDimensionsIds[l]);
			}
		}

		MetaData.push_back(metadata);
		point3D.push_back(P);
		++nbPointsRead;
		return true;
	};

	f.setCallback(ccProcessOne);
	f.prepare(fields);
	f.execute(fields);

	double scale_x = las_header.scaleX();
	double scale_y = las_header.scaleY();
	double scale_z = las_header.scaleZ();
	double offset_x = las_header.offsetX();
	double offset_y = las_header.offsetY();
	double offset_z = las_header.offsetZ();

	bbox.at(0) = Eigen::Vector3d(las_header.minX(), las_header.minY(), las_header.minZ());
	bbox.at(1) = Eigen::Vector3d(las_header.maxX(), las_header.maxY(), las_header.maxZ());

	return point3D;
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

void SaveLas(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<Eigen::VectorXd> MetaData, vector<ExtraDim> extraDims)
{
	Options las_opt;
	las_opt.add("filename", outputfile);

	// if necessary, this values should come from the original file
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
	table.layout()->registerDim(Dimension::Id::Intensity);
	table.layout()->registerDim(Dimension::Id::ReturnNumber);
	table.layout()->registerDim(Dimension::Id::NumberOfReturns);
	table.layout()->registerDim(Dimension::Id::ScanDirectionFlag);
	table.layout()->registerDim(Dimension::Id::EdgeOfFlightLine);
	table.layout()->registerDim(Dimension::Id::Classification);
	table.layout()->registerDim(Dimension::Id::ScanAngleRank);
	table.layout()->registerDim(Dimension::Id::UserData);
	table.layout()->registerDim(Dimension::Id::PointSourceId);
	table.layout()->registerDim(Dimension::Id::Red);
	table.layout()->registerDim(Dimension::Id::Green);
	table.layout()->registerDim(Dimension::Id::Blue);
	table.layout()->registerDim(Dimension::Id::GpsTime);
	std::vector<Id> extraDimensionsIds;
	//StringList extraNames;
	string extraDimsArg;
	if (!extraDims.empty())
	{
		for (unsigned i = 0; i < extraDims.size(); ++i)
		{
			//extraDims[i].m_name
			string name = "Original_cloud_index";
			extraDimsArg += name + "=" + interpretationName(extraDims[i].m_dimType.m_type) + ",";

			//extraNames.push_back(extraDims[i].m_name);

			table.layout()->registerOrAssignDim(name, extraDims[i].m_dimType.m_type);

			extraDimensionsIds.push_back(table.layout()->findDim(name));
		}

		las_opt.add("extra_dims", extraDimsArg);
	}

	PointViewPtr view(new PointView(table));
	for (int i = 0; i < point3D.size(); ++i)
	{
		//cout << point3D[i][0] << endl;
		view->setField(pdal::Dimension::Id::X, i, point3D[i][0]);
		view->setField(pdal::Dimension::Id::Y, i, point3D[i][1]);
		view->setField(pdal::Dimension::Id::Z, i, point3D[i][2]);

		view->setField(pdal::Dimension::Id::Intensity, i, MetaData[i][0]);
		view->setField(pdal::Dimension::Id::ReturnNumber, i, MetaData[i][1]);
		view->setField(pdal::Dimension::Id::NumberOfReturns, i, MetaData[i][2]);
		view->setField(pdal::Dimension::Id::ScanDirectionFlag, i, MetaData[i][3]);
		view->setField(pdal::Dimension::Id::EdgeOfFlightLine, i, MetaData[i][4]);
		view->setField(pdal::Dimension::Id::Classification, i, MetaData[i][5]);
		view->setField(pdal::Dimension::Id::ScanAngleRank, i, MetaData[i][6]);
		view->setField(pdal::Dimension::Id::UserData, i, MetaData[i][7]);
		view->setField(pdal::Dimension::Id::PointSourceId, i, MetaData[i][8]);
		view->setField(pdal::Dimension::Id::Red, i, MetaData[i][9]);
		view->setField(pdal::Dimension::Id::Green, i, MetaData[i][10]);
		view->setField(pdal::Dimension::Id::Blue, i, MetaData[i][11]);
		view->setField(pdal::Dimension::Id::GpsTime, i, MetaData[i][12]);

		if (!extraDims.empty())
		{
			for (int l = 0; l < extraDims.size(); l++)
			{
				view->setField(extraDimensionsIds[l], i, MetaData[i][int(13 + l)]);
			}
		}
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


// only works for M3C2 result saving
bool SaveLasExtraDims(const char* outputfile, vector<Eigen::Vector3d> point3D,
	vector<Eigen::VectorXd> MetaData, vector<ExtraDim> extraDims,
	vector<double> val1, vector<double> val2, vector<double> val3)
{
	Options las_opt;
	las_opt.add("filename", outputfile);

	// if necessary, this values should come from the original file
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
	table.layout()->registerDim(Dimension::Id::Intensity);
	table.layout()->registerDim(Dimension::Id::ReturnNumber);
	table.layout()->registerDim(Dimension::Id::NumberOfReturns);
	table.layout()->registerDim(Dimension::Id::ScanDirectionFlag);
	table.layout()->registerDim(Dimension::Id::EdgeOfFlightLine);
	table.layout()->registerDim(Dimension::Id::Classification);
	table.layout()->registerDim(Dimension::Id::ScanAngleRank);
	table.layout()->registerDim(Dimension::Id::UserData);
	table.layout()->registerDim(Dimension::Id::PointSourceId);
	table.layout()->registerDim(Dimension::Id::Red);
	table.layout()->registerDim(Dimension::Id::Green);
	table.layout()->registerDim(Dimension::Id::Blue);
	table.layout()->registerDim(Dimension::Id::GpsTime);
	std::vector<Id> extraDimensionsIds;
	//StringList extraNames;
	string extraDimsArg;
	if (!extraDims.empty())
	{
		for (unsigned i = 0; i < extraDims.size(); ++i)
		{
			//extraDims[i].m_name
			string name = "Original_cloud_index";
			extraDimsArg += name + "=" + interpretationName(extraDims[i].m_dimType.m_type) + ",";

			//extraNames.push_back(extraDims[i].m_name);

			table.layout()->registerOrAssignDim(name, extraDims[i].m_dimType.m_type);

			extraDimensionsIds.push_back(table.layout()->findDim(name));
		}

		las_opt.add("extra_dims", extraDimsArg);
	}

	las_opt.add("extra_dims", "M3C2dist=double, SigChg=double, LOD=double");
	table.layout()->registerOrAssignDim("M3C2dist",Dimension::Type::Double);
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

		view->setField(pdal::Dimension::Id::Intensity, i, MetaData[i][0]);
		view->setField(pdal::Dimension::Id::ReturnNumber, i, MetaData[i][1]);
		view->setField(pdal::Dimension::Id::NumberOfReturns, i, MetaData[i][2]);
		view->setField(pdal::Dimension::Id::ScanDirectionFlag, i, MetaData[i][3]);
		view->setField(pdal::Dimension::Id::EdgeOfFlightLine, i, MetaData[i][4]);
		view->setField(pdal::Dimension::Id::Classification, i, MetaData[i][5]);
		view->setField(pdal::Dimension::Id::ScanAngleRank, i, MetaData[i][6]);
		view->setField(pdal::Dimension::Id::UserData, i, MetaData[i][7]);
		view->setField(pdal::Dimension::Id::PointSourceId, i, MetaData[i][8]);
		view->setField(pdal::Dimension::Id::Red, i, MetaData[i][9]);
		view->setField(pdal::Dimension::Id::Green, i, MetaData[i][10]);
		view->setField(pdal::Dimension::Id::Blue, i, MetaData[i][11]);
		view->setField(pdal::Dimension::Id::GpsTime, i, MetaData[i][12]);
		view->setField(dist, i, val1[i]);
		view->setField(sgchg, i, val2[i]);
		view->setField(lod, i, val3[i]);
		if (!extraDims.empty())
		{
			for (int l = 0; l < extraDims.size(); l++)
			{
				view->setField(extraDimensionsIds[l], i, MetaData[i][int(13 + l)]);
			}
		}
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

bool is_file_exist(const char* fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

// add extended name to the input file name. 
// e.g the input file is 'abc.las', extendName is '_Volreport.txt'. The returned value
// would be abc_Volreport.txt
const char* CombineFileName(const char* inputfile, const char* extendName)
{
	std::string strname(inputfile);

	std::string base_filename = strname.substr(strname.find_last_of("/\\") + 1);

	std::string::size_type const p(base_filename.find_last_of('.'));

	std::string file_without_extension = base_filename.substr(0, p);

	string extname(extendName);

	std::string groundname = file_without_extension + extname;

	return groundname.c_str();
}
