# include "VolumeCal.hpp"
# include "Delaunator.hpp"
#define NIEGH_RATIO 7.8

// automatical volume calculation for a single point cloud

void SaveResult(const char*SavedReportName, double gridstep, VolumeResults result, bool optimalgrid)
{
	std::ofstream out(SavedReportName);
	std::string resultstr;
	if (optimalgrid)
	{
		resultstr = "the optimal girdstep is " + std::to_string(gridstep);
	}
	else
	{
		resultstr = "the girdstep is " + std::to_string(gridstep);
	}
	 
	out << resultstr;
	out << endl;

	resultstr = "---------------Volume Info---------------"; 
	out << resultstr;
	out << endl;

	resultstr = "the volume difference is " + std::to_string(result.VolumeDiff);
	out << resultstr;
	out << endl;

	resultstr = "the added volume is " + std::to_string(result.addedVolume);
	out << resultstr;
	out << endl;

	resultstr = "the removed volume is " + std::to_string(result.removedVolume);
	out << resultstr;
	out << endl;

	resultstr = "the surface is " + std::to_string(result.Surface);
	out << resultstr;
	out << endl;

	resultstr = "the matching percentage between ground and ceil is " + std::to_string(result.matchPercentage) + "%";
	out << resultstr;
	out << endl;
	
	resultstr = "the non-matching percentage ground is " + std::to_string(result.nonmatchGroundPercentage) + "%";
	out << resultstr;
	out << endl;

	resultstr = "the non-matching percentage ceil is " + std::to_string(result.nonmatchCeilPercentage) + "%";
	out << resultstr;
	out << endl;
}

void Volcal(const char* volfile)
{

	string filestr(volfile);

	if (filestr.substr(filestr.find_last_of(".") + 1) == "las" |
		filestr.substr(filestr.find_last_of(".") + 1) == "laz")
	{
		VolPointCloud(volfile);
	}
	if (filestr.substr(filestr.find_last_of(".") + 1) == "ply" |
		filestr.substr(filestr.find_last_of(".") + 1) == "obj")
	{
		VolMesh(volfile);
	}

	return;
}

void VolMesh(const char* volfile)
{
	cout << "the input is a mesh file!" << endl;
	cout << "Reading the input mesh!" << endl;
	TriangleMesh trimesh;
	open3d::io::ReadTriangleMesh(volfile, trimesh);

	/*if (!trimesh.IsWatertight()) {
		cout << "The mesh is not watertight, and the volume may be not accurate!" << endl;
	}

	if (!trimesh.IsOrientable()) {
		cout << "The mesh is not orientable, and the volume may be not accurate!" << endl;
	}*/

	Eigen::Vector3d minconer = trimesh.GetMinBound();
	auto GetSignedVolumeOfTriangle = [&](size_t tidx) {
		const Eigen::Vector3i& triangle = trimesh.triangles_[tidx];
		const Eigen::Vector3d& vertex0 = trimesh.vertices_[triangle(0)] - minconer;
		const Eigen::Vector3d& vertex1 = trimesh.vertices_[triangle(1)] - minconer;
		const Eigen::Vector3d& vertex2 = trimesh.vertices_[triangle(2)] - minconer;
		return vertex0.dot(vertex1.cross(vertex2)) / 6.0;
	};

	double volume = 0;
	int64_t num_triangles = trimesh.triangles_.size();
	cout << "the number of triangles in mesh: " << num_triangles << endl;
#pragma omp parallel for reduction(+ : volume) num_threads(utility::EstimateMaxThreads())
	for (int64_t tidx = 0; tidx < num_triangles; ++tidx) {
		volume += GetSignedVolumeOfTriangle(tidx);
	}
	cout << "the volume of the mesh is: " << std::abs(volume)
		<< " (The value may not be accurate if the mesh is not Watertight!)" << endl;

}

// N is the number of points
vector<int> findboundary(delaunator::Delaunator d, int N)
{
	vector<int> boundaryIdices;
	vector<vector<int>> labels(N); // store the position indices of each point in d

	for (int i = 0; i < d.triangles.size(); i++)
	{
		int idx = d.triangles[i];

		labels[idx].push_back(i);
	}


	for (int i = 0; i < labels.size(); i++)
	{
		// check whether it is boundary
		int faceno = labels[i].size();
		vector<int> posvec;
		for (int j = 0; j < faceno; j++)
		{
			int pos = labels[i][j];

			if (pos % 3 == 0)
			{
				posvec.push_back(d.triangles[pos + 1]);
				posvec.push_back(d.triangles[pos + 2]);
			}
			else if (pos % 3 == 1)
			{
				posvec.push_back(d.triangles[pos - 1]);
				posvec.push_back(d.triangles[pos + 1]);
			}
			else // pos % 3 == 2
			{
				posvec.push_back(d.triangles[pos - 2]);
				posvec.push_back(d.triangles[pos - 1]);
			}
		}

		sort(posvec.begin(), posvec.end());
		int count = std::distance(posvec.begin(), std::unique(posvec.begin(), posvec.begin() + posvec.size()));

		// it is a boundary vertex
		if (count != faceno)
		{
			boundaryIdices.push_back(i);
		}
	}

	return boundaryIdices;

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

void PaintMesh(open3d::geometry::TriangleMesh& mesh,
	const vector<Eigen::Vector3d>& color) {
	mesh.vertex_colors_.resize(mesh.vertices_.size());
	for (size_t i = 0; i < mesh.vertices_.size(); i++) {
		mesh.vertex_colors_[i] = color[i];
	}
}

bool planefitABC(vector<Eigen::Vector3d> points, double n[4])
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

double dotproduct(Eigen::Vector3d v1, Eigen::Vector3d v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void ProjDisPoint(vector<Eigen::Vector3d> points, double* n, vector<Eigen::Vector3d>& Projpoints, vector<double>& distoplane)
{
	// n is the plane formulation
	Eigen::Vector3d normal(n[0], n[1], n[2]);
	double u = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];

	for (int i = 0; i < points.size(); i++)
	{
		double v = dotproduct(normal, points[i]) + n[3];

		double dis = std::abs(v) / std::sqrt(u);
		distoplane.push_back(dis);

		double t = v / u;

		Projpoints.push_back(Eigen::Vector3d(points[i][0] - n[0] * t, points[i][1] - n[1] * t,
			points[i][2] - n[2] * t));
	}
}

double trianglearea(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
{
	Eigen::Vector3d AB = B - A;
	Eigen::Vector3d AC = C - A;

	Eigen::Vector3d crossproduct;
	crossproduct[0] = AB[1] * AC[2] - AB[2] * AC[1];
	crossproduct[1] = -(AB[0] * AC[2] - AB[2] * AC[0]);
	crossproduct[2] = AB[0] * AC[1] - AB[1] * AC[0];

	double area = 0.5 * std::sqrt(crossproduct[0] * crossproduct[0] +
		crossproduct[1] * crossproduct[1] + crossproduct[2] * crossproduct[2]);

	return area;
}

void VolPointCloud(const char* volfile)
{
	cout << "the input is a point cloud file !" << endl;
	
	PointCloud pointcloud;
	vector<Eigen::Vector3d> bbox(2);
	vector<Eigen::VectorXd> MetaData;
	vector<ExtraDim> extraDims;
	vector<Eigen::Vector3d> point3D = ReadLas(volfile, bbox, MetaData, extraDims);
	PointCloudAI* pCloud = new PointCloudAI(point3D, bbox);

	cout << "the point cloud has " << point3D.size() << " points!" << endl;

	cout << "Meshing the point cloud is running!" << endl;

	/*Delanuary Triangulation*/
	std::vector<double> coords;
	for (int i = 0; i < pCloud->m_point3D.size(); i++)
	{
		coords.push_back(pCloud->m_point3D[i][0]);
		coords.push_back(pCloud->m_point3D[i][1]);
	}

	std::clock_t start;
	double duration;
	start = std::clock(); // get current time
	delaunator::Delaunator d(coords);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "Meshing took " << duration << "seconds" << std::endl;

	/*writing triangulation into a mesh (.*ply)*/
	open3d::geometry::TriangleMesh trimesh;
	for (std::size_t i = 0; i < d.triangles.size(); i += 3) {

		Eigen::Vector3i aa(d.triangles[i], d.triangles[i + 1], d.triangles[i + 2]);
		trimesh.triangles_.push_back(aa);
	}
	trimesh.vertices_ = pCloud->m_point3D;

	const char* SavedTriName = CombineFileName(volfile, "_MESH.ply");

	open3d::io::WriteTriangleMesh(SavedTriName, trimesh);
	
	/*boundary detection of the mesh above*/
	std::vector<int> boundaryindices;

	start = std::clock(); // get current time
	boundaryindices = findboundary(d, coords.size());
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "Boundary detection took " << duration << "seconds" << std::endl;

	/*rewrite a ply file with boundary points as red */
	vector<Eigen::Vector3d> colors;
	for (int i = 0; i < trimesh.vertices_.size(); i++)
	{
		colors.push_back(Eigen::Vector3d(1, 1, 1));
	}

	// give red to the boundarypoints
	for (int i = 0; i < boundaryindices.size(); i++)
	{
		colors[boundaryindices[i]] = Eigen::Vector3d(1, 0, 0);
	}

	PaintMesh(trimesh, colors);
	SavedTriName = CombineFileName(volfile, "_toes.ply");

	vector<Eigen::Vector3d> boundarypoints;
	for (int i = 0; i < boundaryindices.size(); i++)
	{
		boundarypoints.push_back(trimesh.vertices_[boundaryindices[i]]);
	}

	//open3d::io::WriteTriangleMesh(SavedTriName, trimesh);
	// save toes as a las file
	vector<ExtraDim> extraDimsbound;
	vector<Eigen::VectorXd> MetaDatabound;
	for (int i = 0; i < boundaryindices.size(); i++)
	{
		if (!extraDims.empty())
		{
			extraDimsbound.push_back(extraDims[boundaryindices[i]]);
		}
		
		MetaDatabound.push_back(MetaData[boundaryindices[i]]);
	}

	const char* SavedLasName = CombineFileName(volfile, "_toes.las");

	SaveLas(SavedLasName, boundarypoints, MetaDatabound, extraDimsbound);



	cout << "Fitting a plane to the boundary points!" << endl;

	// fit plane with boundary points
	double N[4];
	planefitABC(boundarypoints, N);

	Eigen::Vector3d normals(N[0], N[1], N[2]);

	vector <double> distances;

	for (int j = 0; j < boundarypoints.size(); j++)
	{
		distances.push_back(std::abs(boundarypoints[j].dot(normals) + N[3]));
	}

	double roughness = stdvariance(distances);

	cout << "the std derivation of the distances (boundary points to the fit plane): " << roughness << endl;

	/* draw the fitting plane */
	// plot a mesh plane
	vector<Eigen::Vector3d> Cornerpoints;
	vector<Eigen::Vector3d> CornerProjpoints;
	vector<double> Cornerdistoplane;

	Eigen::Vector3d maxcorner = trimesh.GetMaxBound();
	Eigen::Vector3d mincorner = trimesh.GetMinBound();

	Eigen::Vector3d aa;
	aa[0] = mincorner[0];
	aa[1] = mincorner[1];
	aa[2] = maxcorner[2];
	Cornerpoints.push_back(aa);

	aa[0] = mincorner[0];
	aa[1] = maxcorner[1];
	aa[2] = maxcorner[2];
	Cornerpoints.push_back(aa);

	aa[0] = maxcorner[0];
	aa[1] = mincorner[1];
	aa[2] = maxcorner[2];
	Cornerpoints.push_back(aa);

	aa[0] = maxcorner[0];
	aa[1] = maxcorner[1];
	aa[2] = maxcorner[2];
	Cornerpoints.push_back(aa);

	ProjDisPoint(Cornerpoints, N, CornerProjpoints, Cornerdistoplane);

	// save triangle into a plane mesh
	open3d::geometry::TriangleMesh planemesh;
	Eigen::Vector3i bb;
	bb[0] = 0;
	bb[1] = 1;
	bb[2] = 2;
	planemesh.triangles_.push_back(bb);
	bb[0] = 1;
	bb[1] = 2;
	bb[2] = 3;
	planemesh.triangles_.push_back(bb);
	planemesh.vertices_ = CornerProjpoints;

	SavedTriName = CombineFileName(volfile, "_PlaneMesh.ply");
	open3d::io::WriteTriangleMesh(SavedTriName, planemesh);

	cout << "vol calculation is running!" << endl;

	// compute projection points and distance to plane for all points
	vector<Eigen::Vector3d> Projpoints;
	vector<double> Distoplane;
	vector<Eigen::Vector3d> Vertices = pCloud->m_point3D;
	ProjDisPoint(Vertices, N, Projpoints, Distoplane);

	// compute vol;
	double vol = 0;
	for (std::size_t i = 0; i < d.triangles.size(); i += 3)
	{
		double area = trianglearea(Projpoints[d.triangles[i]], Projpoints[d.triangles[i + 1]], Projpoints[d.triangles[i + 2]]);
		double dissum = Distoplane[d.triangles[i]] + Distoplane[d.triangles[i + 1]] +
			Distoplane[d.triangles[i + 2]];
		vol += 1.0 / 3 * area * dissum;
	}

	cout << "the volume calculated is: " << vol << endl;

	cout << "The mesh file (ply), toes (las) and the fit plane (ply) have been saved in the working directory!" << endl;
}


void Volcal(const char* groundfile, const char* ceilfile)
{
	if (!is_file_exist(groundfile))
	{
		cout << "the first input file not found!" << endl;
		return;
	}

	if (!is_file_exist(ceilfile))
	{
		cout << "the second input file not found!" << endl;
		return;
	}

	// read files
	vector<Eigen::Vector3d> groundbbox(2);
	vector<Eigen::VectorXd> MetaData;
	vector<ExtraDim> extraDims;
	vector<Eigen::Vector3d> groundpoints = ReadLas(groundfile, groundbbox, MetaData, extraDims);

	vector<Eigen::Vector3d> ceilbbox(2);
	vector<Eigen::Vector3d> ceilpoints = ReadLas(ceilfile, ceilbbox, MetaData, extraDims);

	// 1. find the minicorner  where we start to grid point cloud
	Eigen::Vector3d miniCorner;
	miniCorner[0] = std::min(groundbbox[0][0], ceilbbox[0][0]);
	miniCorner[1] = std::min(groundbbox[0][1], ceilbbox[0][1]);
	miniCorner[2] = std::min(groundbbox[0][2], ceilbbox[0][2]);

	Eigen::Vector3d maxiCorner;
	maxiCorner[0] = std::max(groundbbox[1][0], ceilbbox[1][0]);
	maxiCorner[1] = std::max(groundbbox[1][1], ceilbbox[1][1]);
	maxiCorner[2] = std::max(groundbbox[1][2], ceilbbox[1][2]);
	double gridStep = 0.01;
	double eachStep = 0.01;
	double finalgridStep = 0.01;
	VolumeResults result;
	double iniratio;
	VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, gridStep, result, iniratio);

	if (iniratio >= NIEGH_RATIO)
	{
		double ratio = iniratio;
		while (ratio >= NIEGH_RATIO)
		{
			gridStep = gridStep /2;
			VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, gridStep, result, ratio);
		}
		finalgridStep = gridStep * 2;
	}
	else
	{
		double ratio = iniratio;
		while (ratio < NIEGH_RATIO)
		{
			gridStep = gridStep + eachStep;
			VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, gridStep, result, ratio);
		}
		finalgridStep = gridStep;
	}
	
	double neighbourratio = 0;

	VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, finalgridStep, result,
		          neighbourratio, FillEmptyStrategy::INTERPOLATE);


	cout << " the optimal girdstep is " << finalgridStep << endl;

	cout << "---------------Volume Info---------------" << endl;
	cout << " the volume difference is " << result.VolumeDiff << endl;
	cout << " the added volume is " << result.addedVolume << endl;
	cout << " the removed volume is " << result.removedVolume << endl;
	cout << " the surface is " << result.Surface << endl;
	cout << " the average number of neighbours for cells is " << neighbourratio << endl;
	cout << " the matching percentage between ground and ceil is " << result.matchPercentage << "%" << endl;
	cout << " the non-matching percentage ground is " << result.nonmatchGroundPercentage << "%" << endl;
	cout << " the non-matching percentage ceil is " << result.nonmatchCeilPercentage << "%" << endl;
	const char* SavedReportName = CombineFileName(groundfile, "_VolChangeReport.txt");
	SaveResult(SavedReportName, finalgridStep, result, true);
	cout << " the volume change result has been saved in " << SavedReportName << endl;

}


// volume change calculation for two scans with a given gridstep
void Volcal(const char* groundfile, const char* ceilfile, double gridstep)
{
	if (!is_file_exist(groundfile))
	{
		cout << "the first input file not found!" << endl;
		return;
	}

	if (!is_file_exist(ceilfile))
	{
		cout << "the second input file not found!" << endl;
		return;
	}

	// read files
	vector<Eigen::Vector3d> groundbbox(2);
	vector<Eigen::VectorXd> MetaData;
	vector<ExtraDim> extraDims;
	vector<Eigen::Vector3d> groundpoints = ReadLas(groundfile, groundbbox, MetaData, extraDims);

	vector<Eigen::Vector3d> ceilbbox(2);
	vector<Eigen::Vector3d> ceilpoints = ReadLas(ceilfile, ceilbbox, MetaData, extraDims);

	// 1. find the minicorner  where we start to grid point cloud
	Eigen::Vector3d miniCorner;
	miniCorner[0] = std::min(groundbbox[0][0], ceilbbox[0][0]);
	miniCorner[1] = std::min(groundbbox[0][1], ceilbbox[0][1]);
	miniCorner[2] = std::min(groundbbox[0][2], ceilbbox[0][2]);

	Eigen::Vector3d maxiCorner;
	maxiCorner[0] = std::max(groundbbox[1][0], ceilbbox[1][0]);
	maxiCorner[1] = std::max(groundbbox[1][1], ceilbbox[1][1]);
	maxiCorner[2] = std::max(groundbbox[1][2], ceilbbox[1][2]);
	
	double neighbourratio = 0;

	VolumeResults result;
	if (VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, gridstep, result,
		neighbourratio, FillEmptyStrategy::INTERPOLATE, true))
	{
		cout << " the given girdstep is " << gridstep << endl;

		cout << "---------------Volume Info---------------" << endl;
		cout << " the volume difference is " << result.VolumeDiff << endl;
		cout << " the added volume is " << result.addedVolume << endl;
		cout << " the removed volume is " << result.removedVolume << endl;
		cout << " the surface is " << result.Surface << endl;
		cout << " the average number of neighbours for cells is " << neighbourratio << endl;
		cout << " the matching percentage between ground and ceil is " << result.matchPercentage << "%" << endl;
		cout << " the non-matching percentage ground is " << result.nonmatchGroundPercentage << "%" << endl;
		cout << " the non-matching percentage ceil is " << result.nonmatchCeilPercentage << "%" << endl;
		
		const char* SavedReportName = CombineFileName(groundfile, "_VolChangeReport.txt");
		SaveResult(SavedReportName, gridstep, result);
		cout << " the volume change result has been saved in " << SavedReportName << endl;
	}

}


bool VolumeDiffCal(vector<Eigen::Vector3d> groundpoints, vector<Eigen::Vector3d> ceilpoints,
	Eigen::Vector3d miniCorner, Eigen::Vector3d maxiCorner,
	double gridStep, VolumeResults &result, double &ratio, 
	FillEmptyStrategy fillstrategy, bool fixedgrid/*false*/, int projtype/*=0*/)
{

	if (projtype != 0 && projtype != 1 && projtype != 2)
	{
		cout << "the last parameter should be 0 or 1 or 2" << endl;
		cout << "  0 : use the average height for each cell!" << endl; // default 
		cout << "  1 : use the minimum height for each cell!" << endl;
		cout << "  2 : use the maximum height for each cell!" << endl;
		return false;
	}

	ProjectionType projectiontype;

	switch (projtype)
	{
	case 0:
		projectiontype = ProjectionType::AVERAGE;
		break;
	case 1:
		projectiontype = ProjectionType::MINIMUM;
		break;
	case 2:
		projectiontype = ProjectionType::MAXIMUM;
		break;
	default:
		break;
	}
	

	// 2. calculate the grid size for minicorner and maxicorner: gridWidth and gridHeight
	// we only consider X and Y direction
	unsigned int gridWidth = 1 + static_cast<unsigned>((maxiCorner[0] - miniCorner[0]) / gridStep + 0.5);
	unsigned int gridHeight = 1 + static_cast<unsigned>((maxiCorner[1] - miniCorner[1]) / gridStep + 0.5);
	
	//cout << " the grid size is " << gridHeight << "*" << gridWidth << endl;

	Cell initcell = { 0, false, 0, 0, 0, 1.0e+8, -1.0e+8, 0, 0 };

	vector<vector<Cell>> groundCells(gridHeight, vector<Cell>(gridWidth, initcell));
	vector<vector<Cell>> ceilCells(gridHeight, vector<Cell>(gridWidth, initcell));

	vector<vector<Cell>> MergeCells(gridHeight, vector<Cell>(gridWidth, initcell)); // only used for neighbour count

	/*if (gridStep <= 1.0e-6 || gridWidth == 0 || gridHeight == 0)
	{
		cout << "so small grid step!" << endl;
		return;
	}*/

	unsigned int gridSize = gridWidth * gridHeight;
	/*if (gridSize == 1)
	{
		cout << "the grid size is not expected to be 1!" << endl;
		return;
	}
	if (gridSize >= 1.0e+8)
	{
		cout << "the grid size is so big and it may cause memory issue!" << endl;
		return;
	}*/

	// 3. for each point in point1, assign the cell index to it

	for (int p = 0; p < groundpoints.size(); p++)
	{
		Eigen::Vector3d point = groundpoints[p];

		Eigen::Vector3d dist = point - miniCorner;

		int i = static_cast<int>(dist[0] / gridStep + 0.5);
		int j = static_cast<int>(dist[1] / gridStep + 0.5);

		if (i < 0 || i >= static_cast<int>(gridWidth)
			|| j < 0 || j >= static_cast<int>(gridHeight)) 
			continue;
		Cell& ncell = groundCells[j][i];
		// give this point to the vector
		++ ncell.Npoints;
		ncell.Height += point[2]; // z

		if (ncell.Npoints)
		{
			if (point[2] < groundCells[j][i].minHeight)
			{
				groundCells[j][i].minHeight = point[2];
				groundCells[j][i].minPointIndex = p;
			}

			if (point[2] > groundCells[j][i].maxHeight)
			{
				groundCells[j][i].maxHeight = point[2];
				groundCells[j][i].maxPointIndex = p;
			}

		}
		
	}


	// 4. for each point in point2, assign the cell index to it
	for (int p = 0; p < ceilpoints.size(); p++)
	{
		Eigen::Vector3d point = ceilpoints[p];

		Eigen::Vector3d dist = point - miniCorner;

		int i = static_cast<int>(dist[0] / gridStep + 0.5);
		int j = static_cast<int>(dist[1] / gridStep + 0.5);

		if (i < 0 || i >= static_cast<int>(gridWidth)
			|| j < 0 || j >= static_cast<int>(gridHeight))
			continue;
					Cell& ncell = ceilCells[j][i];
		// give this point to the vector
		ncell.Npoints += 1;
		ncell.Height += point[2]; // z


		if (ncell.Npoints)
		{
			if (point[2] < ceilCells[j][i].minHeight)
			{
				ceilCells[j][i].minHeight = point[2];
				ceilCells[j][i].minPointIndex = p;
			}

			if (point[2] > ceilCells[j][i].maxHeight)
			{
				ceilCells[j][i].maxHeight = point[2];
				ceilCells[j][i].maxPointIndex = p;
			}

		}

	}


	// 5. for corresponding cells, compute height difference

	int matchcells = 0;
	int nonmatchcells_ceil = 0;
	int nonmatchcells_ground = 0;
	int cellCount = 0;

	double volume = 0;
	double addedvolume = 0;
	double removedvolume = 0;
	double surface = 0;


	for (int j = 0; j < gridHeight; j++)
	{
		for (int i = 0; i < gridWidth; i++)
		{
			Cell &gCell = groundCells[j][i];
			Cell &cCell = ceilCells[j][i];

			switch (projectiontype)
			{
			case ProjectionType::AVERAGE:
				if (cCell.Npoints !=0)
					cCell.Height = cCell.Height / cCell.Npoints;
				if (gCell.Npoints !=0)
					gCell.Height = gCell.Height / gCell.Npoints;
				break;

			case ProjectionType::MINIMUM:
				cCell.Height = cCell.minHeight;
				gCell.Height = gCell.minHeight;
				break;

			case ProjectionType::MAXIMUM:
				cCell.Height = cCell.maxHeight;
				gCell.Height = gCell.maxHeight;
				break;

			default:
				break;

			}


			double hdiff = 0;
			// it can compute difference
			if (gCell.Npoints != 0 && cCell.Npoints != 0)
			{
				++matchcells;
				++cellCount;
				surface += 1.0;
				
				hdiff = cCell.Height - gCell.Height;
				volume += hdiff;

				MergeCells[j][i].Height = hdiff;
			}
			else
			{
				if (gCell.Npoints != 0)
				{
					++cellCount;
					++nonmatchcells_ground;
				}
				else if (cCell.Npoints != 0)
				{
					++cellCount;
					++nonmatchcells_ceil;
				}
				MergeCells[j][i].Height = std::numeric_limits<double>::quiet_NaN();
			}

			if (hdiff < 0)
			{
				removedvolume -= hdiff;
			}
			else if (hdiff > 0)
			{
				addedvolume += hdiff;
			}
			
		}
		
	}
	
	//count the average number of valid neighbors
	int validNeighborsCount = 0;
	int count = 0;
	for (unsigned j = 1; j < gridHeight - 1; ++j)
	{
		for (unsigned i = 1; i < gridWidth - 1; ++i)
		{
			Cell curCell = MergeCells[j][i];

			if (std::isfinite(curCell.Height))
			{
				for (unsigned k = j - 1; k <= j + 1; ++k)
				{
					for (unsigned l = i - 1; l <= i + 1; ++l)
					{
						if (k != j || l != i)
						{
							Cell otherCell = MergeCells[k][l];
							if (std::isfinite(otherCell.Height))
							{
								++validNeighborsCount;
							}
						}
					}
				}
				++count;
			}

		}
	}

	if (fixedgrid) // using a fixed gridstep
	{
		ratio = static_cast<double>(validNeighborsCount) / count;
		if (ratio < 7.0)
		{
			cout << "the cells are sparse. Please increase the grid step!" << endl;
			return false;
		}
	}

	// interpolate. It may cause the error in the edge. Only used for the final volume calculation 

	if (fillstrategy == FillEmptyStrategy::INTERPOLATE)
	{
		// reset all parameters
		matchcells = 0;
		nonmatchcells_ceil = 0;
		nonmatchcells_ground = 0;
		cellCount = 0;

		volume = 0;
		addedvolume = 0;
		removedvolume = 0;
		surface = 0;

		

		for (unsigned j = 1; j < gridHeight - 1; ++j)
		{
			for (unsigned i = 1; i < gridWidth - 1; ++i)
			{
				Cell &gCell = groundCells[j][i];
				Cell &cCell = ceilCells[j][i];

				if (gCell.Npoints ==0) // no point
				{
					int neighpoint = 0;
					double totalHeight = 0;
					// KNN interpolation (K=8)
					for (unsigned k = j - 1; k <= j + 1; ++k)
					{
						for (unsigned l = i - 1; l <= i + 1; ++l)
						{
							if (k != j || l != i)
							{
								Cell otherCell = groundCells[k][l];
								if (otherCell.Npoints !=0)
								{
									neighpoint++;
									totalHeight += otherCell.Height;
								}
							}
						}
					}

					if (neighpoint != 0)
					{
						gCell.Height = totalHeight / neighpoint;
						gCell.Binterpolate = true;
					}

				}

				if (cCell.Npoints == 0) // no point
				{
					int neighpoint = 0;
					double totalHeight = 0;
					// KNN interpolation (K=8)
					for (unsigned k = j - 1; k <= j + 1; ++k)
					{
						for (unsigned l = i - 1; l <= i + 1; ++l)
						{
							if (k != j || l != i)
							{
								Cell otherCell = ceilCells[k][l];
								if (otherCell.Npoints !=0)
								{
									neighpoint++;
									totalHeight += otherCell.Height;
								}
							}
						}
					}

					if (neighpoint !=0)
					{
						cCell.Height = totalHeight / neighpoint;
						cCell.Binterpolate = true;
					}

				}

			}
		}

		for (int j = 0; j < gridHeight; j++)
		{
			for (int i = 0; i < gridWidth; i++)
			{
				Cell gCell = groundCells[j][i];
				Cell cCell = ceilCells[j][i];

				double hdiff = 0;
				// it can compute difference
				if (gCell.Npoints != 0 && cCell.Npoints != 0)
				{
					++matchcells;
					++cellCount;
					surface += 1.0;
					hdiff = cCell.Height - gCell.Height;
					volume += hdiff;
				}
				else
				{
					if (gCell.Npoints != 0)
					{
						++cellCount;
						if (cCell.Binterpolate)
						{
							hdiff = cCell.Height - gCell.Height;
							volume += hdiff;
							surface += 1.0;
							++matchcells;
							//cout << "match cell here" << endl;
						}
						else
						{
							++nonmatchcells_ground;
						}	
					}
					else if (cCell.Npoints != 0)
					{
						++cellCount;
						if (gCell.Binterpolate)
						{
							hdiff = cCell.Height - gCell.Height;
							++matchcells;
							volume += hdiff;
							surface += 1.0;
						}
						else
						{
							++nonmatchcells_ceil;
						}
					}
					else { // both Npoints==0
						++cellCount;
						if (gCell.Binterpolate && cCell.Binterpolate)
						{
							hdiff = cCell.Height - gCell.Height;
							++matchcells;
							volume += hdiff;
							surface += 1.0;
						}
						else if (gCell.Binterpolate)
						{
							++nonmatchcells_ceil;
						}
						else
						{
							++nonmatchcells_ground;
						}
					}
					
				}

				if (hdiff < 0)
				{
					removedvolume -= hdiff;
				}
				else if (hdiff > 0)
				{
					addedvolume += hdiff;
				}
			}
		}
	}


	//cout << "the ratio is: " << 1.0*validNeighborsCount / count << endl;
	ratio = static_cast<double>(validNeighborsCount) / count;
	double area = gridStep * gridStep;
	result.addedVolume = addedvolume * area;
	result.removedVolume = removedvolume * area;
	result.Surface = surface * area;
	result.VolumeDiff = volume * area;
	result.matchPercentage = static_cast<float>(matchcells * 100) / cellCount;
	result.nonmatchGroundPercentage = static_cast<float>(nonmatchcells_ground *100)/ cellCount;
	result.nonmatchCeilPercentage = static_cast<float>(nonmatchcells_ceil * 100) / cellCount;
	
	return true;
}


