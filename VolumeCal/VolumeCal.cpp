# include "VolumeCal.hpp"

#define NIEGH_RATIO 7.8

bool is_file_exist(const char* fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

void Volcal(const char* pcfile)
{
	return;
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
	vector<Eigen::Vector3d> groundpoints = ReadLas(groundfile, groundbbox);

	vector<Eigen::Vector3d> ceilbbox(2);
	vector<Eigen::Vector3d> ceilpoints = ReadLas(ceilfile, ceilbbox);

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
	//finalgridStep = 0.1;
	VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, finalgridStep, result,
		          neighbourratio, FillEmptyStrategy::INTERPOLATE);
	/*VolumeDiffCal(groundpoints, ceilpoints, miniCorner, maxiCorner, finalgridStep, result,
		neighbourratio);*/


	cout << " the final girdstep is " << finalgridStep << endl;

	cout << "---------------Volume Info---------------" << endl;
	cout << " the volume difference is " << result.VolumeDiff << endl;
	cout << " the added volume is " << result.addedVolume << endl;
	cout << " the removed volume is " << result.removedVolume << endl;
	cout << " the surface is " << result.Surface << endl;
	cout << " the average number of neighbours for cells is " << neighbourratio << endl;
	cout << " the matching percentage between ground and ceil is " << result.matchPercentage << "%" << endl;
	cout << " the non-matching percentage ground is " << result.nonmatchGroundPercentage << "%" << endl;
	cout << " the non-matching percentage ceil is " << result.nonmatchCeilPercentage << "%" << endl;

}


void VolumeDiffCal(vector<Eigen::Vector3d> groundpoints, vector<Eigen::Vector3d> ceilpoints,
	Eigen::Vector3d miniCorner, Eigen::Vector3d maxiCorner,
	double gridStep, VolumeResults &result, double &ratio, 
	FillEmptyStrategy fillstrategy, int projtype/*=0*/)
{

	if (projtype != 0 && projtype != 1 && projtype != 2)
	{
		cout << "the last parameter should be 0 or 1 or 2" << endl;
		cout << "  0 : use the average height for each cell!" << endl; // default 
		cout << "  1 : use the minimum height for each cell!" << endl;
		cout << "  2 : use the maximum height for each cell!" << endl;
		return;
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
				cCell.Height = cCell.Height / cCell.Npoints;
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
	// interpolate. It may cause the error of the edge. Only used for the final volume calculation 
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
						if (gCell.Binterpolate && cCell.Binterpolate)
						{
							hdiff = cCell.Height - gCell.Height;
							++matchcells;
							volume += hdiff;
							surface += 1.0;
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
	
	return;
}
