#include <iostream>
#include "VolumeCal.hpp"

#include <ctime>
// @ author: Cheng Li 


int main(int argc, char* argv[])
{

    if ((argc!=3) && (argc !=2) && (argc != 4))
    {
        cout << "the correct input format is 'VolumeCal file1 file2' or 'VolumeCal file1 file2 gridtsep' or ''VolumeCal file1'" << endl;
        cout << "VolumeCal file1 file2 -- calculate volume change between two input files." << endl;
        cout << "VolumeCal file1 file2 gridtsep -- calculate volume change between two input files with a given gridstep." << endl;
        cout << "VolumeCal file1 -- calculate volume for a point cloud or a mesh" << endl;
        return 0;
    }

    if (argc == 4)
    {
        cout << "the volumemetric change algorithm with the given gridstep is running!" << endl;
        const char* inputfile1 = argv[1];
        const char* inputfile2 = argv[2];
        double gridstep = atof(argv[3]);

        std::clock_t start;
        double duration;
        start = std::clock(); // get current time

        Volcal(inputfile1, inputfile2, gridstep);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        cout << "Operation took " << duration << " seconds" << endl;

        system("pause");
    }

    if (argc == 3)
    {
        cout << "the volumemetric change algorithm with the optimal gridstep is running!" << endl;
        const char* inputfile1 = argv[1];
        const char* inputfile2 = argv[2];

        std::clock_t start;
        double duration;
        start = std::clock(); // get current time

        Volcal(inputfile1, inputfile2);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        cout << "Operation took " << duration << " seconds" << endl;

        system("pause");
    }

    if (argc == 2)
    {
        
        const char* inputfile1 = argv[1];

        if (!is_file_exist(inputfile1))
        {
            cout << "the input file not found!" << endl;
            return 0;
        }

        string filestr(inputfile1);

        if ((filestr.substr(filestr.find_last_of(".") + 1) != "las") &&
            (filestr.substr(filestr.find_last_of(".") + 1) != "laz") &&
            (filestr.substr(filestr.find_last_of(".") + 1) != "ply") &&
            (filestr.substr(filestr.find_last_of(".") + 1) != "obj"))
        {
            cout << "the input file format is not supported! (only support .las, .laz, .ply, .obj)" << endl;
            return 0;
        }

        std::clock_t start;
        double duration;
        start = std::clock(); // get current time
        Volcal(inputfile1);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        cout << "Operation took " << duration << " seconds" << endl;

        system("pause");
    }

}