#include <iostream>
#include "VolumeCal.hpp"

// @ author: Cheng Li 


int main(int argc, char* argv[])
{

    if (argc != 4)
    {
        cout << "the number of the input parameters is not 3!" << endl;
        cout << "the correct form is 'VolumeCal inputcloud gridstep'" << endl;
        cout << "The gridstep 0.1 or 0.01 is recommended!" << endl;
        return 0;
    }

    cout << "the volumemetric change algorithm is running!" << endl;

    const char* inputfile1 = argv[1];
    const char* inputfile2 = argv[2];
    double gridStep = std::atof(argv[3]);
    double volDiff = 0;
    VolumeResults result;
    int projType = 0;

    VolumeDiffCal(inputfile1, inputfile2, gridStep, projType);

    system("pause");
}