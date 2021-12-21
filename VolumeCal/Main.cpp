#include <iostream>
#include "VolumeCal.hpp"

// @ author: Cheng Li 


int main(int argc, char* argv[])
{

    if (argc != 3)
    {
        cout << "the number of the input parameters is not 2!" << endl;
        cout << "the correct form is 'VolumeCal referenceCloud targetCloud'" << endl;
        
        return 0;
    }

    cout << "the volumemetric change algorithm is running!" << endl;

    const char* inputfile1 = argv[1];
    const char* inputfile2 = argv[2];
    double volDiff = 0;
    VolumeResults result;
    int projType = 0;

    Volcal(inputfile1, inputfile2);

    system("pause");
}