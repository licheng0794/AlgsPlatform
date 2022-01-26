#include <iostream>
#include "VolumeCal.hpp"


// @ author: Cheng Li 


int main(int argc, char* argv[])
{

    if ((argc!=3) && (argc !=2) && (argc != 4))
    {
        cout << "the correct input format is 'VolumeCal file1 file2' or 'VolumeCal file1 file2 gridtsep' or ''VolumeCal file1'" << endl;
        
        return 0;
    }

    if (argc == 4)
    {
        cout << "the volumemetric change algorithm with the given gridstep is running!" << endl;
        const char* inputfile1 = argv[1];
        const char* inputfile2 = argv[2];
        double gridstep = atof(argv[3]);
        Volcal(inputfile1, inputfile2, gridstep);

        system("pause");
    }

    if (argc == 3)
    {
        cout << "the volumemetric change algorithm with the optimal gridstep is running!" << endl;
        const char* inputfile1 = argv[1];
        const char* inputfile2 = argv[2];
        Volcal(inputfile1, inputfile2);

        system("pause");
    }

    if (argc == 2)
    {
        cout << "the single volume calculation is running!" << endl;
        const char* inputfile1 = argv[1];
        Volcal(inputfile1);

        system("pause");
    }


   

   
}