@author:Cheng Li

Three packages are required:
    pdal is used to read and write las file
    open3d is used for neighbour search and normal calculation
    armadillo is used to find a fitted plane for a group of points

1. Note that in CMakeList.txt, add_compile_definitions(NOMINMAX) must be there. I already added there.
2. Configure cmake gui
a. download pdal. I have installed OSGeo4W (https://www.usna.edu/Users/oceano/pguth/md_help/html/osgeo4w.htm). 
   Next configure three pdal variables in cmake gui. Please see the screenshot cmakeguide.png
b. download open3d opensource from github (https://github.com/isl-org/Open3D) and compile open3d to generate a installation folder. 
   For example, I have installed it in 'D:\work\open3d\Open3D\open3d_install'
   Next configure three open3d variables in cmake gui. Please see the screenshot cmakeguide.png
   An alternative method is to find open3d path in anaconda env rather than compiling open3d from scratch
c. download armadillo-10.7.3 (https://fossies.org/linux/misc/armadillo-10.7.3.tar.xz/). 
   NExt configure two related variables in cmake gui. Please see the screenshot cmakeguide.png
3. Open sln to build the project. It only works at 'Release' version. 
   

Examples:
1. CSF for ground filtering. An example is:
 - CSF ../../data/CSFutm.laz
2. VolumeCal for the volume calculation between two time epochs. 
   An example with the automatic gridstep is:
 - VolumeCal ../../data/F5_March_2019_1cm_crop_crop_small_small.las ../../data/F5_October_2019_1cm_crop_crop_register_small_small.las
   An example with a given gridstep is:
 - VolumeCal ../../data/F5_March_2019_1cm_crop_crop_small_small.las ../../data/F5_October_2019_1cm_crop_crop_register_small_small.las 0.1
3. M3C2 for change detection between two time epochs. An example is: 
 - M3C2 ../../data/F5_March_2019_1cm_crop_crop_small_small.las ../../data/F5_October_2019_1cm_crop_crop_register_small_small.las M3C2result.las
3. M3C2-EP for change detection between two time epochs. An example is: 
 - M3C2EP ../../data/2017_652700_5189300_gnd.laz ../../data/2018A_652700_5189300_gnd.laz M3C2EPresult.las
