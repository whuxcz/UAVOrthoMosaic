-------- UAVOrthoMosaic v0.1 --------

-- What is UAVOrthoMosaic --

UAVOrthoMosaic provides an open source solution for generating orthomosaic from Unmanned Aerial Vehicle(UAV) images. The tool allows users to get the 3D information(point cloud), DEM and final orthomosaic from 2D UAV images without any ground control points. And the World Wind is used to visualize the result. The process of generating orthomosaic is written in C and C++, and JAVA is used to chained the whole mosaic workflow and provide a user-friend UI.

Currently, UAVOrthoMosaic has been primarily compiled and tested under Linux

-- Web Resource --

Homepage: http://geos.whu.edu.cn:8080/UAVOrthoMosaic/ <br/>
Github: https://github.com/whuxcz/UAVOrthoMosaic

-- Requirement --

Note that the process of point cloud and dense point cloud generation is based on open source software Bundler and CMVS/PMVS, so you have to have them installed on your machine! <br/>
-SIFT: http://www.cs.ubc.ca/~lowe/keypoints/ <br/>
-Bundler: http://phototour.cs.washington.edu/bundler/ <br/>
-CMVS/PMVS2: http://www.di.ens.fr/cmvs/ <br/>
-ImageMagick http://www.imagemagick.org/script/index.php

-- INSTALL --

UAVOrthoMosaic relies on several libraries and to compile the software you have to install them before compiling the software. The required libraries are: <br/>
-Eigen http://eigen.tuxfamily.org/index.php?title=Main_Page <br/>
-Flann http://www.cs.ubc.ca/research/flann/ <br/>
-Boost http://www.boost.org/ <br/>
-gzip http://www.gzip.org/ <br/>
-libjpeg http://libjpeg.sourceforge.net/ <br/>
-libtiff http://www.libtiff.org/ <br/>
-libgeotiff http://trac.osgeo.org/geotiff/

After install those libraries, modify the Makefile to link those libraries. Then cd UAVOrthoMosaic make After that, you can run the launch.sh to run the UAVOrthoMosaic!

-- Reference --

[1] Noah Snavely, Steven M. Seitz, Richard Szeliski. Photo Tourism: Exploring image collections in 3D. ACM Transactions on Graphics (Proceedings of SIGGRAPH 2006), 2006.

[2] Noah Snavely, Steven M. Seitz, Richard Szeliski. Modeling the World from Internet Photo Collections. International Journal of Computer Vision, 2007.

[3] Yasutaka Furukawa, Brian Curless, Steven M. Seitz, and Richard Szeliski Towards Internet-scale Multi-view Stereo Computer Vision and Pattern Recognition, 2010.

[4] Yasutaka Furukawa and Jean Ponce Accurate, Dense, and Robust Multi-View Stereopsis IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 32, Issue 8, Pages 1362-1376, August 2010.
