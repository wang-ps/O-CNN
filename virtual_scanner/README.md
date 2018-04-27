# Virtual scanner for converting 3D model to point cloud

This folder contains the code for converting the 3D models to dense point clouds with normals (\*.points). As detailed in our paper, we build a virtual scanner and shoot rays to calculate the intersection point and oriented normal. 

The code is based on [Boost](https://www.boost.org/), [CGAL](http://www.cgal.org/) and the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) libraries. After configuring these three libraries properly, the code can be built with visual studio easily.

`Note`: Sometimes, the executive file might collapse when the scale of the mesh is very large. This is one bug of CGAL. In order to mitigate this you can run VirtualScanner with the normalize flag set to 1.

## Running Virtual Scanner
The pre-built executive file is contained in the folder `exe`, which has been test on the Win10 x64 system. 

    Usage:  
	cout << "Usage: VirtualScanner.exe <file name/folder name> "
		"[view_num] [flags] [normalize]" << endl;
        VirtualScanner.exe <file_name> [nviews] [flags] [normalize]
            file_name: the name of the file (*.obj; *.off) to be processed.
            nviews: the number of views for scanning. Default: 6
            flags: Indicate whether to output normal flipping flag. Default: 0
            normalize: Indicate whether to normalize input mesh. Default: 0
    Example:
        VirtualScanner.exe input.obj 14         // process the file input.obj
        VirtualScanner.exe D:\data\ 14          // process all the obj/off files under the folder D:\Data

The result is in the format of `points`, which can be parsed with the following code:

```cpp
void load_points(vector<float>& pts, vector<float>& normals, const string& filename)
{
    ifstream infile(filename, ios::binary);

    int n;
    infile.read((char*)(&n), sizeof(int));
    pts.resize(3 * n);      // x_1, y_1, z_1, ..., x_n, y_n, z_n
    infile.read((char*)pts.data(), sizeof(float)*3*n);
    normals.resize(3 * n);  // nx_1, ny_1, nz_1, ..., nx_n, ny_n, nz_n
    infile.read((char*)normals.data(), sizeof(float)*3*n);

    infile.close();
}
```
## Building On Windows
To build in Windows you can use [Vcpkg](https://github.com/Microsoft/vcpkg) to install/build all the dependencies. Note this takes a long time.
```
git clone https://github.com/Microsoft/vcpkg
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install
.\vcpkg install cgal eigen3 boost --triplet x64-windows
```
Then to build, you can use the supplied solution file VirtualScanner.sln

## Building On Ubuntu
To build with ubuntu, you can use apt for the dependencies.
```
apt-get install -y --no-install-recommends libboost-all-dev libcgal-dev libeigen3-dev
```
Then you can use g++ to build the executable
```
g++ main.cpp VirtualScanner.cpp  -O3 -DNDEBUG --std=c++11 -fPIC -Wall -Wno-sign-compare -Wno-uninitialized -fopenmp -lboost_filesystem -lboost_system -lCGAL -I /usr/include/eigen3 -o virtualscanner
```


If you use our code, please cite our paper.

    @article {Wang-2017-OCNN,
        title     = {O-CNN: Octree-based Convolutional Neural Networks for 3D Shape Analysis},
        author    = {Wang, Peng-Shuai and Liu, Yang and Guo, Yu-Xiao and Sun, Chun-Yu and Tong, Xin},
        journal   = {ACM Transactions on Graphics (SIGGRAPH)},
        volume    = {36},
        number    = {4},
        year      = {2017},
    }
