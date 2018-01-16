# Virtual scanner for converting 3D model to point cloud

This folder contains the code for converting the 3D models to dense point clouds with normals (\*.points). As detailed in our paper, we build a virtual scanner and shoot rays to calculate the intersection point and oriented normal. 

The code is based on the [CGAL](http://www.cgal.org/) and the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. After configuring these two libraries properly, the code can be built with vs2015 easily.

`Note`: Sometimes, the executive file might collapse when the scale of the mesh is very large. This is one bug of CGAL. Just scale the mesh into the unit bounding box, then it will be OK.

The pre-built executive file is contained in the folder `exe`, which has been test on the Win10 x64 system. 

    Usage:  
        VirtualScanner.exe <file_name> [nviews]
            file_name: the name of the file (*.obj; *.off) to be processed.
            nviews: the number of views for scanning. suggested value: 14
    Example:
        VirtualScanner.exe input.obj 14         // process the file input.obj
        VirtualScanner.exe D:\data\*.obj 14     // process all the obj files under the folder D:\Data

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

If you use our code, please cite our paper.

    @article {Wang-2017-OCNN,
        title     = {O-CNN: Octree-based Convolutional Neural Networks for 3D Shape Analysis},
        author    = {Wang, Peng-Shuai and Liu, Yang and Guo, Yu-Xiao and Sun, Chun-Yu and Tong, Xin},
        journal   = {ACM Transactions on Graphics (SIGGRAPH)},
        volume    = {36},
        number    = {4},
        year      = {2017},
    }

