# Virtual scanner for converting 3D model to point cloud

This folder contains the code for converting the 3D models to dense point clouds with normals. As detailed in our paper, we build a virtual scanner and shoot rays to calculate the intersection point and oriented normal. 

The code is based on the [CGAL](http://www.cgal.org/) and the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. After configuring these two libraries properly, the code can be built with vs2015 easily.


If you use our code, please cite our paper.

    @article {Wang-2017-OCNN,
        title     = {O-CNN: Octree-based Convolutional Neural Networks for 3D Shape Analysis},
        author    = {Wang, Peng-Shuai and Liu, Yang and Guo, Yu-Xiao and Sun, Chun-Yu and Tong, Xin},
        journal   = {ACM Transactions on Graphics (SIGGRAPH)},
        volume    = {36},
        number    = {4},
        year      = {2017},
    }

