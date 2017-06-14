# DenseCRF for the refinement of 3D shape part segmentation result

This folder contains the code for the refinement of 3D shape part segmentation in our [O-CNN](http://wang-ps.github.io/O-CNN.html) paper. The mathematical details are presented in the Section 5.3 of our paper.

The instructions to build the code are as follows:

- Download the original [DenseCRF](http://graphics.stanford.edu/projects/drf) library first. 
- Download the code in current folder and use it to override the original [DenseCRF](http://graphics.stanford.edu/projects/drf) code.
- Follow the instructions in the [DenseCRF](http://graphics.stanford.edu/projects/drf) to finish the build, and produce the executive files `dense_learning.exe` and `dense_inference.exe`.


If you use our code, please cite our paper.

    @article {Wang-2017-OCNN,
        title     = {O-CNN: Octree-based Convolutional Neural Networks for 3D Shape Analysis},
        author    = {Wang, Peng-Shuai and Liu, Yang and Guo, Yu-Xiao and Sun, Chun-Yu and Tong, Xin},
        journal   = {ACM Transactions on Graphics (SIGGRAPH)},
        volume    = {36},
        number    = {4},
        year      = {2017},
    }

