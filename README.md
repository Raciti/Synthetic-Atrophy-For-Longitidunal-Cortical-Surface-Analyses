# Requirements

Advice for installations following the use of [cmake](https://cmake.org/download/). <\br>
After downloading the .tar file,.
``` 
tar -xvf cmake-<version>.tar
cd cmake-<version>
./bootstrap
make -j$(nproc)
sudo make install
```


To use this project we need to:
- [VTK](https://vtk.org): The Visualization Toolkit, ([download link](https://vtk.org/download/)).
- [ITK](https://itk.org): Insight Toolkit, ([download link](https://docs.itk.org/en/latest/download.html)).


For esecute the example, the file named *Synthetic-Atrophy_example.sh*, we need external tools:
-  [c3d](https://sourceforge.net/p/c3d/git/ci/master/tree/): Convert3D Medical Image Processing Tool Code
-  [slicer](https://download.slicer.org/): 3D Slicer
-  [greedy](https://sites.google.com/view/greedyreg/): Fast Deformable Registration for 3D Medical Images, at the followin [link](https://sites.google.com/view/greedyreg/installation) the installation guide
-  [imagemath](https://github.com/NIRALUser/niral_utilities/tree/master/ImageMath): a tool of NITRC (Neurologic Tools & Resources Collaboratory)
-  [holedetection]
