## bcClean

This project try to solve a practical problem many encounters when doing FEM over CAD models. It is ambiguous how to set the boundary condition after mesh repair/ mesh refine. This package provide a lib and UI to transport the boundary condition labels by either nearest neighbor or Euclidean projection followed by a clean up.

### Installation
To run the, a python interpreter (>3.6) is required and python module yaml, wildmeshing, igl are necessary. We recommend using a separate anaconda env. To set up the prerequisite
```shell
conda activate env_name
conda install yaml
conda install wildmeshing
conda install igl
```
To build the binary, run the following code
```shell
git clone 
cd bcclean
mkdir build
cmake -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/dir/to/your/python/interpreter/python -DPython3_LIBRARIES=/dir/to/your/python/lib/libpython3.7m.dylib ..
make
```

### how to use
The bcclean_bin get the initial parameters from a .json file for example
```shell
cd build
./bcclean_bin -j ../run.json
```
where a standard run.json file is provided as
```jsonc
{
    "lambda_refine":1, // control the graph cut parameter lambda
    "stop_energy": 10, // control the tetwild stopping criterion
    "data_root": "../data/0", // directory to the data folder
    "num_subdiv": 2, // number of sub-division
    "label_num": 6 // label num (set to the number of patches you will use)
}
```
The detailed Usage looks like below
![](README_support/blur.gif)
