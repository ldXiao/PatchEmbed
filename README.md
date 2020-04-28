## PatchEmbed

This project try to solve a practical problem many encounter when tring to do simulations over CAD models. It is ambiguous how to set the boundary condition after mesh repair/ mesh refine. This package provide a lib and UI to transport the boundary condition labels or the patchlayout which is gauranteed to preserve topology.

### Installation
To build the binary, run the following code
```shell
git clone 
cd PatchEmbed
mkdir build
cmake ..
make
```

### how to use
The bcclean_bin get the initial parameters from a .json file for example
```shell
cd build
./patchembed -d dir/you/put/obj,yml/in -u upsample_level -b backtrackthreshold, -t tracing_type
```
where a standard set of arguments is provided as
```shell
./PatchEmbed -d ../data/2 -u 0 -b 0.4 -t loop
```


