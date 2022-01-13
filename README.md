# Central California 3D Seismic Velocity Model software


## Use with the UCVM framework

This package is intended to be installed as part of the UCVM framework (version 15.3.0 or later) and by the UCVM framework.
While it is possible to link to this library using your own C or Fortran code, we recommend using the UCVM suite of utilities.
Most common functions, such as mesh generation and query capabilities, are included with the UCVM framework.

## Use outside of the UCVM framework

The cca library may be linked into any user application.
The header file defining the API is located in src/cca.h and will be installed to $prefix/include/cca.

Note: Even when the cca library is used outside of the UCVM framework, it requires the UCVM Etree Vs30 map file `model/ucvm/ucvm.e`.
