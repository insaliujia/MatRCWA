# RCWA
An electromagnetic simulation tool programmed in Matlab by using the method Rigorous Coupled-Wave analysis (RCWA), developed by [Jia LIU].

For the detailed derivation please consult this: 

Moharam, M. G., & Gaylord, T. K. (1981). Rigorous coupled-wave analysis of planar-grating diffraction. JOSA, 71(7), 811-818.

Rumpf, R. C. (2011). Improved formulation of scattering matrices for semi-analytical methods that is consistent with convention. Progress In Electromagnetics Research B, 35, 241-261.

Liscidini, M., Gerace, D., Andreani, L. C., & Sipe, J. E. (2008). Scattering-matrix analysis of periodically patterned multilayers with asymmetric unit cells and birefringent media. Physical Review B, 77(3), 035324.

# User manual
The program use object oriented programming paradigm in Matlab without using any toolboxes. So only a valide version of Matlab is needed to use it. 

## Installation
No installation is need but the path has to be added in Matlab as shown in the examples.

## How to use 

There examples  have been given in the exmaple file. Basic simulation needs will be satisfied by changing parameters in the examples. In case, some explanations are given in this file.

Four Objects are mainly used: RCWA, Source, Device and Material.

  - RCWA is the main object to control the RWA calculation engine.
  - Source is used to define the illumination source. Here only plane wave source can be used. The wavelength, polarization and illumination angle can be controled.
  - Device is used to define the details of the simulated structure. In this program different shapes (cylinder, rectangle etc.) can be included.
  - Material is a seperate object to control materials used in simulation. Users have to include reflective index according to a certain format as shown in the material folder. In addition simple reflective index can be defined directly as the example shown in defining the reflective index of Air.

License
----

MIT

   [Jia LIU]: <http://l-j.xyz/>

