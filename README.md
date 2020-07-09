# irregular-wave-for-OpenFOAM-v5.x
The irregular wave codes for OpenFOAM v5.x is based on the Airy's theory. The aim is achieved by combining multiple Airy waves together.

The irregularWaveProperties should locate in the constant folder

// All the examples are ran in the "tutorials/Multiphase/interFoam/laminar/wave"

// The first validation-- only one component (The Airy wave)
wave length 300m, wave amplitude 2.5m, phase 0 
![image](https://github.com/hhkbob/irregular-wave-for-OpenFOAM-v5.x/blob/master/RemeImage/Onecomponent.jpg)

// The four components of wave
   length = [300 250 300 350], amplitude = [2.5 2.0 1.5 3.0]
![image]
