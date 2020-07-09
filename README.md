# irregular-wave-for-OpenFOAM-v5.x
The irregular wave codes for OpenFOAM v5.x is based on the Airy's theory. The aim is achieved by combining multiple Airy waves together.

The irregularWaveProperties should locate in the constant folder

// All the examples are ran in the "tutorials/Multiphase/interFoam/laminar/wave"

// The first validation-- only one component (The Airy wave)
wave length 300m, wave amplitude 2.5m, phase 0 

![image](https://github.com/hhkbob/irregular-wave-for-OpenFOAM-v5.x/blob/master/RemeImage/Onecomponent.jpg)

// The four components of wave
   length = [300 250 300 350], amplitude = [2.5 2.0 1.5 3.0], phase =[0 0 0 0]

![image](https://github.com/hhkbob/irregular-wave-for-OpenFOAM-v5.x/blob/master/RemeImage/componet4.jpg)

// the 16 components of wave (mesh doest not change)

length = [300 250 300 350 100 150 400 420 320 360 423 380 311.2 341.3 281.7 265.9];

amplitude = [0.51 0.25 0.62 0.21 0.35 0.08 0.14 1.22 1.021 1.05 0.86 0.879 1.043 1.24 1.30 0.66];

phase = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
![image](https://github.com/hhkbob/irregular-wave-for-OpenFOAM-v5.x/blob/master/RemeImage/component16.jpg)
