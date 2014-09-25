# hinet_decon

A fortran2003 code to deconvolve Hi-net velocity record by its seismometer response by using inverse filtering technique.

## System Requirements

Linux/Mac/Unix with Fortran2003 compiler. I confirmed that this program successfully works with gfortran@4.8.2 and intel fortran@14.0.1 on Mac OSX.

## Compile

* Edit makefile
* make


## Usage

    hinet_decon.x sacfile [-o outfile] [-f0 freq] [-h0 damp] [-f1 freq]
                          [-h1 damp] [-decon] [-int] [-fint freq] [-hint freq]

      '-o  outfile'   ( indicate output sac filename <sacfile.out>          )
      '-f0 freq'      ( eigen frequency of input signal <1.0>               )
      '-h0 damp'      ( damping constant of input signal <0.7>              )
      '-f1 freq'      ( eigen frequency of output signal <0.00833333>       )
      '-h1 damp'      ( damping constant of output signal <0.707>           )
      '-decon'        ( deconvolve without simulation seismometer           )
      '-int'          ( output displacement record by numerical integration )
      '-fint freq'    ( corner frequency of low-cut filter for integration  )
      '-hint damp'    ( damping constant of low-cut filter for integration  )
  
### Hint

    ./hinet_decon.x (input.sac) -o (output.sac)

 This default option will convert Hi-net type velocity seismometer response to that of (approximated) STS-2. 
 By adding "-int" option, deconvolved displacement waveform will be obtained.

 
We confirm that this program successfully works with gfortran@4.6.0 and intel fortran@13.0.1 .


## Copyright and License

Copyright (C) 2009-2014, Takuto Maeda, All rights reserved. All source codes included in this archive are released under the MIT License. 

If you write a scientific paper describing research that made substantive use of this program, please cite the following paper.

Maeda, T., K. Obara, T. Furumura, and T. Saito, 
Interference of long-period seismic wavefield observed by dense Hi-net array in Japan,
J. Geophys. Res., 116, B10303, doi:10.1029/2011JB008464, 2011.
<http://dx.doi.org/10.1029/2011JB008464>
