COMMENT

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
//
// Copyright 2007, The University Of Pennsylvania
// 	School of Engineering & Applied Science.
//   All rights reserved.
//   For research use only; commercial use prohibited.
//   Distribution without permission of Maciej T. Lazarewicz not permitted.
//   mlazarew@seas.upenn.edu
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENDCOMMENT
					       
NEURON {
  POINT_PROCESS fvpre
  NONSPECIFIC_CURRENT i
  RANGE gmax, g, i
  GLOBAL a, b, th, e
  POINTER vpre
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  gmax=1e-4 (uS)
  a=12 (/ms)
  b=0.1 (/ms)
  e=-75	   (mV)
  th=0 (mV) 
}

ASSIGNED { vpre (mV) v (mV) i (nA)  g (uS)}

STATE { s }

INITIAL {  s =  a*F(vpre)/(a*F(vpre)+b) }

BREAKPOINT {
  SOLVE state METHOD cnexp
  g = gmax * s
  i = g*(v - e)
}

DERIVATIVE state { s' = a*F(vpre)*(1-s) - b*s }

FUNCTION F (v1 (mV)) {  F = 1/(1 + exp(-(v1-th)/2(mV))) }  
