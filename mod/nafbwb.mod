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
	SUFFIX Nafbwb
	USEION na WRITE ina
	GLOBAL phih
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gna  = 35 (mS/cm2)
    ena  = 55 (mV)
    phih = 5
}
    
ASSIGNED {
    v       (mV)
    ina     (mA/cm2)
    minf    (1)
    hinf    (1)
    taoh    (ms)
	celsius (degC)
}

STATE { h }

INITIAL {
    rates(v)
    h    = hinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ina = (1e-3) * gna * minf^3 * h * (v-ena)
}

DERIVATIVE states { 
    rates(v)
    h' = (hinf-h)/taoh
}

PROCEDURE rates(v(mV)) { LOCAL am, bm, ah, bh, q10
    
    q10  = phih:^((celsius-27.0(degC))/10.0(degC))	
    
    am   = fun3(v,  -35, -0.1,    -10)
    bm   = fun1(v,  -60,  4,      -18) 
    minf = am/(am+bm)
 
    ah   = fun1(v,  -58,    0.07,  -20)
    bh   = fun2(v,  -28,    1,     -10)
    hinf = ah/(ah+bh)
    taoh = 1./((ah+bh)*q10)
}

INCLUDE "aux_fun.inc"