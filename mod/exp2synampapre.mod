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
	POINT_PROCESS Exp2SynAMPApre
	RANGE e, i
	NONSPECIFIC_CURRENT i
	POINTER pre

	RANGE gmax
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {

	e    = 0 (mV)
	gmax = 1 (uS)
}

ASSIGNED {

    pre (mV)
	v   (mV)
	i   (nA)
}

STATE { W }

INITIAL { W = 0 }

BREAKPOINT {

	SOLVE state METHOD cnexp
	
	i = gmax * W * (v - e)
}

DERIVATIVE state {

    W' = H(pre)/1(ms)-W/2(ms)    
}

FUNCTION H(x(mV)) {

    if (x>-40) {
    
        H = 1
        
    }else{
        
        H = 0
        
    }
}
