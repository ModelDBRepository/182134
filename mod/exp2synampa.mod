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
	POINT_PROCESS Exp2SynAMPA
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE S, total
	GLOBAL c
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {

	tau1  =   1 (ms) <1e-9,1e9>
	tau2  =   2 (ms) <1e-9,1e9>
	e     = 0	(mV)
	c     = 16.5
}

ASSIGNED {

	v      (mV)
	i      (nA)
	S      (uS)
	factor (1)
	total  (1)
}

STATE {

	A (uS)
	B (uS)
}

INITIAL {

	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
    tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	total = 0
}

BREAKPOINT {

	SOLVE state METHOD cnexp
	S = B - A
	
	i = S * (v - e)
}

DERIVATIVE state {

	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(w (uS) ) {

	: the factor c comes from the fact that the synaptic mechanism described in the paper
	: is using the presynaptic membrane potential continuously to define the differential equation
	: and in NEURON these synapses are implemented by event based synapses
	state_discontinuity(A, A + factor*c*w)
	state_discontinuity(B, B + factor*c*w)
	total = total + 1
}
