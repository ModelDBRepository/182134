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
	POINT_PROCESS Exp2SynNMDA
	RANGE tau1, tau2, tau12, tau22, e, i, i2
	NONSPECIFIC_CURRENT i, i2

	RANGE S, S2, Smax, gmax, total, r
	GLOBAL c
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {

	tau1  =  15 (ms) <1e-9,1e9>
	tau2  = 150 (ms) <1e-9,1e9>	
:	tau12 =  0.05  (ms) <1e-9,1e9>
	tau12 =  0.2  (ms) <1e-9,1e9>
:	tau22 =  5.3  (ms) <1e-9,1e9>
	tau22 =  1.0  (ms) <1e-9,1e9>
	e     = 0	(mV)
	mg    = 1
	Smax  = 125
	c     = 1
:	c     = 1.53775
	r     = 1
}

ASSIGNED {

	v       (mV)
	i       (nA)
	i2      (nA)
	S       (1)
	S2      (1)
	gmax    (uS)
	mgblock (1)
	factor  (1)
	factor2 (1)
	total   (1)
}

STATE {

	A  (1)
	B  (1)
	A2 (1)
	B2 (1)
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
	
	if (tau12/tau22 > .9999) {
		tau12 = .9999*tau22
	}
	A2 = 0
	B2 = 0	
	tp = (tau12*tau22)/(tau22 - tau12) * log(tau22/tau12)
	factor2 = -exp(-tp/tau12) + exp(-tp/tau22)
	factor2 = 1/factor2
	total = 0
}

BREAKPOINT {
	
	SOLVE state METHOD cnexp
	
	S  = B  - A
	S2 = B2 - A2
	
	:if (S>=Smax) {S = Smax}     
	
	mgblock = 1.0 / (1.0 + 0.28 * exp(-0.062(/mV) * v) )
	i  = gmax *  r *  S  * (v - e) * mgblock
	i2 = gmax *       S2 * (v - e)
}

DERIVATIVE state {

	A'  = -A/tau1
	B'  = -B/tau2	
	A2' = -A2/tau12
	B2' = -B2/tau22
}

NET_RECEIVE(w(uS)) {
	: the factor c comes from the fact that the synaptic mechanism described in the paper
	: is using the presynaptic membrane potential continuously to define the differential equation
	: and in NEURON these synapses are implemented by event based synapses
	if (S<Smax) { 
	
		state_discontinuity(A, A + factor*c)
		state_discontinuity(B, B + factor*c)
	}
		
	state_discontinuity(A2, A2 + factor2)
	state_discontinuity(B2, B2 + factor2)
	gmax = w
	total = total + 1
}
