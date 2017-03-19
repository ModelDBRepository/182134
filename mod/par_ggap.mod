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

This mod file is taken from the NEURON port of Traub's thalamocortical model
added into Maciej's model by William Stacey 6/08



ENDCOMMENT

: par_ggap.mod
: This is a conductance based gap junction model rather
: than resistance because Traub occasionally likes to 
: set g=0 which of course is infinite resistance.
NEURON {
	POINT_PROCESS gGapPar
	RANGE g, i, vgap
	ELECTRODE_CURRENT i
}
PARAMETER { g = 1e-10 (1/megohm) }
ASSIGNED {
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
}
BREAKPOINT { i = (vgap - v)*g }
