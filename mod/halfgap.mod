NEURON {
	POINT_PROCESS HalfGap
	ELECTRODE_CURRENT i
	RANGE g, i, vgap
}

PARAMETER {
	g = 1 (nanosiemens)
}

ASSIGNED {
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
}

BREAKPOINT {
	i = .001 * g * (vgap - v)
}

:the .001 is to fix the units to nA