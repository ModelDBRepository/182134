TITLE IKDR CA1

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
	SUFFIX kdrcurrent
	NONSPECIFIC_CURRENT ik
	RANGE g, e
}
 
PARAMETER {
        v		(mV)
        celsius		(degC)
        g= 0.010		(mho/cm2)
        e= -90		(mV)
}
 
STATE {
	n 
}
 
ASSIGNED {
	ik		(mA/cm2) 
 	ninf
	ntau    (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ik=g*n*(v-e)      
}
 
DERIVATIVE states { 
       rates(v)
       n'= (ninf- n)/ ntau
}

INITIAL { 
	rates(v)
	n= ninf
}


PROCEDURE rates(v (mV)) {
LOCAL  a, b
UNITSOFF
a = exp(-0.11*(v-13))
b = exp(-0.08*(v-13)) 	
	ntau=50*b/(1+a)
	if (ntau<2) {ntau=2}
	ninf=1/(1+a)
UNITSON
}

