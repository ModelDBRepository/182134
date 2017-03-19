COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS noisesyn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	POINTER ptr
	RANGE g,start,spikedur,spikefreq,weight,nospike_tau,spike_tau,normalmean,normalstd,poisson_mean,gid,syn_index
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.05 (ms) <1e-9,1e9>
	tau2 = 5.3 (ms) <1e-9,1e9>
	e=0	(mV)
	: additional parameters from proc shortspikes()
	start = 10 (ms)
	nospike_tau=3.3333 (ms) :mean time between synaptic events in between network spikes (corresponds to thresh=97 with time step of 0.1 ms); this is the tau value of a negative exponential distribution
	spike_tau = 0.6667 (ms) :mean time (ms) between synaptic events during network spikes (corresponds to thresh=85 with time step of 0.1 ms)
	spikedur = 150 (ms)
	spikefreq = 2 (hz)
	normalmean = 0 (ms)
	normalstd = 4.4721 (ms) :4.4721 = sqrt(20)
	weight = 0.00053407075(uS) :standard weight of synaptic events(from shortspikes)
	poisson_mean = 0.8 :this is the mean of the poisson distribution used to modulate the weight of specific synaptic events (from shortspikes)
	gid = 0 :need to remember to change this when the noisesyn is created in hoc
	syn_index = 0 :may change this when the noisesyn is created in hoc, but don't have to because this is just used by the random number generator to generate different streams for different synapses on the same cell; right now, we only have one noisy synapse on each cell
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	t_master (ms) :this is the master time which determines when you go in and out of a network spike
	interspike_int (ms)
	t_out (ms) :this is the 'local' time which deviates from t_master by normally distributed jitter for each cell
	temp_time (ms) :this will be used to check whether the new event time over-runs t_out
	ptr
	cachedNormal :for use with generation of normal distribution
	haveCachedNormal :for use with generation of normal distribution
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	:printf("weight = %g \n",weight)
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1) :from exp2syn
	factor = -exp(-tp/tau1) + exp(-tp/tau2) :from exp2syn
	factor = 1/factor :from exp2syn
	
	setrand(gid,syn_index) :this should be set from hoc, so that different cells have different seeds; note that if this is not set, a segmentation fault will result
	cachedNormal=0 :pertains to generation of normally distributed data
	haveCachedNormal=0 :pertains to generation of normally distributed data
	
	t_master=start
	interspike_int = 1000/spikefreq-spikedur :time between the end of one network spike and the beginning of the next
	net_send(t_master + normal(normalmean,normalstd),1) :at time 't_master + normal(normalmean,normalstd)' from now (t=0), go into stat 1, which is the beginning of a network spike
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(dummy (uS)) {
	if(flag==1) { :flag==1 means you just entered the INTRAspike interval
		chweight(my_poisson(poisson_mean)*weight*factor) 							:this initiates a synaptic event..........why did I have to make this a function?
		t_master=t_master+spikedur 							:re-set t_master to the time at which synapse will exit the network spike
		t_out=t_master + normal(normalmean,normalstd) 		:this specific cell will exit the spike close to t_master, but with some stochastic deviation (assumed that normalstd is small compared to durations of spike and interspike interval)
		:printf("t=%8.1f, t_master=%8.1f, t_out=%8.1f in State 1\n",t,t_master,t_out)
		net_send(negexp(spike_tau),3) 						:go to state 3 after a certain time chosen from neg. exp. distribution
	} else if(flag==2) { :flag==2 means you are just entered the INTERspike interval
		:printf("t=%8.1f, t_master=%8.1f, t_out=%8.1f in State 2\n",t,t_master,t_out)
		t_master=t_master+interspike_int 					:re-set t_master to the time at which synapse will next enter the network spike
		t_out=t_master + normal(normalmean,normalstd)
		temp_time=negexp(nospike_tau)  :choose  a random value from the negative exponential distribution, then see if it over-runs the next 'up' state
		if(t + temp_time<t_out) { :if the next event time occurs before the beginning of the next 'up' state, then stay within the inter-spike interval (aka 'down' state) by going to state 4
			net_send(temp_time,4) 					
		} else{ :make sure you don't over-run the next 'up' state by going directly to state 1 at the time the next 'up' state is supposed to occur
			net_send(t_out-t,1) 
		}
	} else if(flag==3) { :flag==3 means you are in the midst of the INTRAspike interval, applying synaptic stimulation at a high frequency (sent here from flag==1)
		:printf("t=%8.1f, t_master=%8.1f, t_out=%8.1f in State 3\n",t,t_master,t_out)
		if(t >= t_out) { :if you have gone past t_out, then that means you are no longer in the middle of the spike, and so you must transition to state 2
			net_send(0,2) 
		} else{ :while you are still in the middle of the spike, intiate a synaptic event, and determine the time at which the next synaptic event will occur
			chweight(my_poisson(poisson_mean)*weight*factor) :initiate synaptic event
			net_send(negexp(spike_tau),3)					:stay within network spike by staying in state 3 after a certain amount of (relatively short) time
		}	
	} else if(flag==4) { :flag==4 means you are in the midst of the INTERspike interval, applying synaptic stimulation at a low frequency (sent here from flag==2)
		:printf("t=%8.1f, t_master=%8.1f, t_out=%8.1f in State 4\n",t,t_master,t_out)
		:printf("factor = %g, weight = %g \n",factor,weight)
		:if(t >= t_out) {  :if you have gone past t_out, then that you means it's time to transition to a network spike, so you must go to state 1,(the 0 tells you to do this right away)
		:	net_send(0,1) 									
		:} else{
		chweight(my_poisson(poisson_mean)*weight*factor) 						:initiate synaptic event
		temp_time=negexp(nospike_tau)  :choose  a random value from the negative exponential distribution, then see if it over-runs the next 'up' state
		if(t + temp_time<t_out) { :if the next event time occurs before the beginning of the next 'up' state, then stay within the inter-spike interval (aka 'down' state) by staying within state 4
			net_send(temp_time,4) 
		} else{ :make sure you don't over-run the next 'up' state by going directly to state 1 at the time the next 'up' state is supposed to occur
			net_send(t_out-t,1)					
		}
	}
}

:this code uses Random123, which requires NEURON 7.3
:uses nrnran123.c and nrnran123.h from http://www.neuron.yale.edu/hg/neuron/nrn/file/9d4ab20927bc/src/oc/
VERBATIM
#define VOIDCAST void** vp = (void**)(&(_p_ptr))
extern void * nrnran123_newstream(int,int);
extern void nrnran123_deletestream(void *);
extern double nrnran123_dblpick(void *);
ENDVERBATIM

PROCEDURE setrand(id1,id2) {
	VERBATIM
	VOIDCAST;
	if(*vp) {
		nrnran123_deletestream(*vp);
	} 
	*vp = nrnran123_newstream((int) _lid1,(int) _lid2);
	ENDVERBATIM
} 

FUNCTION pick() {
	VERBATIM
	VOIDCAST;
	_lpick = nrnran123_dblpick(*vp);
	ENDVERBATIM
}

PROCEDURE chweight(delta) {
	:printf("Weight change = %g \n",delta)
	A = A + delta
	B = B + delta
}

FUNCTION normal(pMean,pStdDev) { :adapted from http://www.neuron.yale.edu/hg/neuron/nrn/file/9d4ab20927bc/src/gnu/Normal.cpp
VERBATIM
if (haveCachedNormal == 1) {
	haveCachedNormal = 0;
	_lnormal = cachedNormal * _lpStdDev + _lpMean ;
} else {

		for(;;) {
			double u1 = pick();
			double u2 = pick();
			double v1 = 2 * u1 - 1;
			double v2 = 2 * u2 - 1;
			double w = (v1 * v1) + (v2 * v2);
 	    
	//
	//	We actually generate two IID normal distribution variables.
	//	We cache the one & return the other.
	// 
			if (w <= 1) {
				double y = sqrt( (-2 * log(w)) / w);
				double x1 = v1 * y;
				double x2 = v2 * y;
	
				haveCachedNormal = 1;
				cachedNormal = x2;
				_lnormal = x1 * _lpStdDev + _lpMean;
				return _lnormal;
			}
		}
    }
ENDVERBATIM
}

FUNCTION negexp(mean) { :adapted from http://www.neuron.yale.edu/hg/neuron/nrn/file/9d4ab20927bc/src/gnu/NegExp.cpp
	negexp = -mean*log(pick())
}

FUNCTION my_poisson(mean) { :adapted from http://www.neuron.yale.edu/hg/neuron/nrn/file/9d4ab20927bc/src/gnu/Poisson.cpp
	LOCAL bound,count,product
	bound = exp(-1.0 * mean)
	count = 0
	product = 1
	while(product >= bound) {
		product = product * pick()
		count = count + 1
	}
	my_poisson = count - 1
}
