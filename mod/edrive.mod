TITLE Drive

NEURON { 
  SUFFIX Edrive 
  NONSPECIFIC_CURRENT edrive,edrivenoise 
  RANGE drive,drivenoise 
} 
PARAMETER { 
drive = 0		(milliamp/cm2)
drivenoise = 0		(milliamp/cm2)
} 
ASSIGNED { 
edrive (milliamp/cm2) 
edrivenoise (milliamp/cm2) 
} 
BREAKPOINT { 
edrive = -drive
edrivenoise = -drivenoise
 } 

