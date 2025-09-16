#include <stdio.h>
#include <math.h>

const float n= 0.68, E= 1.0, R= 100.0, Is= 1e-15, V0= 0.025; 

// f(U) = E - U - R * Is * (exp(U*(n/V0)) - 1) NEWTON
