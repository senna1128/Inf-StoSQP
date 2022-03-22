
StoSQPSet = Parameter.StoSQP(true,
	                  1e5,                           # Max_Iter
	                  1,                             # Rep
	                  50,                            # tau
	                  [0.5,0.6,0.7,0.8,0.9,1],       # Decay stepsize
	                  [1e-8, 1e-4, 1e-2, 1e-1, 1])   # Sigma
