
StoSQPSet = Parameter.StoSQP(true,
	                  1e5,                           # Max_Iter
	                  200,                           # Rep
	                  20,                            # tau
					  1,                             # c_1
	                  [0.501],                       # c_2
					  2,                             # c_3
	                  [1e-4,1e-2,1e-1,1])       # Sigma
