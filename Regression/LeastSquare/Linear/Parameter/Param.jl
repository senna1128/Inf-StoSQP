
StoSQPSet = Parameter.StoSQP(true,
	                  1e5,                           # Max_Iter
	                  200,                           # Rep
	                  [0,20,40,60],                  # tau
					  1,                             # c_1
	                  [0.501],                       # c_2
					  2,                             # c_3
					  [0,0.4,0.5,0.6],               # RToe
					  [0.1,0.2,0.3],                 # REqu
					  [5,20,40,60],                  # dimension
					  [0.5])                         # constraints

ConstMSet = Parameter.ConstM(
	                  1e5,                           # Sample Size
	                  200)                           # Rep