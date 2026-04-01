# ML-TSP
This repository contains the tested instances and code for algorithms discussed in the paper "Machine-learning enhanced strategies to generate subtour elimination constraints for the symmetric traveling salesman problem"  by Fatih Burak Akcay, Maxence Delorme and Ulrich Pferschy.

All algorithms (except ML training) are coded in C++.
ILP models require the commercial solver Gurobi (we used version 11.0.3).
 
Almost all folders have the following contents:
-main.h/cpp   						| code for calling our methods to solve a given instance
-model.h/cpp						| code implementing the main optimization model(s)
-helper_functions.h/cpp				| code containing instance reading, distance computations, and miscellaneous supportive functions
-cycleGen_helper_functions.h/cpp 	| code for generating mst and reduced-cost cycles procedures

Moreover, folder "INSTANCES" contains a txt-file for each of our test instances. These are spread over 4 subfolders, each corresponding to a different data set: TSP_LIB, HardTSP_LIB, PORTGEN, PORTCGEN.

For questions on the code, please send an e-mail to f.b.akcay@tilburguniversity.edu

Folders:
1. B&C framework implementations: contains the implementations of the B&C frameworks
	1.1. CTx_y_OS: contains the models CT1-SM-OS, CT1-50-OS, CT1-100-OS, CT2-SM-OS, CT2-50-OS, CT2-100-OS and CT3-100-OS.
	1.2. CTx_y_ASEF: contains the models CT1-SM-ASEF, CT1-50-ASEF, CT1-100-ASEF, CT2-SM-ASEF, CT2-50-ASEF, CT2-100-ASEF and CT3-100-ASEF.
	1.3. CTx_y_AS: contains the models CT1-SM-AS, CT1-50-AS, CT1-100-AS, CT2-SM-AS, CT2-50-AS, CT2-100-AS and CT3-100-AS.
	1.4. CTx_y_OTF: contains the models CT1-SM-OTF, CT1-50-OTF, CT1-100-OTF, CT2-SM-OTF, CT2-50-OTF, CT2-100-OTF and CT3-100-OTF.

2. Minimum and maximum number of SECs: contains the ILP models of the c_min and c_max
	2.1. cmin: contains the c_min model
	2.2. cmax: contains the c_max model

3. ML proof of concept: contains the codes to generate Table 7 and Table 8. (Note the size of the bitset has to be set before compilation. We used 2000.)
	3.1. Random: contains the code to generate Table 7.
	3.2. MST_RC: contains the code to generate Table 8. 
	
4. Cycle generation algorithms: contains the codes to generate Table 9 and 10.
	4.1. P1_P2_strategies: contains the code for Table 9 and 10. To use, change the values in the cycleGen_helper_functions.cpp under TESTING CONFIGURATION.
	4.2. RC_MSTdist: contains the code for Table 10 on strategy K=3, MST=L
	4.3. RC_MSTrc: contains the code for Table 10 on strategy K=3, MST=R
