# Machine-learning-enhanced strategies to generate subtour elimination constraints for the symmetric traveling salesman problem

This repository contains the tested instances and code for the algorithms discussed in the paper:

**"Machine-learning enhanced strategies to generate subtour elimination constraints for the symmetric traveling salesman problem"**  
by *Fatih Burak Akcay, Maxence Delorme, and Ulrich Pferschy*.

---

## Implementation Details

- All algorithms (except ML training) are implemented in **C++**.
- ILP models require the commercial solver **Gurobi** (version 11.0.3 was used).

---

## Repository Structure

Most folders contain the following files:

- `main.h / main.cpp`  
  → Code for calling the implemented methods to solve a given instance  

- `model.h / model.cpp`  
  → Implementation of the main optimization model(s)  

- `helper_functions.h / helper_functions.cpp`  
  → Instance reading, distance computation, and utility functions  

- `cycleGen_helper_functions.h / cycleGen_helper_functions.cpp`  
  → Procedures for generating MST and reduced-cost cycles  

---

## Instances

The folder `INSTANCES` contains `.txt` files for all test instances.

These are divided into four datasets:

- `TSP_LIB`
- `HardTSP_LIB`
- `PORTGEN`
- `PORTCGEN`

---

## Folder Overview

### 1. Branch & Cut Framework Implementations

Contains implementations of the B&C frameworks:

#### 1.1 `CTx_y_OS`
Models:
- CT1-SM-OS, CT1-50-OS, CT1-100-OS  
- CT2-SM-OS, CT2-50-OS, CT2-100-OS  
- CT3-100-OS  

#### 1.2 `CTx_y_ASEF`
Models:
- CT1-SM-ASEF, CT1-50-ASEF, CT1-100-ASEF  
- CT2-SM-ASEF, CT2-50-ASEF, CT2-100-ASEF  
- CT3-100-ASEF  

#### 1.3 `CTx_y_AS`
Models:
- CT1-SM-AS, CT1-50-AS, CT1-100-AS  
- CT2-SM-AS, CT2-50-AS, CT2-100-AS  
- CT3-100-AS  

#### 1.4 `CTx_y_OTF`
Models:
- CT1-SM-OTF, CT1-50-OTF, CT1-100-OTF  
- CT2-SM-OTF, CT2-50-OTF, CT2-100-OTF  
- CT3-100-OTF  

---

### 2. Minimum and Maximum Number of SECs

Contains ILP models:

- `cmin` → Implementation of the **c_min** model  
- `cmax` → Implementation of the **c_max** model  

---

### 3. ML Proof of Concept

Contains code to generate **Table 7** and **Table 8**.

> Note: The size of the bitset must be set before compilation (used value: `2000`).

- `Random` → Code for Table 7  
- `MST_RC` → Code for Table 8  

---

### 4. Cycle Generation Algorithms

Contains code to generate **Table 9** and **Table 10**:

- `P1_P2_strategies`  
  → Code for Tables 9 and 10  
  → Modify parameters in `cycleGen_helper_functions.cpp` under **TESTING CONFIGURATION**

- `RC_MSTdist`  
  → Table 10 (strategy: K=3, MST=L)

- `RC_MSTrc`  
  → Table 10 (strategy: K=3, MST=R)

---

### 5. ML Enhancement

#### 5.1 `ML_training`
- Training of the ML agent (**Python**)

#### 5.2 `ML_Cpp`
- Final model version:
  - Generates cycles  
  - Calls trained ML agent for predictions  
  - Incorporates predictions as **pre-cuts**

---

## Contact

For questions regarding the code, please contact:

f.b.akcay@tilburguniversity.edu
