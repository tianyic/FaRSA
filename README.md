<img src=https://github.com/tianyic/FaRSA/blob/master/farsa.jpg?raw=true width=135/> Fast Reduced Space Algorithm
=====
---
[Installation](https://tianyic.github.io/)|
[Documentation](https://tianyic.github.io/)|
[License](https://tianyic.github.io/)|

FaRSA is an optimized library for solving l1-regularized convex optimization problem designed to be **efficient**, **flexible** and **user-friendly**.  l1-regularized convex optimization problem plays an important role in machine learning applications. An optimization framework [FaRSA](http://www.optimization-online.org/DB_FILE/2016/02/5331.pdf) has been proposed to efficiently solve such problem.  FaRSA implements highly optimized solutions for several widely used loss functions, like logistic loss, and supports generic loss function optimization.  

---

##Installation Guide:

#### 1. Download source file:

- Download zip file from our [homepage](https://tianyic.github.io/)

- Clone it from our github repository:
	
		git clone https://github.com/tianyic/FaRSA
		
- For Unix users, you can also get our source files via: 

		wget someurl.FaRSA.zip
	
After obtaining our source file, please uncompress zip file. 

#### 2. Build on OSX 

On OSX, make sure gcc has been installed. If not, to install gcc compiler on Mac OS X, you need to download and install "Command Line Tools for Xcode", which is available in Apple’s developer page.

Then jump to uncompressed FaRSA directory, run the below command in terminal to compile:

	make

If compiling failed due to permission issue, please compile by 

	sudo make

To test whether installation is successful, run the below command:

	farsa -t test

#### 3. Build on Ubuntu or other Unix Machines

Make sure gcc complier has been installed on your machine, if not, please run the below commands to install necessary compliers:

	sudo apt-get update
	sudo apt-get install build-essential

FaRSA directory, run the below command in terminal to compile:

	make

If compiling failed due to permission issue, please compile by 

	sudo make

To test whether installation is successful, run the below command:

	farsa -t test
	
---

## Documentation:

### 1. Quick Start 

In general, FaRSA contains special optimized routine and a generic rountine. Then 
special optimized routine provides highly efficient solutions for solving the below popular optimization problem:

- logistic loss 

- ...


Generic routine can also support solving personalized loss functions by a generic routine.

In order to successfully call FaRSA to solve aimed problem, users are expected to finish the below steps:

- Configure FaRSA profile

- Prepare readable dataset if necessary

- If aim problems do not belong to FaRSA's existed problem solutions, for C users, modifying **client.c** based on aimed problem.

Once the above steps are complete, run the below command in terminal to call FaRSA for solving aimed problem:

	farsa -p path_to_param_profile

Results would be displayed in terminal.



### 2. Configure FaRSA Profile

To employ FaRSA on aimed problem correctly, you need to define a FaRSA profile. Two demos of profiles have been provided in the FaRSA source folder, named FaRSA1.profile, FaRSA2.profile. 

There are two type parameters in general: basic and advanced. 

- Basic parameters are necessary for successfully running the FaRSA. 

- Advanced parameters are designed for experienced users to control parameters in optimization algorithm.

The format in profile is 

	parameter1 : value1
	parameter2 : value2
	parameter3 : value3
	...


#### 2.1 Basic Parameters 

- **objective function type** : Objective function type for choosing whether use special optimized routine  or generic rountine. Feasible values are 1 or 3.

	- 1 : call logistic loss rountine
	
	- 3 : call generic rountine

- **lambda coefficient**: Designed for special optimized routine. **Feasible range: [ 0, infinity )**, the **lambda** in aimed optimization problem would be **lambda coefficient / number of samples**, where the number of sample is data file's sample set size.

- **lambda**: Set lambda in aimed optimization problem. Complementary for lambda coefficient. **Feasible range: [ 0, infinity )**. If lambda is nonnegative, then FaRSA would ignore lambda coefficient, otherwise FaRSA would use lambda coefficient to calculate lambda.

- **data file**: Necessary for special optimized routine, not required for generic routine. Provide absolute path to access data file. 

- **data format**: Necessary for special optimized routine, not required for generic routine. Declare data file format. Current feasible values: libsvm. 


#### 2.2 Advanced Parameters

- **maximum iteration**: Maximum iteration of FaRSA optimization algorithm. Feasible range: positive integer set. Default: 1000

- **Gamma**: Coefficient for switching between phi step or beta step. Feasible range: ( 0, 1 ]. Default: 1

- **eta_r**: Coefficient to determine target residual for conjugate gradient method. Feasible range: ( 0, 1 ]. Default: 0.1 

- **eta**: Coefficient in sufficient descrease. Feasible values: ( 0, 1 ]. Default: 0.01.

- **xi**: Shrinkage coefficient for backtracking line search. Feasible values: ( 0, 1 ]. Default: 0.5

- **tolerance**: Tolerance for termination. Feasible range: nonnegative float set. Suggest to set it from [ 1E-6, 1E-1 ]. Default: 1E-6.

- **betaFrac**: Beta fraction to set the ratio of variables in beta step's calculation. Feasible range: ( 0, 1 ]. Default: 1.0

- **phiFrac**: Phi fraction to set the ratio of variables in phi step's calculation. Feasible range: ( 0, 1 ]. Default: 1.0

- **max_CG_iter**: Maximum conjugate gradient method iteration. Feasible range: positive integer set. Default: 1000

- **TRradiusPhi**: Trust-region-like adaptive strategy radius for preventing the size of direction in reduced space solver from too large. Feasible range: positive float set. Default: 1000

- **TRradiusBeta**: Trust-region-like adaptive strategy radius for controlling the size of descent direction in beta step. Feasible range: 
positive float set. Default: 0.1

- **maxback**: Maximum allowed backtracking line search time. Feasible range: positive integer set. Default: 100

- **fracViol**: Maximum allowed violating variable ratio for reduced space solver. Feasible range: [ 0, 1 ]. Default: 0.1. Real maximum allowed violating variable also depends on the next parameter: **nVmaxAllow**.

- **nVmaxAllow**: The number of maximum allowed violating variable for reduced space solver. Feasible range: positive integer. Default: 1000. 

### 3 Data Set Prepration

Users are expected to offer data file in specific data format if called FaRSA routine needs. Data set format that current FaRSA supports is libsvm format, which is suitable for sparse data set.  

**libsvm format** let each row represent each sample. The first entry of each row is the label of corresponding sample. Features with value as zero are ignored. Features with nonzero values are aligned as feature index, colon delimiter, feature value, as below:

	sample_1_label feature_1_idx:feature_1_value feature_2_idx:feature_2_value feature_5_idx:feature_5_value
	sample_2_label feature_2_idx:feature_2_value feature_3_idx:feature_3_value feature_4_idx:feature_4_value
	sample_3_label feature_1_idx:feature_1_value feature_3_idx:feature_3_value feature_5_idx:feature_5_value
	sample_4_label feature_1_idx:feature_1_value feature_3_idx:feature_3_value feature_5_idx:feature_5_value
	...


### 4 Construct Personalized Problem

FaRSA's generic routine can support any l1-regularized convex optimization problem. This attractive feature allows users to design their own objective function, instead of being limited by the problems in special optimized routine. To achieve calling generic routine, users are required to fill in **structure Input_FaRSA** to let FaRSA understand what the personalized objective function is.

```
struct Input_FaRSA{
	int n;
	double ( *func )( double * );
	double *( *grad_f )( double * );
	double *( *hessVec )( double * );
};

```

In general, there are three functions to be implemented:

- A function to calculate objective function value given current x

- A function to calculate gradient of f that is objecitve function without l1-regularized term given current x.

- A function to calculate reduced hessian vector product given a vector. 




#### 4.1 Function To Caculate Objecitve Function Value Given Current x

The template of this function is the following:

```
double personalized_func( double *x ){

	double func = 0.0;
	
	/*  Users fill in based on aimed problem */
	
	return func;
}

```



#### 4.2 Function To Caculate Gradient Of f Given Current x

The template of this function is as below:

```
double *personalized_grad_f( double *x ){

	double *grad_f;
	
	/* Users fill in based on aimed problem */
	
	return grad_f;
}

```
**Note:** f is the objective function without l1-regularization term. 


#### 4.3 Function To Calculate Reduced Hessian Vector Product Given A Vector.


The template is here:

```
double *personalized_hessVec( double *v ){

	double *hv;
	
	/*  
		Users fill in based on aimed problem
		Return reduced_hessian_matrix * v
		
	 */
	
	return hv;
}
```

To fill in the above template, three global variables maybe useful and can be directly used without any declaration: 

- **x**: A double pointer points to a double array to saves current **x**.

- **S**: A integer pointer points to a integer array that saves free variables' indexes at current **x**.

- **nS**: A integer save the number of free vairables in **S**.


 --- 
 
## License 

© Contributors, 2016.

FaRSA is free to academic use. 

Contributors:

- Tianyi Chen, PhD Candidate, Department of Applied Mathematics and Statistic, Johns Hopkins University. Email: <tchen59@jhu.edu>

- Daniel P. Robinson, Assistant Professor, Department of Applied Mathematics and Statistic, Johns Hopkins University. Email: <daniel.p.robinson@gmail.com>

- Frank E. Curtis, Associate Professor, Department of Industrial and Systems Engineering, Lehigh University. Email: <frank.e.curtis@gmail.com>


 
 
 
 




























