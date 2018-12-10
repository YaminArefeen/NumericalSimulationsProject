# NumericalSimulationsProject
Github repository for the MIT Numerical Simulations Class Project written by Zijing Dong, Yamin Arefeen, and Victoria Preston.  

Found in this repository is all code generated for each project milestone check-in meeting in addition to the final deliverable for the course. An outline of the course structure follows. The code in this repository is representative of all project milestones. Fo;der `Simulation_Tumor_Thermal_Therapy` contains the code and resources necessary for running the demonstration presented at the end of the course.

The project pursued by this team was simulating heat diffusion through brain tissue during the process of thermal therapy. The thermal input is treated as a black-box function, and the brain matter is defined by classifications available on a segmented image of an MRI scan. The resolution of the model is directly related to the resolution of the segmentation. 

By applying nonlinar model order reduction techniques and dynamic trapezoidal integration the speed at which a proposed thermal treatment can be assessed in time and in steady-state is vastly sped-up for acceptable accuracy. 


### Project Milestone 1
Initial formulation for the project is proposed.

### Project Milestone 2
Formal project proposal is submitted.

### Project Milestone 3 
Implement an initial integration or solving scheme by utilizing an iterative method discussed in class (Newton Solve, GCR, etc.)

### Project Milestone 4
Introduce or address nonlinearity in the proposed system.

### Project Milestone 5
Consider dynamical inputs and model order reduction.

### Project Milestone 6 (optional)
Consider any other method discussed in class. We elected to implement alternative model order reduction tehcniques and linearization techniques.

## Repository Structure
`Misc_PM_Material` contains all code for project check-in that were supplementary to the main project repository.

`2DHeating_Brain` contains an initial fully-linear version of our project implementation, prior to PM4. `demo_2DHeat_process_Brain.m` is the primary script to run in this folder.

`2DHeating_Brain_nonLinear` contains the final version of the main code used to generate the demosntration material. `exploringModelReduction.m` demonstrates multiple model order reduction techniques (except POD), `exploringInegration.m` demonstrates the dynamic trapezoidal method. Several folders contain information and results from PM5 and PM6 meetings. 

`2DHeating_Brain_POD` contain code to perform proper othogonal decomposition (POD) analysis on the linear system; a strategy that was ultimately was used to generate demonstration material.

`Simulation_Tumor_Thermal_Therapy` contains demonstration material. `demo_Tumor_Heating_final.m` is the relevant script.

