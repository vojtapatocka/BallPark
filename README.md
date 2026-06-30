# BallPark
Numerical model for building scoria cones by ballistic deposition of particles

To run the code, download the folder and then type
"./gun.bat here"
on a Linux machine (gfortran is the default compiler). Use 
"./gun.bat myrun"
if you want the simulation to run on background. 

Check OMP_NUM_THREADS in gun.bat to set the number of OpenMP threads.
The code uses param.in on input, check param.info for the description of most numerical parameters.
param_Heaney.in, param_NWCaloris.in, and param_NECaloris.in are the best-fit models from 
The Solar System’s largest potential scoria cone is on Mercury

When the simulation is finished, type
"./postpro.bat cone"
to plot the obtained solution (cone morphology)
