README
Samantha Cohen
CIS 497 DMD Senior Project
Code Version 3

Accomplishments:
	Programming:
	1) Debugged stress and friction calculations such that the simulation now runs without breaking

	Research:
	N/A

	Miscellaneous:
	N/A

Issues Encountered:
	Programming:
	1)The simulation started breaking due to stress and friction/cohesion calculations stemming from the squaring of a kernel value (I am not sure whether or not it should be squared at this point but as of now the value is not squared)

	Research:
	1) I am unsure of which kernel calculations I should be using for the stress and friction calculations
	2) Many equations use the transpose of different kernels and I am not sure exactly how that translates into code
	
Goals for Version 4:
	Programming:
	1) I would like to start debugging discrete particle force calculations to improve the simulation
	2) I would like to change the initial position of the particles such that I can more easily compare their behavior to that of other sand simulations
	3) I would like to start implementing the boundary forces
	2) I would like to have the normal/coarse/low-resolution particle simulation working pretty much correctly as soon as possible, so I expect to go through and debug until I achieve this (without accomplishing this much I will not be able to proceed with the rest of my project, so I will have to do this for as long as it takes)