README
Samantha Cohen
CIS 497 DMD Senior Project
Code Version 1

Accomplishments:
	Programming:
	1) All interaction force calculations between normal/coarse/low-resolution particles coded
	2) Primary solver loop for normal/coarse/low-resolution particle physics at each time step created
	3) Segmentation of solver loop for water using SPH (completed in a previous course) and granular solver loop
	4) Have a basic playable simulation of interacting particles
	5) Particles colored based upon their velocity (for debugging purposes)
	6) Added boundary particles but not the corresponding forces they will enact on granular particles

	Research:
	1) Gathered most papers I think I will need for this project
	2) Created preliminary bibliography of said papers
	3) Read through and tried to implement the math in all of the papers on paper and in code

	Miscellaneous:
	1) Created a basic schedule of when I would like different parts of my project to be finished

Issues Encountered:
	Programming:
	1) Though running without error, the simulation does not look like it is working correctly, but I am not sure if this is a math/physics issue relating specifically to particle forces, boundary condition handling, or the height at which the initial collection of particles (in a cube shape) are dropped
	2) Granular material simulations all seem to be based upon Predictive-Corrective Smoothed Particle Hydrodynamics (PCISPH), and I had previously implemented a different version of SPH called Position Based Fluids (PBF) which uses SPH but discards the pressure term in favor of constraint calculations that directly affect particle position rather than velocity, and also uses a fake viscosity term, whereas PCISPH calculates pressure and viscosity and applies them instead, and reconciling these two methods has proven to be a bit of a challenge

	Research:
	1) The papers I have been reading all speak of similar methods but have differences in certain equations, solver loop set up, etc. and it has been difficult to discern which methods work best
	
Goals for Version 2:
	Programming:
	1) I am going to implement boundary particle force calculations and use that to deal with resolving boundaries-the method I used for PBF in my base code does not seem to be working out as well/I cannot tell if it is working so I am going to follow the papers I have been reading and try this new method of collision handling where particles are assigned to the surface of a rigid body and then the corresponding forces they enact on granular particles will be added
	2) I would like to have the normal/coarse/low-resolution particle simulation working pretty much correctly in the next version, so I expect to go through and debug until I achieve this (without accomplishing this much I will not be able to proceed with the rest of my project, so I will have to do this for as long as it takes)