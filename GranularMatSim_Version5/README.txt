README
Samantha Cohen
CIS 497 DMD Senior Project
Code Version 5

Accomplishments:
	Programming:
	1) Added yield and cohesion conditions, so now the simulation exhibits clumping such that the particles are mostly confined to a large mass instead of covering the floor and continuously spreading out

	Research:
	1) Was finally able to make some sense of the yield and cohesion constraints

	Miscellaneous:
	N/A

Issues Encountered:
	Programming:
	1)The particles kept spreading all around the floor and never stuck together, whereas a real pile of sand would ideally form a cone and have the particles stop spreading due to friction and cohesion
	2) I am unsure as to whether or not using my old implementation of pressure calculations is actually ok to do, because friction/cohesion force calculations rely on stress calculations, which in turn rely on pressure calculations. However in my previous SPH water simulation I bypassed the pressure calculation by solving certain constraints and directly updating the position of the particles without using pressure interpolation (though I did use the SPH pressure kernel deriv). Right now I am using both sets of calculations, and the sim seems to look ok so I am going to leave it for now but I do feel conflicted about it

NOTE: The corrective pressure calculations are updating pressure, not pressure force, whereas my method interpolating the pressure actually updates the pressure force, but I never use the pressure force or pressure for my final pressure bypass calculation, however all use the same kernels which I am confident are correct, so perhaps there is actually nothing to worry about

	Research:
	N/A
	
Goals for Version 6:
	Programming:
	1) I would like to test different values for the yield and cohesion constraints to try and get the simulation to behave more similarly to the video I found that directly corresponds to one of the papers I am working off of (by which I mean I want the simulation to look more cone-like)
	2) I would like to start working on the high resolution particle calculations