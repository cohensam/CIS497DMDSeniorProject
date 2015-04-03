README
Samantha Cohen
CIS 497 DMD Senior Project
Code Version 9

Accomplishments:
	Programming:
	1) I regressed to Version 5 of my code before I added the high-resolution particles since there was a lot of code involved in that and I needed a simpler version with only the initial 1000 particles so I could debug the yield and cohesion conditions. 
	2) Was able to add rigidity as an attribute of the particle objects and create clusters of rigid particles, and averaged each particle's velocity contribution to the whole cluster's movement and then applied the averaged velocity to all particles.
	3) Simulation now exhibits more realistic clumping.

	Research:
	1) Was able to actually understand what was supposed to happen if a yield condition was met (the particles are designated as rigid).
	2) Realized that when particles are rigid, they are added to "clusters", or collections of rigid particles that are connected to other rigid particles through rigid particle neighbors. The clusters then move as units with the same velocity.

	Miscellaneous:
	N/A

Issues Encountered:
	Programming:
	N/A

	Research:
	N/A
	
Goals for Version 10:
	Programming:
	1) I would like to create a MEL script to load the granular particle positions into Maya, assign their positions to some geometry, and render out a test of the sumulation.
	2) I would like to add the high resolution particle math back into this version of the code from the previous version, Version 7.