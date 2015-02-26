README
Samantha Cohen
CIS 497 DMD Senior Project
Code Version 6

Accomplishments:
	Programming:
	1) Added high resolution particles to the simulation at 7 points around/in the low resolution particles
	2) Added high resolution particle math such that though the high resolution particles generally conform to the normal external forces physics that all low-resolution particles are subject to, they have a different velocity calculation for when they have at least one neighbouring low-resolution particle

	Research:
	1) I have found other papers that talk about modeling granular particles as one course particle with smaller particles attached to it, effectively making non-spherical particles, such that the group of 4 or so small particles and one large one move as a single unit, as one low-resolution particle. I am considering adding this to my implementation.

	Miscellaneous:
	N/A

Issues Encountered:
	Programming:
	1)Originally, I set the high-resolution particles to move only corresponding to separate calculations that I had thought were supposed to always be used for them, and this caused odd behaviour such as all of the high resolution particles hovering above the low-resolution particles. I ended up realizing that the high resolution particles should also be using the low resolution particle physics and only use the new separate high resolution calculations when they had at least one or two low resolution neighbours. After fixing this the high resolution particles behaved normally.

	Research:
	1) The paper I am going off of says that I should have a random sampling of particles at 7 points in/around the low-resolution particles. However, I decided to just make those 7 points smaller high resolution particles instead.
	
Goals for Version 7:
	Programming:
	1) I would like to create non-spherical low-resolution particles by attaching smaller particles to the low-resolution particles and having them move as a unit.
	2) I would like to add more high-resolution particles but of varying size to make the simulation more realistic with more particles and more variance.