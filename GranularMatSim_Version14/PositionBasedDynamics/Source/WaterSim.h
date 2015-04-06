#ifndef WATER_SIM_H_INCLUDED
#define WATER_SIM_H_INCLUDED


#include "openGL_headers.h"
#include "math_headers.h"
#include "particlelist.h"
#include "vec.h"
#include "triangle.h"
#include "constraint.h"
#include "scene.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
using namespace std;
#include "sphere.h"
#include "array1.h"
#include "array2.h"
#include "array3.h"
#include "vec.h"
#include <math.h>
//#include "marching_cubes.h"

class WaterSim
{
public:
	//Water-Granular Interaction
	int granularParticleNum;
	void interpolateBridge(glm::vec3 R, float h, int i, vector<int> nbrs);
	void disperseWaterWetness(glm::vec3 R, float h, int i, vector<int> nbrs);
	void disperseGranularWetness(glm::vec3 R, float h, int i, vector<int> nbrs);
	vector<int> dispersedWaterParticles;
	float Wthreshold;
	float Wmax;

	//Granular Simulation
	void assignRigidClusters();
	void assignRigidBodyVelocity(float dt);
	vector<vector<int>> rigidClusters;
	vector<glm::vec3> rigidClustersForces;
	vector<glm::vec3> gravityForces;
	float rhoMax;
	float angleRepose;
	void interpolateFriCoh(glm::vec3 R, float h, int i, vector<int> nbrs);
	void calculateCorrectivePressure(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void calculateDiscreteParticleForces(int i, vector<int> nbrs);
	void calculateStrainRate(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void calculateCorrectiveStress(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void testYieldandCohesion(int i, float dt);

    WaterSim();
    WaterSim(unsigned int n);
    virtual ~WaterSim();
	vector<vec3> objMeshTriangles;
	vector<int> objMeshFaces;
	vector<float> objMeshVertices;

	void initialize();

	//Load .obj file
	void objLoader(string filename);

	//Number of particles
	int particleNum;
	float particleRad; 
	float rho_0;
	//List of particles
	ParticleList particleList;
	//Grid of particles (hash table)
	int gridX; 
	int gridY; 
	int gridZ; 
	std::hash<glm::vec3> hash; 
	unordered_map<int,std::vector<int>> table; 
	std::vector<int> keys;
	std::vector<std::vector<int>> vals;
	std::vector<float> phiCell;

	//Main algorithm loop
	void update(const Scene* const scene, float dt);
	int solverIterations;

	//Visualize particles
	void particleToSphere();
	vector<Sphere*> dots;
	void updateSpheres();

	//List of triangles making up the obj mesh
	vector<triangle*> triangles;
	float tri_xmax, tri_xmin, tri_ymax, tri_ymin, tri_zmax, tri_zmin; 

	void resolve_constraints(const Scene* const scene); 
	//void collision_detection(const Scene* const scene);
	void collision_detection(int id, const Scene* scene);
	//void resolve_collisions(); 
	bool box_intersection(const glm::vec3& p1, const glm::vec3& p2, float threshold, glm::vec3& intersect, glm::vec3& normal) const; 
	void compute_predicted_position(float dt); 

	void obj_collision(int id); 
	bool triangle_collision(/*const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& v,*/ int tri, int id); 
	void compute_mesh_bounds();
	bool WaterSim::check_obj_bound(const glm::vec3 p);

	//List of collision constraints-currently used to deal with particle-floor collisions
	 std::vector<CollisionConstraint> m_constraints_ext;
	 float m_thick;

	 //Neighbor search 
	 void populateGrid(); 
	 std::vector<int> neighborSearch(int id);
	 bool isValidCell(int i, int j, int k); 
	 int hash_function(int i, int j, int k); 

	 float nbrRadius; 

	 //Constraint Projection
	 float calculateLambda(int id); 
	 glm::vec3 calculateDeltaP(int id); 

	 //Kernel calculations
	 float kernelGeneral(float r, float h);
	 float kernelPressure(float r, float h);
	 float kernelViscosity(float r, float h);
	 float kernelCohesion(float r, float h);

	 //Kernel derivs
	 glm::vec3 kernelPressureDeriv(glm::vec3 R, float h);
	 glm::vec3 kernelVorticityDeriv(glm::vec3 R, float h);
	 glm::vec3 kernelGeneralDeriv(glm::vec3 R, float h);

	 //Kernel Second Deriv
	 float kernelViscositySecondDeriv(glm::vec3 R, float h);
	 float kernelGeneralSecondDeriv(glm::vec3 R, float h);

	 //Anisotropic kernel calculations for phi level set calculations
	 float kernelAnisotropic(glm::vec3 R, float h);
	 float kernelAnisotropicSpline(glm::vec3 R, float h);
	 float interpolatePhi(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void marchingCube(const VBO& vbos);

	 //Modified from https://raw.githubusercontent.com/christopherbatty/Fluid3D/master/main.cpp
	 void export_particles(string path, int frame, const std::vector<glm::vec3>& particles, float radius);

	 //SPH interpolation
	 float interpolateDensity(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void interpolatePressure(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void applyXSPH(glm::vec3 R, float h, int i, vector<int> nbrs);

	 //Box Bounds
	 float xmin, xmax;
	 float ymin, ymax;
	 float zmin, zmax;

	 // Grid Bounds
	 float gridXmin, gridXmax; 
	 float gridYmin, gridYmax; 
	 float gridZmin, gridZmax;
	 vector<vector<vector<float>>> grid;//[133][133][133];
	 vector<vector<vector<float>>> phis;//[133][133][133];

};

#endif