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
	//Granular Simulation
	float tolP;
	float tolS;
	float rhoMax;
	float angleRepose;
	void interpolateFriCoh(glm::vec3 R, float h, int i, vector<int> nbrs);
	void calculateCorrectivePressure(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void calculateDiscreteParticleForces(int i, vector<int> nbrs);
	void calculateStrainRate(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	glm::vec3 interpolateU(glm::vec3 R, float h, int i, vector<int> nbrs);
	void calculateCorrectiveStress(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void testYieldandCohesion(int i);
	float calculateFerr();
	float calculateStressTensor(float dt);
	int boundaryParticleNum;
	int xy;
	int xz;
	int yz;
	void interpolateBoundaryDensity(glm::vec3 R, float h, int i, vector<int> nbrs);
	float calculateLambda(glm::vec3 R, float h, int i, vector<int> nbrs);
	void interpolateBoundaryPressure(glm::vec3 R, float h, int i, vector<int> nbrs);
	glm::vec3 interpolateBoundaryU(glm::vec3 R, float h, int i, vector<int> nbrs);
	glm::vec3 interpolateBoundaryFriction(glm::vec3 R, float h, int i, int j);
	glm::vec3 interpolateBoundaryViscosity(glm::vec3 R, float h, int i, int j, float dt);
	float calculatePi(float h, int i, int b, float dt);
	void calculateBoundaryForce(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	void interpolateCohesion(glm::vec3 R, float h, int i, vector<int> nbrs);

	//High Resolution Calculations
	float hrCalcWeight(float r);
	glm::vec3 hrCalcAvgVelocity(glm::vec3 R, float h, int i, vector<int> nbrs);
	void hrCalcVelocity(glm::vec3 R, float h, int i, vector<int> nbrs, float dt);
	float hrCalcAlpha(glm::vec3 R, float h, int i, vector<int> nbrs);
	void hrCalcPosition(glm::vec3 R, float h, int i, vector<int> nbrs, float dt); 
	std::vector<int> hrNeighborSearch(int i);



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
	ParticleList waterParticleList;
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
	void updateWater(const Scene* const scene, float dt);
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
	 //Modified from https://raw.githubusercontent.com/christopherbatty/Fluid3D/master/fluidsim.cpp
	 void compute_phi();
	 Array3f nodal_solid_phi;  
	 Array3f liquid_phi;
	 int ni,nj,nk;
	 float phi_dx;

	 //SPH interpolation
	 float interpolateDensity(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void interpolatePressure(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs);
	 void interpolateVorticity(glm::vec3 R, float h, int i, vector<int> nbrs);
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

	 void draw(const VBO& vbos); 

	 //Rendering Visualization with Level Sets
	 vector<int> surfaceSampling(vector<int> levelSet, float r, float t, float e);

    /*void initialize(unsigned int dim_x, unsigned int dim_z, const glm::vec3& cloth_min, const glm::vec3& cloth_max);
    void update(const Scene* const scene, float dt);
    void draw(const VBO& vbos);
    void flip_draw_mode()
    {
        m_draw_wire = !m_draw_wire;
    }
protected:
    struct Edge
    {
        unsigned int m_v1, m_v2;
        unsigned int m_tri1, m_tri2;
    };
protected:
    unsigned int m_dimx, m_dimz;
    float m_thick;
    unsigned int m_solver_iterations;
    // vertices and estimated position.
    ParticleList m_vertices;
    // internal and external constraints.
    std::vector<Constraint*> m_constraints_int;
    std::vector<CollisionConstraint> m_constraints_ext;
    std::vector<SelfCollisionConstraint> m_self_collision;
    //std::vector<Constraint*> m_constraints_ext;
    // for visualize the cloth.
    bool m_draw_wire;
    std::vector<glm::vec3> m_normals;
    std::vector<glm::vec3> m_colors;
    std::vector<unsigned int> m_triangle_list;
    // for generating constraints.
    std::vector<Edge> m_edge_list;
private:
    // generate edge list from the geometry representation.
    void generate_edge_list();
    // generate all the internal constraints based on the edge list. 
    void generate_internal_constraints();
    // update the normal per frame for visualization.
    void compute_normal();
    // apply external force to the system.
    void apply_external_force(const glm::vec3& force, float dt);
    // damp velocity for all vertices.
    void damp_velocity(float k_damp);
    // compute predicted position based on velocity.
    void compute_predicted_position(float dt);
    // collision detection, generating external constraints. need to generate constraints per frame.
    void collision_detection(const Scene* const scene);
    // self collision.
    void self_collision_detection();
    // resolve all the constraints, both internal and external.
    void resolve_constriants();
    // update the position and velocity.
    void integration(float dt);
    void update_velocity(float friction, float restitution);
    void clean_collision_constraints();*/

	 // These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

// a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
unsigned int a2iEdgeConnection[12][2];

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
float a2fEdgeDirection[12][3];

// aiCubeEdgeFlags lists flags for all possible cases
// each flag has 12 bits, corresponding to 12 edges. 0 for no intersection, 1 for intersection.
int aiCubeEdgeFlags[256];

// a2iTriangleConnectionTable lists for each case, how many triangles are there and which edges they resides on.
// -1 for no triangle any more, 0 ~ 11 indicating the edge id.
// for any case, up to 5 triangles are possible, and we need another value in the end to indicate the table item is done.
// so each item has 5 * 3 + 1 = 16 elements.
int a2iTriangleConnectionTable[256][16];

};

#endif