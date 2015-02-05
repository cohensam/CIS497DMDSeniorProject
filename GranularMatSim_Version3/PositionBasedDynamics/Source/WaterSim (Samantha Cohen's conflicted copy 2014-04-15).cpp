#include "WaterSim.h"

#define PI 3.1415926535897

WaterSim::WaterSim() {
	particleNum = 1000;
	particleRad = 0.1; 
	particleList.resize(particleNum);
	m_thick = 0.15f;
	solverIterations = 100;
	length = 200;
	width = 150;
	height = 100;
	nbrRadius = 0.5;
}

WaterSim::~WaterSim() {
	m_constraints_ext.clear();
}

void WaterSim::initialize() {
    for(int i = 0; i < particleNum; ++i)
	{
		//Set original velocity of each particle
        particleList.vel(i) = glm::vec3(0.0f);
		//Set orignal inv mass of each particle
		particleList.set_inv_mass(i, 1.0f);
    }

	//Set original position of each particle to be a cube of 20x20x20
	int pos = 0;
	for (float i = 0; i < 2; i+=(2*particleRad)) 
	{
		for (float j = 0; j < 2; j+=(2*particleRad))
		{
			for (float k = 0; k < 2; k+=(2*particleRad))
			{
				particleList.pos(pos) = glm::vec3(i,j+0.2f,k);
				pos++;
			}
		}
	}

	//particleToSphere();

}

void WaterSim::compute_predicted_position(float dt)
{
    for(unsigned int i = 0; i < particleNum; ++i)
    {
		// TODO: compute predicted position for all vertices.
        // this is just an example line and assign a initial value so that the program won't be too slow.
        particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);
    }
}

void WaterSim::collision_detection(const Scene* const scene)
{
	// detect collision between water and the scene, and generate external constraints.
	m_constraints_ext.clear();
    unsigned int i; 
    glm::vec3 x, p, q, n;
    for(i = 0; i < particleNum; ++i)
    {
        x = particleList.pos(i);
        p = particleList.predicted_pos(i);
        if(scene->line_intersection(x, p, m_thick, q, n))
        {
            CollisionConstraint c(&particleList, i, q, n);
            m_constraints_ext.push_back(c);
        }
    }
}

void WaterSim::resolve_collisions()
{
	// resolve all collisions
    bool all_solved = true;
    bool reverse = false;
    int size = m_constraints_ext.size();
    for(int i = reverse ? (size - 1) : 0; (i < size) && (i >= 0); reverse ? --i : ++i)
    {
        all_solved &= m_constraints_ext[i].project_constraint();
    }

    reverse = !reverse;
    particleList.unlock_pos_all();
}

void WaterSim::objLoader(string filename) {	
	triangles.clear();
	ifstream inFile(filename, ifstream::in);//"gourd.obj", ifstream::in);
	string line;
	char *a = NULL;
	while (inFile.good()) {
		getline(inFile, line);
		if (line.size() == 0) {
			continue;
		}
		char* tokens = strtok_s(&line[0], " ", &a);
		//Storing vertices as floats
		if (tokens != nullptr && tokens[0] == 'v') {
			tokens = strtok_s(NULL, " ", &a);
			while (tokens != NULL) {
				objMeshVertices.push_back((float)atof(tokens));
				tokens = strtok_s(NULL, " ", &a);
			}
		}
		//Storing faces
		if (tokens != nullptr && tokens[0] == 'f') {
			tokens = strtok_s(NULL, " ", &a);
			while (tokens != NULL) {
				objMeshFaces.push_back((int)atoi(tokens)-1);
				tokens = strtok_s(NULL, " ", &a);
			}
		}
	}
	vector<glm::vec3> ver;
	for (int i = 0; i < objMeshVertices.size(); i+=3) {
		ver.push_back(glm::vec3(objMeshVertices[i],objMeshVertices[i+1], objMeshVertices[i+2]));
	}
	for (int i = 0; i < objMeshFaces.size(); i+=3) {
		triangle* t = new triangle();
		t->color = glm::vec3(1.0,0.0,0.0);
		t->p0 = glm::vec3(ver[objMeshFaces[i]].x, ver[objMeshFaces[i]].y, ver[objMeshFaces[i]].z);
		t->p1 = glm::vec3(ver[objMeshFaces[i+1]].x, ver[objMeshFaces[i+1]].y, ver[objMeshFaces[i+1]].z);
		t->p2 = glm::vec3(ver[objMeshFaces[i+2]].x, ver[objMeshFaces[i+2]].y, ver[objMeshFaces[i+2]].z);
		t->norm = t->normal();
		triangles.push_back(t);
	}
}

void WaterSim::particleToSphere() {
	dots.clear(); 
	Sphere* s;
	for (int i = 0; i < particleList.size(); i++) {
		s = new Sphere(particleList.pos(i),particleRad);
		dots.push_back(s);
	}
}

/*void WaterSim::populateGrid() {
	for (int i = 0; i < particleNum; i++) 
	{
		int ix = particleList.grid(i).x;  
		int iy = particleList.grid(i).y; 
		int iz = particleList.grid(i).z; 
	}
}
	 
std::vector<int> WaterSim::getNeighbors(const int id) {
	vector<int> nbrs; 
	for (int i = -1; i <= 1; i++)
    {
		for(int j = -1; j <= 1; j++)
        {
			for (int k = -1; k <= 1; k++) 
			{

			}
        }
    }
    return nbrs;
}*/

std::vector<int> WaterSim::neighborSearch(const int id) {
	vector<int> nbrs;
	glm::vec3 particlePos = particleList.predicted_pos(id); 
	for (int i = 0; i < particleNum; i++) 
	{
		int nbrid = i; 
		//if (nbrid != id)
		//{
			glm::vec3 nbrPos = particleList.predicted_pos(i); 
			float r = glm::length(nbrPos - particlePos); 
			if (r < nbrRadius) 
			{
				nbrs.push_back(nbrid);
			}
		//}
	}
	return nbrs; 
}

float WaterSim::calculateLambda(const int id) {
	// compute rho_i 
	std::vector<int> Ni = particleList.nbr(id); 
	float rho_0 = 1000; 
	float rho_i; 
	float Wj; 
	std::vector<float> W; 

	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		float dist = glm::length(particleList.predicted_pos(nbrid)-particleList.predicted_pos(id)); 
		Wj = kernelGeneral(dist,nbrRadius); 
		W.push_back(Wj); 
	}
 
	rho_i = 0.0f; 
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		float mj = 1.0f/particleList.inv_mass(nbrid); 
		rho_i += (mj * W[j]); 
	}

	// calculate constraint Ci and lambda_i 
	float epsilon = 0.5f; 
	float Ci = rho_i/rho_0 - 1; 
	float lambda_i = 0.0f; 
	std::vector<float> delW; 
	float delWj; 
	// spike kernel for gradient calculation 
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		float dist = glm::length(particleList.predicted_pos(nbrid)-particleList.predicted_pos(id)); 
		delWj = kernelGeneral(dist,nbrRadius); 
		delW.push_back(delWj); 
	}
	float sumCiGrad = 0.0f; 
	float nbrterm = 0.0f; 
	float particleterm = 0.0f; 
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		particleterm += delW[j]; 
		if (nbrid != id) 
		{
			nbrterm += pow(abs(-1.0f/rho_0*delW[j]),2); 
		}
	}
	particleterm = particleterm / rho_0; 
	sumCiGrad = nbrterm + pow(abs(particleterm),2); 
	lambda_i = -Ci/(sumCiGrad+epsilon); 
	return lambda_i; 
}

float WaterSim::calculateDeltaP(const int id) {
	//
	std::vector<int> Ni = particleList.nbr(id); 
	float rho_0 = 1000; 
	float Wj, delWj, Wq, delq; 
	std::vector<float> W; 
	std::vector<float> delW; 
	std::vector<float> pij; 
	std::vector<float> scorr;  
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		float dist = glm::length(particleList.predicted_pos(nbrid)-particleList.predicted_pos(id)); 
		Wj = kernelGeneral(dist,nbrRadius); 
		W.push_back(Wj); 
	} 

	// scorr calculation for tensile stability 
	float k = 0.1f; 
	int n = 4; 
	delq = 0.1f*nbrRadius; 
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		Wq = kernelGeneral(delq,nbrRadius); 
		float val = -k*pow(W[j]/Wq,n); 
	}

	// spike kernel for gradient calculation 
	float deltaP_i = 0.0f;
	for (int j = 0; j < Ni.size(); j++) 
	{
		int nbrid = Ni[j]; 
		deltaP_i += ((particleList.lambda(id)+particleList.lambda(nbrid)+scorr[j])*delW[j]); 
	}
	deltaP_i = deltaP_i/rho_0; 
	return deltaP_i; 
}

float WaterSim::kernelGeneral(float r, float h) {
	if (0 <= r && r <= h) {
		return 315/(64*PI*pow(h,9)) * pow((pow(h,2) - pow(r,2)),3);
	} else {
		return 0;
	}
}

float WaterSim::kernelPressure(float r, float h) {
	if (0 <= r && r <= h) {
		return 15/(PI*pow(h,6)) * pow((h-r),3);
	} else {
		return 0;
	}
}

float WaterSim::kernelViscosity(float r, float h) {
	if (0 <= r && r <= h) {
		return 15/(2*PI*pow(h,3)) * (-(pow(r,3)/(2*pow(h,3))) + (pow(r,2)/(pow(h,2))) + (h/(2*r)) - 1);
	} else {
		return 0;
	}
}

float WaterSim::interpolateDensity(glm::vec3 R, float h, int i, vector<int> nbrs) {	
	float As = 0.0f;
	glm::vec3 temp;
	float norm = 0;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		norm = sqrt(pow(temp.x,2)+pow(temp.y,2)+pow(temp.z,2));
		As += particleList.inv_mass(nbrs[j])*kernelGeneral(norm, h);
	}
	return As; 
}

void WaterSim::interpolatePressure(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		As += particleList.inv_mass(nbrs[j])*((particleList.press(i)+particleList.press(nbrs[j]))/(2.0f*particleList.rho(nbrs[j]))) * kernelPressureDeriv(temp,h);
	}
		particleList.press(i) = -As;
}

void WaterSim::interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		As += particleList.inv_mass(nbrs[j])*((particleList.vel(nbrs[j])-particleList.vel(i))/(particleList.rho(nbrs[j])))*kernelViscositySecondDeriv(temp, h); 
	}
		particleList.visc(i) = 1.f*As;
		particleList.vel(i) += particleList.visc(i);
}

//NOT DONE!!!!!!!!!!
void WaterSim::interpolateVorticity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	float As = 0.0f;
	glm::vec3 temp;
	glm::vec3 N;
	glm::vec3 omega;
	glm::vec3 vel;
	glm::vec3 pcm; 
	glm::vec3 eta; 
	glm::vec3 fvort = glm::vec3(0.0f, 0.0f, 0.0f); 
	float epsilon = 0.00001;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		vel = particleList.vel(nbrs[j]) - particleList.vel(i);
		omega = glm::cross(vel,kernelGeneralDeriv(temp, h));
		pcm = (particleList.inv_mass(i)*particleList.predicted_pos(i)
			+particleList.inv_mass(nbrs[j])*particleList.predicted_pos(nbrs[j]))
			/(particleList.inv_mass(i)+particleList.inv_mass(nbrs[j])); 
		eta = pcm-particleList.predicted_pos(i); 
		N = glm::normalize(eta); 
		fvort += (epsilon*glm::cross(N,omega)); 
		//omega += glm::cross(vel,kernelGeneralDeriv(temp, h));
	}
	
	particleList.vort(i) = fvort; 
	particleList.vel(i) += particleList.vort(i);
}

void WaterSim::update(const Scene* const scene, float dt) {
	//For all particles
	for (int i = 0; i < particleList.size(); i++) {
		//Apply forces Vi <- Vi + dt*fext(Xi) 
		//TODO: Are there other external forces??? Does mass come into play???
		particleList.vel(i) = particleList.vel(i) + dt * glm::vec3(0,-9.81,0);

		//Predict position Xi* <- Xi + dt*Vi 
		particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);
	}
	//For all particles
	for (int i = 0; i < particleList.size(); i++) {
		//Find neighboring particles
		particleList.nbr(i) = neighborSearch(i);
	}
	resolve_constraints(scene); 
	for (int i = 0; i < particleList.size(); i++) {
		//Update velocity v = (1/dt) * (Xi* - Xi)
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));

		//Apply vorticity confinement and XSPH viscosity
		interpolateVorticity(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));

		//Update position Xi <- Xi*
		particleList.pos(i) = particleList.predicted_pos(i);
	}
	particleToSphere();
	//updateSpheres();
    m_constraints_ext.clear();
}

void WaterSim::updateSpheres() {
	for (int i = 0; i < dots.size(); i++) {
		dots[i]->m_center = particleList.pos(i);
	}
}

void WaterSim::resolve_constraints(const Scene* const scene) {
	for (unsigned int iter = 0; iter < solverIterations; ++iter) {
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Calculate lambda
			particleList.lambda(i) = calculateLambda(i); 
		}
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Calculate deltaP
			particleList.deltaP(i) = calculateDeltaP(i); 

			//Perform collision detection and response
			collision_detection(scene); 
			resolve_collisions(); 
		}
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Update position Xi* = Xi* + deltaP
			particleList.predicted_pos(i) = particleList.predicted_pos(i) + particleList.deltaP(i);
		}
	}
}

float WaterSim::kernelViscositySecondDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	return (45/(PI*pow(h,6)))*(h-r);
}

glm::vec3 WaterSim::kernelPressureDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	float temp = (45/(PI*pow(h,6)*r)) * pow((h-r),2);
	return -R * temp;
	//-r should be the actual R)
}

glm::vec3 WaterSim::kernelVorticityDeriv(glm::vec3 R, float h) {

}

glm::vec3 WaterSim::kernelGeneralDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	float temp = (945/(32*PI*pow(h,9))) * pow((pow(h,2)-pow(r,2)),2);
	return -R * temp;
}

float WaterSim::kernelGeneralSecondDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	return (945/(8*PI*pow(h,9))) * (pow(h,2) - pow(r,2)) * (pow(r,2) - ((3.f/4.f) * (pow(h,2) - pow(r,2))));
}

void WaterSim::applyXSPH(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	float norm = 0;
	float c = 0.01;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		norm = sqrt(pow(temp.x,2)+pow(temp.y,2)+pow(temp.z,2));
		As += (particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscosity(norm, h); 
	}
		particleList.vel(i) = particleList.vel(i) + c*As;
}