#include "WaterSim.h"
#include <fstream>

#define PI 3.1415926535897

WaterSim::WaterSim() {
	particleNum = 1000;//4096;//8000;//2197;//1728;//
	particleRad = 0.1f; 
	particleList.resize(particleNum);
	m_thick = 0.15f;
	solverIterations = 10;
	xmin = -3.0f;//-5.0f;
	xmax = 3.0f;//5.0f;
	ymin = 0.0f;
	ymax = 6.0f;//10.0f;
	zmin = -3.0f;//-5.0f; 
	zmax = 3.0f;//5.0f; 

	if (particleNum == 8000)
	{
		xmin = -5.0f;
		xmax = 5.0f;
		ymin = 0.0f;
		ymax = 10.0f;
		zmin = -5.0f; 
		zmax = 5.0f; 
	}

	nbrRadius = 0.5f;
	gridXmin = -20.0f;
	gridXmax = 20.0f;
	gridYmin = -20.0f;
	gridYmax = 20.0f;
	gridZmin = -20.0f;
	gridZmax = 20.0f; 
	gridX = int ((gridXmax-gridXmin)/nbrRadius); 
	gridY = int ((gridYmax-gridYmin)/nbrRadius); 
	gridZ = int ((gridZmax-gridZmin)/nbrRadius); 

	// density of water is 1000 kg/m^3 = 1 g/cm^3
	// particle has volume: (4/3)pi*0.1^3 = 0.004189
	// mass = density * volume = 4.189
	rho_0 = 1000; 
}

WaterSim::~WaterSim() {
	m_constraints_ext.clear();
}

void WaterSim::initialize() {
    for(int i = 0; i < particleNum; ++i)
	{
		//Set original velocity of each particle
        particleList.vel(i) = glm::vec3(0.0f);
		//Set original mass of each particle
		particleList.set_mass(i, /*rho_0**/(4.0f/3.0f)*PI*pow(particleRad,3));
    }

	//Set original position of each particle to be a cube of 20x20x20
	int pos = 0;
	float lim = 2; 
	if (particleNum == 8000) 
		lim = 4; 
	for (float i = 0; i < /*4*//*3.2*/lim; i+=(2*particleRad)) 
	{
		for (float j = 0; j < /*4*//*3.2*/lim; j+=(2*particleRad))
		{
			for (float k = 0; k < /*4*//*3.2*/lim; k+=(2*particleRad))
			{
				particleList.pos(pos) = glm::vec3(i+0.5f,j+0.5f,k+0.5f);
				pos++;
			}
		}
	}

	particleToSphere();
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

bool WaterSim::box_intersection(const glm::vec3& p1, const glm::vec3& p2, float threshold, glm::vec3& intersect, glm::vec3& normal) const
{
	float v1, v2;
	// check xmin
	normal = glm::vec3(1.0f, 0.0f, 0.0f); 
    v1 = glm::dot(p1, normal) - xmin; //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = glm::dot(p2, normal) - xmin;
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	// check xmax
	normal = glm::vec3(-1.0f, 0.0f, 0.0f); 
    v1 = xmax - glm::dot(p1, normal); //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = xmax - glm::dot(p2, normal);
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	// check ymin
	normal = glm::vec3(0.0f, 1.0f, 0.0f); 
    v1 = glm::dot(p1, normal) - ymin; //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = glm::dot(p2, normal) - ymin;
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	// check ymax
	normal = glm::vec3(0.0f, -1.0f, 0.0f); 
    v1 = ymax - glm::dot(p1, normal); //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = ymax - glm::dot(p2, normal);
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	// check zmin
	normal = glm::vec3(0.0f, 0.0f, 1.0f); 
    v1 = glm::dot(p1, normal) - zmin; //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = glm::dot(p2, normal) - zmin;
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	// check zmax
	normal = glm::vec3(0.0f, 0.0f, -1.0f); 
    v1 = zmax - glm::dot(p1, normal); //Is basically the plane formula if p1 is on the plane then v1 = 0
    v2 = zmax - glm::dot(p2, normal);
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
            intersect = ((v1 - threshold) * p2 - (v2 - threshold) * p1) / (v1 - v2);
        }
        else
        {// static collision handling.
            intersect = p2 - (v2 - threshold) * normal;
        }
        return true;
    }
	return false; 
}

/*void WaterSim::collision_detection(int id, const Scene* const scene)
{
	// detect collision between water and the scene, and generate external constraints.
	m_constraints_ext.clear();
    glm::vec3 x, p, q, n;
    x = particleList.pos(id);
    p = particleList.predicted_pos(id);
    //if(scene->line_intersection(x, p, m_thick, q, n))
    //{
    //    CollisionConstraint c(&particleList, id, q, n);
    //    m_constraints_ext.push_back(c);
    //}
	if(box_intersection(x, p, m_thick, q, n))
    {
        //CollisionConstraint c(&particleList, id, q, n);
        //c.project_constraint();
		particleList.predicted_pos(id) = q;
    }
}*/

void WaterSim::collision_detection(int id, const Scene* scene)
{
	// detect collision between water and the bounding box, and generate external constraints.
	float thresh = 0.15f; 
	float rest = 0.8f; 
    glm::vec3 x, p, v, temp;
    x = particleList.pos(id);
    p = particleList.predicted_pos(id);
	v = particleList.vel(id); 

    if (p[0]-xmin < thresh)
    {
        particleList.predicted_pos(id) = glm::vec3(xmin+thresh+0.00001f,p[1],p[2]); 
		//particleList.vel(id) = glm::vec3(-rest*v[0],v[1],v[2]);
    }
	if (xmax-p[0] < thresh)
	{
		particleList.predicted_pos(id) = glm::vec3(xmax-thresh-0.00001f,p[1],p[2]); 
		//particleList.vel(id) = glm::vec3(-rest*v[0],v[1],v[2]);
	}
	if (p[1]-ymin < thresh) 
	{ 
		particleList.predicted_pos(id) = glm::vec3(p[0],ymin+thresh+0.00001f,p[2]); 
		//particleList.vel(id) = glm::vec3(v[0],-rest*v[1],v[2]);
	}
	if (ymax-p[1] < thresh) 
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],ymax-thresh-0.00001f,p[2]);
		//particleList.vel(id) = glm::vec3(v[0],-rest*v[1],v[2]);
	}
	if (p[2]-zmin < thresh)
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],p[1],zmin+thresh+0.00001f);
		//particleList.vel(id) = glm::vec3(v[0],v[1],-rest*v[2]);
	}
	if (zmax-p[2] < thresh) 
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],p[1],zmax-thresh-0.00001f);
		//particleList.vel(id) = glm::vec3(v[0],v[1],-rest*v[2]);
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
	for (int i = 0; i < dots.size(); i++) {
		dots[i]->~Sphere();
	}
	dots.clear(); 
	Sphere* s;
	for (int i = 0; i < particleList.size(); i++) {
		s = new Sphere(particleList.pos(i),particleRad);
		dots.push_back(s);
	}
}

void WaterSim::populateGrid() {
	table.clear(); 
	int key, ix, iy, iz;
	glm::vec3 particlePos;
	for (int i = 0; i < particleNum; i++) 
	{
		particlePos = particleList.predicted_pos(i); 
		ix = (int) ceil((particlePos.x-gridXmin)/nbrRadius)-1; 
		iy = (int) ceil((particlePos.y-gridYmin)/nbrRadius)-1; 
		iz = (int) ceil((particlePos.z-gridZmin)/nbrRadius)-1; 

		if (ix < 0) ix = 0; 
		if (iy < 0) iy = 0; 
		if (iz < 0) iz = 0; 

		if (ix > gridX-1) ix = gridX-1; 
		if (iy > gridY-1) iy = gridY-1; 
		if (iz > gridZ-1) iz = gridZ-1;

		if (isValidCell(ix, iy, iz)) 
		{
			key = ix*(gridY*gridZ) + iy*gridZ + iz; 
		//glm::vec3 key = glm::vec3(ix,iy,iz); 
			table[key].push_back(i); 
			
		}
	}
}

std::vector<int> WaterSim::neighborSearch(int id) {
	vector<int> nbrs;
	glm::vec3 particlePos = particleList.predicted_pos(id); 
	int idx = (int) ceil((particlePos.x-gridXmin)/nbrRadius)-1; 
	int idy = (int) ceil((particlePos.y-gridYmin)/nbrRadius)-1; 
	int idz = (int) ceil((particlePos.z-gridZmin)/nbrRadius)-1; 

	if (idx < 0) idx = 0; 
	if (idy < 0) idy = 0; 
	if (idz < 0) idz = 0; 

	if (idx > gridX-1) idx = gridX-1; 
	if (idy > gridY-1) idy = gridY-1; 
	if (idz > gridZ-1) idz = gridZ-1; 

	int indx, indy, indz; 
	float r;
	glm::vec3 nbrPos;
	int nbrid, key;
	std::vector<int> nbr_potential;
	for (int i = -1; i <= 1; i++)
    {
		for(int j = -1; j <= 1; j++)
        {
			for (int k = -1; k <= 1; k++) 
			{
				indx = idx+i; 
				indy = idy+j; 
				indz = idz+k; 

				if (isValidCell(indx,indy,indz)) 
				{
					key = indx*(gridY*gridZ) + indy*gridZ + indz; 
					//glm::vec3 key = glm::vec3(indx,indy,indz);
					nbr_potential = table[key]; 
					for (int a = 0; a < nbr_potential.size(); a++)
					{
						nbrid = nbr_potential[a];
						nbrPos = particleList.predicted_pos(nbrid); 
						r = glm::length(nbrPos - particlePos);  
						if (r < nbrRadius) 
						{
							nbrs.push_back(nbrid);
						}
					}
				}
			}
        }
    }
	if (nbrs.empty()) 
		nbrs.push_back(id); 
    return nbrs;
}

/*std::vector<int> WaterSim::neighborSearch(int id) {
	vector<int> nbrs;
	int nbrid;
	float r;
	glm::vec3 nbrPos;
	glm::vec3 particlePos = particleList.predicted_pos(id); 
	for (int i = 0; i < particleNum; i++) 
	{
		nbrid = i; 
		//if (nbrid != id)
		//{
			nbrPos = particleList.predicted_pos(i); 
			r = glm::length(nbrPos - particlePos);  
			if (r < nbrRadius) 
			{
				nbrs.push_back(nbrid);
			}
		//}
	}  
	return nbrs; 
}*/

bool WaterSim::isValidCell(int i, int j, int k) {
	if (i < 0)
		return false; 
	if (i >= gridX) 
		return false; 
	if (j < 0)
		return false; 
	if (j >= gridY) 
		return false; 
	if (k < 0) 
		return false; 
	if (k >= gridZ) 
		return false; 
	return true; 
}

float WaterSim::calculateLambda(int id) {
	// compute rho_i 
	std::vector<int> Ni = particleList.nbr(id); 
	float rho_i; 
	float Wj; 
	std::vector<float> W; 
	int nbrid;
	float dist2;
	float mj; 
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		dist2 = glm::length(particleList.predicted_pos(id)-particleList.predicted_pos(nbrid)); 
		Wj = kernelGeneral(dist2,nbrRadius); 
		W.push_back(Wj); 
	}
	rho_i = interpolateDensity(particleList.predicted_pos(id), nbrRadius,id,particleList.nbr(id));
	/*rho_i = 0.0f; 
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		mj = particleList.mass(nbrid); 
		rho_i += (mj * W[j]);
	}*/

	// calculate constraint Ci and lambda_i 
	float epsilon = 50.0f; 
	float Ci = rho_i/rho_0 - 1.0f; 
	float lambda_i = 0.0f; 
	std::vector<glm::vec3> delW; 
	glm::vec3 delWj; 
	// spike kernel for gradient calculation 
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		glm::vec3 dist = particleList.predicted_pos(id)-particleList.predicted_pos(nbrid); 
		delWj = kernelPressureDeriv(dist,nbrRadius); 
		delW.push_back(delWj); 
	}
	float sumCiGrad = 0.0f; 
	float nbrterm = 0.0f; 
	glm::vec3 particleterm = glm::vec3(0.0f); 
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		mj = particleList.mass(nbrid); 
		if (nbrid != id) 
		{
		//particleterm = glm::vec3(0.0f);
		//for (int k = 0; k < Ni.size(); k++) 
		//{
			particleterm += mj*(delW[j]); 
			nbrterm += pow(-1.0f*mj/rho_0*glm::length(delW[j]),2); 
		//}
		}
		//sumCiGrad += pow(1.0f/rho_0*glm::length(particleterm),2);
	}
	sumCiGrad = nbrterm + pow(1.0f/rho_0*glm::length(particleterm),2); 
	lambda_i = -Ci/(sumCiGrad+epsilon); 
	return lambda_i; 
}

float WaterSim::calculateDeltaP(int id) {
	//
	std::vector<int> Ni = particleList.nbr(id); 
	float Wj, delWj, Wq, delq; 
	std::vector<float> W; 
	std::vector<float> delW; 
	std::vector<float> pij; 
	std::vector<float> scorr;  
	int nbrid;
	float dist;
	float mj; 
	/*for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		dist = glm::length(particleList.predicted_pos(nbrid)-particleList.predicted_pos(id)); 
		Wj = kernelGeneral(dist,nbrRadius); //kernelPressure(dist,nbrRadius); 
		W.push_back(Wj); 
	} */

	// scorr calculation for tensile stability 
	float k = 0.1f; 
	int n = 4; 
	delq = 0.1f*nbrRadius; 
	Wq = kernelGeneral(delq,nbrRadius);
	//Added by Sam
	glm::vec3 dist2;


	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j];
		//Added by Sam
		dist = glm::length(particleList.predicted_pos(nbrid)-particleList.predicted_pos(id)); 
		Wj = kernelGeneral(dist,nbrRadius);
		/////
		dist/*val*/ = -k*pow(Wj/*W[j]*//Wq,n); 
		scorr.push_back(dist/*val*/);
		dist2 = particleList.predicted_pos(id)-particleList.predicted_pos(nbrid); 
		delWj = glm::length(kernelPressureDeriv(dist2,nbrRadius)); 
		delW.push_back(delWj); 
	}

	/*glm::vec3 dist2;

	// spike kernel for gradient calculation 
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		dist2 = particleList.predicted_pos(id)-particleList.predicted_pos(nbrid); 
		delWj = glm::length(kernelPressureDeriv(dist2,nbrRadius)); 
		delW.push_back(delWj); 
	}*/

	// spike kernel for gradient calculation 
	float mi = particleList.mass(id);
	float deltaP_i = 0.0f;
	for (int j = 0; j < Ni.size(); j++) 
	{
		nbrid = Ni[j]; 
		mj = particleList.mass(nbrid); 
		deltaP_i += ((particleList.lambda(id)+particleList.lambda(nbrid)+scorr[j])*delW[j]); 
	}
	deltaP_i = deltaP_i*mi/rho_0; 
	return deltaP_i; 
}

float WaterSim::kernelGeneral(float r, float h) {
	if (0 <= r && r <= h) {
		return 315.0f/(64.0f*PI*pow(h,9)) * pow((pow(h,2) - pow(r,2)),3);
	} else {
		return 0;
	}
}

float WaterSim::kernelPressure(float r, float h) {
	if (0 <= r && r <= h) {
		return 15.0f/(PI*pow(h,6)) * pow((h-r),3);
	} else {
		return 0;
	}
}

float WaterSim::kernelViscosity(float r, float h) {
	if (0 <= r && r <= h) {
		return 15.0f/(2.0f*PI*pow(h,3)+0.00001) * (-(pow(r,3)/(2.0f*pow(h,3)+0.00001)) + (pow(r,2)/(pow(h,2)+0.00001)) + (h/(2*r+0.00001)) - 1.0f);
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
		As += particleList.mass(nbrs[j])*kernelGeneral(norm, h);
	}
	return As; 
}

void WaterSim::interpolatePressure(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		As += particleList.mass(nbrs[j])*((particleList.press(i)+particleList.press(nbrs[j]))/(2.0f*particleList.rho(nbrs[j]))) * kernelPressureDeriv(temp,h);
	}
		particleList.press(i) = -As;
}

void WaterSim::interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	float mu = 1.0f; 
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		As += particleList.mass(nbrs[j])*((particleList.vel(nbrs[j])-particleList.vel(i))/(particleList.rho(nbrs[j])))*kernelViscositySecondDeriv(temp, h); 
	}
		particleList.visc(i) = mu*As;
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
	float epsilon = 0.000001;
	float temp2;
	for(int j = 0; j < nbrs.size(); j++) {
		temp = R-particleList.predicted_pos(nbrs[j]);
		vel = particleList.vel(nbrs[j]) - particleList.vel(i);
		omega = glm::cross(vel,kernelGeneralDeriv(temp, h));
		if (particleList.mass(i)+particleList.mass(nbrs[j]) != 0) {
			pcm = (particleList.mass(i)*particleList.predicted_pos(i)
				+particleList.mass(nbrs[j])*particleList.predicted_pos(nbrs[j]))
				/(particleList.mass(i)+particleList.mass(nbrs[j])+ glm::vec3(0.000000001,0.000000001,0.000000001));
		} else {
			pcm = glm::vec3(0,0,0);
		}
		eta = pcm-particleList.predicted_pos(i);
		temp2 = (glm::length(eta)+0.000001);
		N = eta/(temp2+0.000001f);//glm::normalize(eta); 
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
	populateGrid(); 
	for (int i = 0; i < particleList.size(); i++) {
		//Find neighboring particles
		particleList.nbr(i) = neighborSearch(i);
	}
	resolve_constraints(scene); 
	for (int i = 0; i < particleList.size(); i++) {
		//Update velocity v = (1/dt) * (Xi* - Xi)
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));

		//Apply vorticity confinement and XSPH viscosity
		//With just this nothing happens due to vorticity, particles hit floor and stay there
		//interpolateVorticity(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		//With just this column of bubbles happens
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));

		//collision_detection(i, scene);

		//Update position Xi <- Xi*
		particleList.pos(i) = particleList.predicted_pos(i);
	}
	particleToSphere();
	//updateSpheres();
    //m_constraints_ext.clear();
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
			//With this, the particles explode
			particleList.deltaP(i) = calculateDeltaP(i); 

			//Perform collision detection and response
			//collision_detection(scene); 
			//resolve_collisions(); 
			/*std::ostringstream stream;	
			stream << i << std::endl;	
			std::cout << stream.str();	
			fflush(stdout); */
			//collision_detection(i, scene);
			//resolve_collisions(); 
		}
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Update position Xi* = Xi* + deltaP
			particleList.predicted_pos(i) = particleList.predicted_pos(i) + particleList.deltaP(i);
		}
		for (int i = 0; i < particleList.size(); i++) {
			collision_detection(i, scene);
		}
	}
}

float WaterSim::kernelViscositySecondDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	return (45/(PI*pow(h,6)))*(h-r);
}

glm::vec3 WaterSim::kernelPressureDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	float temp = (45.0f/(PI*pow(h,6)*(r+0.0000001))) * pow((h-r),2);
	return -R * temp;
	//-r should be the actual R)
}

glm::vec3 WaterSim::kernelVorticityDeriv(glm::vec3 R, float h) {
	return glm::vec3(0.0f);
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
		norm = glm::length(temp);//sqrt(pow(temp.x,2)+pow(temp.y,2)+pow(temp.z,2));
		//As += (particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscosity(norm, h); 
		As += /*(particleList.mass(j)/(0.00001f + particleList.rho(j)))**/(particleList.vel(nbrs[j])-particleList.vel(i))*kernelGeneral(norm,h);//SecondDeriv(temp, h); 
	}
		particleList.vel(i) = particleList.vel(i) + c*As;
}

vector<int> WaterSim::surfaceSampling(vector<int> levelSet, float r, float t, float e) {
//for all grid cells C where  changes sign do
	//for (int i = 0; i < ; i++) {
//	for t attempts do
//		Generate random point p in C
//		Project p to surface of phi
//		if p meets the Poisson Disk criterion in S then
//			S <- S U {p}
//			Break
//	if no point was found in C then
//		Continue
//	while new samples are found do
//		Generate random tangential direction d to surface at p
//		q <- p + d  e  r
//		Project q to surface of 
//		if q meets the Poisson Disk criterion in S then
//			S <- S U {q}
//			p <- q
}