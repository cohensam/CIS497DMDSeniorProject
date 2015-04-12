/*TO DO:
1) Initialize all particleList member functions you just added - 1/1
2) Check pressure/spiky vs. general kernels for all equations - ???
3) Debug and check parameters! - ???
4) Once working, start incorporating high and low res particle paper - in school
5) Once working, do water and sand interaction paper - in school
6) Download LaTex
7) Type all equations used with LaTex
*/


#include "WaterSim.h"
#include <fstream>


#define PI 3.1415926535897

WaterSim::WaterSim() {

	Wthreshold = 0.05f;//0.5f;
	Wmax = 0.1f;//1.0f;
	rhoMax = 2000.0f;
	//Degrees, for dry loose beach sand
	angleRepose = 38.0f;

	particleNum = 1000;//6000;//2000+7000;//2000;//1000;//2197;//1728;//1000;//4096;//
	granularParticleNum = 1000;
	particleRad = 0.1f; //0.05f; 
	m_thick = 0.05f+particleRad;
	solverIterations = 10;
	xmin = -3.0f;//-5.0f;
	xmax = 3.0f;//5.0f;
	ymin = 0.0f;
	ymax = 7.0f;//10.0f;
	zmin = -3.0f;//-5.0f; 
	zmax = 3.0f;//5.0f; 

	/*ni = 60;
	nj = 60;
	nk = 60;
	phi_dx = 1 / (float)ni;   
	nodal_solid_phi.resize(ni+1,nj+1,nk+1); 
	liquid_phi.resize(ni,nj,nk);*/

	/*if ((particleNum == 8000 && particleRad == 0.1f) || particleNum == 8000+3375 || particleNum == 16000)
	{
		xmin = -3.0f;
		xmax = 4.0f;
		ymin = 0.0f;
		ymax = 8.0f;
		zmin = -3.0f; 
		zmax = 7.0f; 
	} else if (particleNum == 10648) {
		xmin =  -6.0f;//-10.0f;
		xmax =  6.0f;//10.0f;
		ymin =  0.0f;//0.0f;
		ymax =  10.0f;//20.0f;
		zmin =  -6.0f;//-10.0f; 
		zmax =  6.0f; //10.0f; 
	}*/

	nbrRadius = 3*particleRad;
	gridXmin = -20.0f;
	gridXmax = 20.0f;
	gridYmin = -20.0f;
	gridYmax = 20.0f;
	gridZmin = -20.0f;
	gridZmax = 20.0f; 
	gridX = int ((gridXmax-gridXmin)/nbrRadius); 
	gridY = int ((gridYmax-gridYmin)/nbrRadius); 
	gridZ = int ((gridZmax-gridZmin)/nbrRadius); 
	
	//1000 for water
	//2000 for dry sand
	rho_0 = 2000; 

	particleList.resize(particleNum+7*particleNum+1000);

	 //density of water is 1000 kg/m^3 = 1 g/cm^3
	 //particle has volume: (4/3)pi*0.1^3 = 0.004189
	 //mass = density * volume = 4.189
}

WaterSim::~WaterSim() {
	m_constraints_ext.clear();
}

void WaterSim::initialize() {
    for(int i = 0; i < /*particleNum*/1000; ++i)
	{
		//Set original velocity of each particle
        particleList.vel(i) = glm::vec3(0.0f);
		//Set original mass of each particle
		// density of water is 1000 kg/m^3 = 1 g/cm^3
		// particle has volume: (4/3)pi*0.1^3 = 0.004189
		// mass = density * volume = 4.189
		
		//particleList.set_mass(i, (4.0f/3.0f)*PI*pow(particleRad,3));
		particleList.cluster_num(i) = -1;
		particleList.type(i) = 1;
		rho_0 = 2000;
		particleList.set_mass(i, rho_0*(4.0f/3.0f)*PI*pow(particleRad,3));
    }

	//Set original position of each particle to be a cube of 20x20x20
	int pos = 0;
	float lim = 2*10*particleRad; 
	if (particleNum == 8000 || particleNum == 8000+3375 || particleNum == 16000) {
		lim = 4*10*particleRad; 
	}
	for (float i = 0; i < /*4*//*3.2*/lim; i+=(2*particleRad)) 
	{
		for (float j = 0; j < /*4*//*3.2*/lim; j+=(2*particleRad))
		{
			for (float k = 0; k < /*4*//*3.2*/lim; k+=(2*particleRad))
			{
				particleList.pos(pos) = glm::vec3(i-0.5f,j+0.2f,k+1.0f);//glm::vec3(i-0.5f,j+4.15f,k-0.5f);0.2
				pos++;
			}
		}
	}

	//Set up high res particles
	int id = granularParticleNum;
	for (int i = 0; i < 1000; i++) { 
		for (int j = 0; j < 7; j++) {
			particleList.vel(id+j) = glm::vec3(0.0f);
			particleList.set_mass(id+j, rho_0*(4.0f/3.0f)*PI*pow(particleRad,3));
			particleList.type(id+j) = 5;
			particleList.cluster_num(id+j) = -1;
		}
		particleList.pos(id) = glm::vec3(particleList.pos(i).x, particleList.pos(i).y, particleList.pos(i).z-particleRad);
		particleList.pos(id+1) = glm::vec3(particleList.pos(i).x+particleRad, particleList.pos(i).y, particleList.pos(i).z);
		particleList.pos(id+2) = glm::vec3(particleList.pos(i).x-particleRad, particleList.pos(i).y, particleList.pos(i).z);
		particleList.pos(id+3) = glm::vec3(particleList.pos(i).x, particleList.pos(i).y+particleRad, particleList.pos(i).z);
		particleList.pos(id+4) = glm::vec3(particleList.pos(i).x, particleList.pos(i).y-particleRad, particleList.pos(i).z);
		particleList.pos(id+5) = glm::vec3(particleList.pos(i).x, particleList.pos(i).y, particleList.pos(i).z+particleRad);
		particleList.pos(id+6) = particleList.pos(i);

		id+=7;
	}

	//SET UP WATER PARTICLES
	for (int i = id; i < particleList.size(); i++) {
		particleList.type(i) = 4;
		rho_0 = 1000;
		particleList.set_mass(i, rho_0*(4.0f/3.0f)*PI*pow(particleRad,3));
		particleList.wet(i) = 1.0f;
	}

	lim = 2*10*particleRad; 
	for (float i = 0; i < /*4*//*3.2*/lim; i+=(2*particleRad)) 
	{
		for (float j = 0; j < /*4*//*3.2*/lim; j+=(2*particleRad))
		{
			for (float k = 0; k < /*4*//*3.2*/lim; k+=(2*particleRad))
			{
				particleList.pos(id) = glm::vec3(i-0.5f,j+4.15f,k-2.0f);//glm::vec3(i-0.5f,j+4.15f,k-0.5f);0.2
				id++;
			}
		}
	}

	particleToSphere();
}

void WaterSim::compute_mesh_bounds()
{
	tri_xmax = triangles[0]->p0.x; 
	tri_xmin = triangles[0]->p0.x; 
	tri_ymax = triangles[0]->p0.y; 
	tri_ymin = triangles[0]->p0.y; 
	tri_zmax = triangles[0]->p0.z; 
	tri_zmin = triangles[0]->p0.z; 

	glm::vec3 v0, v1, v2; 
	float temp; 
	// calculate axis-aligned mesh limits 
	for (int i = 0; i < triangles.size(); i++) 
	{
		v0 = triangles[i]->p0; 
		v1 = triangles[i]->p1; 
		v2 = triangles[i]->p2; 
		temp = max(v0.x, max(v1.x, v2.x)); 
		if (temp > tri_xmax) 
			tri_xmax = temp; 
		temp = min(v0.x, min(v1.x, v2.x)); 
		if (temp < tri_xmin) 
			tri_xmin = temp; 

		temp = max(v0.y, max(v1.y, v2.y)); 
		if (temp > tri_ymax) 
			tri_ymax = temp; 
		temp = min(v0.y, min(v1.y, v2.y)); 
		if (temp < tri_ymin) 
			tri_ymin = temp; 

		temp = max(v0.z, max(v1.z, v2.z)); 
		if (temp > tri_zmax) 
			tri_zmax = temp; 
		temp = min(v0.z, min(v1.z, v2.z)); 
		if (temp < tri_zmin) 
			tri_zmin = temp; 
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

void WaterSim::collision_detection(int id, const Scene* scene)
{
	// object collision 
	//obj_collision(id); 

	// detect collision between water and the bounding box
	float thresh = m_thick; 
	float rest = 0.5f; 
    glm::vec3 x, p, v, temp;
    x = particleList.pos(id);
    p = particleList.predicted_pos(id);
	v = particleList.vel(id); 

    if (p[0]-xmin < thresh)
    {
        particleList.predicted_pos(id) = glm::vec3(xmin+thresh+0.00001f,p[1],p[2]); 
		particleList.vel(id) = glm::vec3(-rest*v[0],v[1],v[2]);
    }
	v = particleList.vel(id); 
	p = particleList.predicted_pos(id);
	if (xmax-p[0] < thresh)
	{
		particleList.predicted_pos(id) = glm::vec3(xmax-thresh-0.00001f,p[1],p[2]); 
		particleList.vel(id) = glm::vec3(-rest*v[0],v[1],v[2]);
	}
	v = particleList.vel(id); 
	p = particleList.predicted_pos(id);
	if (p[1]-ymin < thresh) 
	{ 
		particleList.predicted_pos(id) = glm::vec3(p[0],ymin+thresh+0.00001f,p[2]); 
		particleList.vel(id) = glm::vec3(v[0],-rest*v[1],v[2]);
	}
	v = particleList.vel(id); 
	p = particleList.predicted_pos(id);
	if (ymax-p[1] < thresh) 
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],ymax-thresh-0.00001f,p[2]);
		particleList.vel(id) = glm::vec3(v[0],-rest*v[1],v[2]);
	}
	v = particleList.vel(id); 
	p = particleList.predicted_pos(id);
	if (p[2]-zmin < thresh)
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],p[1],zmin+thresh+0.00001f);
		particleList.vel(id) = glm::vec3(v[0],v[1],-rest*v[2]);
	}
	v = particleList.vel(id); 
	p = particleList.predicted_pos(id);
	if (zmax-p[2] < thresh) 
	{
		particleList.predicted_pos(id) = glm::vec3(p[0],p[1],zmax-thresh-0.00001f);
		particleList.vel(id) = glm::vec3(v[0],v[1],-rest*v[2]);
	}
	//obj_collision(id); 
}

void WaterSim::obj_collision(int id) {
	glm::vec3 x, p, v;
	x = particleList.pos(id);
	p = particleList.predicted_pos(id);
	v = particleList.vel(id); 
	if (check_obj_bound(x) || check_obj_bound(p))
	{
		for (int i = 0; i < triangles.size(); i++) 
		{
			//triangle_collision(i, id); 
			// assume only one collision with obj mesh
			if (triangle_collision(i, id))
				break; 
		}
	}
}

bool WaterSim::check_obj_bound(const glm::vec3 p)
{	
	float thresh = m_thick; 
	if (tri_xmin-p.x < thresh && p.x-tri_xmax < thresh)
	{
		if (tri_ymin-p.y < thresh && p.y-tri_ymax < thresh)
		{
			if (tri_zmin-p.z < thresh && p.z-tri_zmax < thresh)
				return true; 
		}
	}
	return false; 
}

bool WaterSim::triangle_collision(int tri, int id) {
	float thresh = m_thick; 
	float rest = 0.2f; 
	float fric = 0.3f; 
	glm::vec3 p1 = particleList.pos(id);
	glm::vec3 p2 = particleList.predicted_pos(id);
	glm::vec3 v = particleList.vel(id); 
	glm::vec3 p12 = glm::normalize(p2-p1); 
	float t2 = glm::length(p2-p1);  
	glm::vec3 v0 = triangles[tri]->p0; 
	glm::vec3 v1 = triangles[tri]->p1; 
	glm::vec3 v2 = triangles[tri]->p2; 
	glm::vec3 n = triangles[tri]->norm; 
	//std::ostringstream stream;	
	//stream << v0.x << std::endl;	
	//std::cout << stream.str();	
	//fflush(stdout); 

	// determine intersection point of plane
	float d = glm::dot(n, v0); 
	float t; 
	/*if (glm::dot(n, p12) == 0)
		t = -1.0f; 
	else */
	t = (d - glm::dot(n, p1)) / (glm::dot(n, p12)+0.000001f); 

	float a1, a2;
	glm::vec3 q; 
	bool intersect = false; 
    a1 = glm::dot(p1, n) - d; //Is basically the plane formula if p1 is on the plane then v1 = 0
    a2 = glm::dot(p2, n) - d;
    if(a2 < thresh)
    {
        if(a1 >= thresh)
        {// continuous collision handling.
            q = ((a1 - thresh) * p2 - (a2 - thresh) * p1) / (a1 - a2);
			intersect = true; 
        }
        else
        {// static collision handling.
            q = p2 - (a2 - thresh) * n;
			intersect = true; 
        }
        //return true;
    }
    //else
    //    return false;

	// point does not intesect plane of triangle
	/*if (t < 0.0f || t > t2+thresh) 
		return false; */

	// check barycentric coordinates to see if point falls inside triangle
	//glm::vec3 q = p1 + t*p12; 
	glm::vec3 temp; 

	temp = glm::cross(v1-v0, q-v0); 
	if (glm::dot(n, temp) < 0) 
		return false; 

	temp = glm::cross(v2-v1, q-v1); 
	if (glm::dot(n, temp) < 0) 
		return false; 

	temp = glm::cross(v0-v2, q-v2); 
	glm::vec3 vt, vn; 
	if (glm::dot(n, temp) < 0) 
		return false; 
	else 
	{
		if (intersect && (t > 0.0f && t <= t2+thresh))
		{
			particleList.predicted_pos(id) = p1 + (t-thresh)*p12;
			vn = glm::dot(v, n)*n; 
			vt = v-vn; 
			particleList.vel(id) = -rest*vn + (1.0f-fric)*vt; 
			//particleList.vel(id) = rest*(-2.0f*glm::dot(v, n)*n+v);
			return true; 
		}
		else 
			return false; 
	}
}

void WaterSim::objLoader(string filename) {	
	triangles.clear();
	ifstream inFile(filename, ifstream::in);
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
	glm::vec3 tempy = glm::vec3(0.0f);
	dots.resize(particleList.size());//particleNum+7*particleNum+1000);
	for (int i = 0; i < dots.size(); i++) {
		
		tempy = glm::vec3(abs(particleList.vel(i).x/glm::length(particleList.vel(i))), abs(particleList.vel(i).y/glm::length(particleList.vel(i))), abs(particleList.vel(i).z/glm::length(particleList.vel(i))));
		//glm::vec3 tempy = glm::vec3(0.0f, 0.0f, 1.0f); 
		if (particleList.type(i)==5) { //if type == 5
			s = new Sphere(particleList.pos(i),particleRad/2.0f, tempy);
		} else {
			s = new Sphere(particleList.pos(i),particleRad, tempy);
		}
		dots[i] = s;//.push_back(s);
	}
	/*Sphere* s2 = new Sphere(glm::vec3(0,1,0),1.75f,glm::vec3(0,0,1));
	dots.push_back(s2);
	Sphere* s3 = new Sphere(glm::vec3(0,3.75f,1),0.5f,glm::vec3(0,1,0));
	dots.push_back(s3);
	Sphere* s4 = new Sphere(glm::vec3(0,3.75f,1.5f),0.5f,glm::vec3(0,1,0));
	dots.push_back(s4);*/
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
			//grid will contain all keys for different grid points, look-up-able in table
			//grid[ix][iy][iz] = key;
			
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
	std::vector<int> tempv;
	tempv.resize(particleList.size());
	int counter = 0;
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
							tempv[counter] = nbrid;
							counter++;
							//nbrs.push_back(nbrid);
						}
					}
				}
			}
        }
    }
	nbrs.resize(counter);
	for (int i = 0; i < counter; i++) {
		nbrs[i] = tempv[i];
	}
	if (nbrs.empty()) {
		nbrs.resize(1);
		nbrs[0] = id;
		//nbrs.push_back(id);
	}
    return nbrs;
}

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
	if (particleList.type(id) == 4) {
		rho_0 = 1000;
	} else {
		rho_0 = 2000;
	}

	// compute rho_i 
	std::vector<int> Ni = particleList.nbr(id); 
	float rho_i = 0.0; 
	float Wj = 0.0; 
	std::vector<float> W; 
	int nbrid = 0;
	float dist2 = 0.0;
	float mj = 0.0; 
	rho_i = particleList.rho(id);//interpolateDensity(particleList.predicted_pos(id), nbrRadius,id,particleList.nbr(id));

	// calculate constraint Ci and lambda_i 
	float epsilon = 40000.0f; 
	float Ci = rho_i/rho_0 - 1.0f; 
	float lambda_i = 0.0f; 
	std::vector<glm::vec3> delW; 
	glm::vec3 delWj; 
	// spike kernel for gradient calculation 
	float sumCiGrad = 0.0f; 
	float nbrterm = 0.0f; 
	glm::vec3 particleterm = glm::vec3(0.0f); 
	glm::vec3 dist = glm::vec3(0.0f);
	glm::vec3 tempy = glm::vec3(0.0f);
	for (int j = 0; j < Ni.size(); j++) 
	{
		dist = particleList.predicted_pos(id)-particleList.predicted_pos(nbrid);
		nbrid = Ni[j]; 
		mj = particleList.mass(nbrid); 
		if (nbrid != id) 
		{ 
			tempy = kernelPressureDeriv(dist,nbrRadius);
			particleterm += mj*(tempy); 
			nbrterm += pow(-1.0f*mj/rho_0*glm::length(tempy),2); 
		}
	}
	sumCiGrad = nbrterm + pow(1.0f/rho_0*glm::length(particleterm),2); 
	lambda_i = -Ci/(sumCiGrad+epsilon); 
	return lambda_i; 
}

glm::vec3 WaterSim::calculateDeltaP(int id) {
	if (particleList.type(id) == 4) {
		rho_0 = 1000;
	} else {
		rho_0 = 2000;
	}
	//
	std::vector<int> Ni = particleList.nbr(id); 
	float Wj, Wq, delq; 
	glm::vec3 delWj;
	std::vector<float> W; 
	std::vector<glm::vec3> delW; 
	delW.resize(Ni.size());
	std::vector<float> scorr;  
	scorr.resize(Ni.size());
	int nbrid;
	float dist;
	float mj; 

	// scorr calculation for tensile stability 
	float k = 0.01f; 
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
		dist = -k*pow(Wj/(Wq+0.00001),n); 
		scorr[j] = dist;
		dist2 = particleList.predicted_pos(id)-particleList.predicted_pos(nbrid); 
		delWj = kernelPressureDeriv(dist2,nbrRadius);
		delW[j] = delWj;
	}

	// spike kernel for gradient calculation 
	float mi = particleList.mass(id);
	glm::vec3 deltaP_i = glm::vec3(0,0,0);//0.0f;
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
		return 315.0f/((64.0f*PI*pow(h,9.0f))+0.000001f) * pow((pow(h,2.0f) - pow(r,2.0f)),3.0f);
	} else {
		return 0;
	}
}

float WaterSim::kernelPressure(float r, float h) {
	if (0 <= r && r <= h) {
		return 15.0f/((PI*pow(h,6))+0.000001f) * pow((h-r),3.0f);
	} else {
		return 0;
	}
}

float WaterSim::kernelViscosity(float r, float h) {
	if (0 <= r && r <= h) {
		return 15.0f/((2.0f*PI*pow(h,3))+0.000001f) * (-(pow(r,3)/((2.0f*pow(h,3))+0.000001f)) + (pow(r,2)/(pow(h,2)+0.000001f)) + (h/((2*r)+0.000001f)) - 1.0f);
	} else {
		return 0;
	}
}

float WaterSim::interpolateDensity(glm::vec3 R, float h, int i, vector<int> nbrs) {	
	float As = 0.0f;
	glm::vec3 temp;
	float norm = 0;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			norm = glm::length(temp);
			As += particleList.mass(nbrs[j])*kernelGeneral(norm, h);
		}
	}
	particleList.rho(i) = As; 
	return As; 
}

void WaterSim::interpolatePressure(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As = glm::vec3(0,0,0);
	glm::vec3 temp;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			//As += particleList.mass(nbrs[j])*((particleList.press(i)+particleList.press(nbrs[j]))/(2.0f*particleList.rho(nbrs[j])+0.000001f)) * kernelPressureDeriv(temp,h);
			if ((particleList.type(nbrs[j]) != 4 && particleList.type(i) != 4) || (particleList.type(nbrs[j]) == 4 && particleList.type(i) == 4)) {
				//IF THE PARTICLES ARE BOTH WATER OR BOTH GRANULAR
				As += particleList.mass(nbrs[j])*((particleList.press(i)/(pow(particleList.rho(i),2.0f)+0.000001f)) + 
					(particleList.press(nbrs[j])/(pow(particleList.rho(nbrs[j]),2.0f)+0.000001f)))* kernelPressureDeriv(temp,h);
			} else { //WATER AND GRANULAR PARTICLES INTERACTING
				if (particleList.type(i) == 4) {
					As += ((particleList.mass(i)/particleList.rho(i))*(particleList.mass(nbrs[j])/particleList.rho(nbrs[j]))/(4.0f))*(particleList.press(i)+particleList.press(nbrs[j]))*kernelPressureDeriv(temp,h);
				} else {
					As += -1.0f * ((particleList.mass(i)/particleList.rho(i))*(particleList.mass(nbrs[j])/particleList.rho(nbrs[j]))/(4.0f))*(particleList.press(i)+particleList.press(nbrs[j]))*kernelPressureDeriv(temp,h);
				}
			}
		}	
	}
	
	particleList.press_force(i) = -particleList.mass(i)*As;
	//for water only we did:
	//particleList.press(i) = -As;
}

void WaterSim::interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As = glm::vec3(0,0,0);
	glm::vec3 temp;
	//viscosity = shear stress/shear rate
	//shear stress = I HAVE THIS
	//shear rate = velocity/distance

	//mu between 0 and 1
	//mu = viscosity of water
	//0.894 = viscosity of water, 
	float mu = 0.894f;//0.894;//6f;//1.0f; 
	if (particleList.type(i) != 4) {
		mu = 0.6f;//particleList.stress(i);
	}
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			if ((particleList.type(nbrs[j]) != 4 && particleList.type(i) != 4) || (particleList.type(nbrs[j]) == 4 && particleList.type(i) == 4)) {
				As += particleList.mass(nbrs[j])*((particleList.vel(nbrs[j])-particleList.vel(i))/(particleList.rho(nbrs[j])+0.000001f))*kernelViscositySecondDeriv(temp, h); 
			} else {
				if (particleList.type(i) == 4) {
					As += (particleList.mass(i)/(particleList.rho(i)+0.000001f))*(particleList.mass(nbrs[j])/(particleList.rho(nbrs[j])+0.000001f))/(2.0f)*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscositySecondDeriv(temp, h);
				} else {
					As += -1.0f * (particleList.mass(i)/(particleList.rho(i)+0.000001f))*(particleList.mass(nbrs[j])/(particleList.rho(nbrs[j])+0.000001f))/(2.0f)*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscositySecondDeriv(temp, h);
				}
			}
		}
	}
	//viscous force
	particleList.visc(i) = mu*As;
}

void WaterSim::interpolateFriCoh(glm::vec3 R, float h, int i, vector<int> nbrs) {
	//CAUSES ERRORS
	glm::vec3 As = glm::vec3(0,0,0);
	glm::vec3 temp = glm::vec3(0,0,0);
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j]) * (particleList.stress(i)/(pow(particleList.rho(i),2.0f)+0.000001f) + particleList.stress(nbrs[j])/(pow(particleList.rho(nbrs[j]),2.0f)+0.000001f)) * kernelGeneralDeriv(temp,h);
		}
	}
	particleList.fricoh(i) = particleList.mass(i)*As;
}

void WaterSim::calculateCorrectivePressure(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	if (particleList.type(i) == 4) {
		rho_0 = 1000;
	} else {
		rho_0 = 2000;
	}
	//This is what is making it crazy when we just have gravity, visc, and pressure
	//You need to check to see if we are using the general or the pressure kernel
	glm::vec3 As1;
	float As2 = 0;
	glm::vec3 temp;
	glm::vec3 kernelCalc;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			kernelCalc = kernelPressureDeriv(temp,h);
			As1 += kernelCalc;
			As2 += glm::dot(kernelCalc,kernelCalc);
		}
	}
	particleList.corrective_press(i) = (particleList.rho(i)-rhoMax)*(pow(rho_0,2.0f)/((2.0f*pow(particleList.mass(i),2.0f)+0.000001f)*pow(dt,2.0f)*(glm::dot(As1,As1)+As2)));//-As1*As1-As2);
}

void WaterSim::calculateDiscreteParticleForces(int i, vector<int> nbrs) {
	glm::vec3 Vn = glm::vec3(0,0,0);
	glm::vec3 Fs = glm::vec3(0,0,0);
	glm::vec3 Fd = glm::vec3(0,0,0);
	float Kt = 0.6f;
	float Ks = 0.000001f;

	float R1 = particleRad;
	float R2 = particleRad;
	
	float d = R1+R2;

	glm::vec3 X1 = particleList.predicted_pos(i);
	glm::vec3 X2 = glm::vec3(0,0,0);
	float zeta = 0.0f;
	glm::vec3 N = glm::vec3(0,0,0);
	glm::vec3 V = glm::vec3(0,0,0);
	glm::vec3 V1 = particleList.predicted_vel(i);
	glm::vec3 V2 = glm::vec3(0,0,0);
	glm::vec3 zetaDot = glm::vec3(0,0,0);
	glm::vec3 Vt = glm::vec3(0,0,0);
	glm::vec3 fn = glm::vec3(0,0,0);
	float Meff = 0;
	//Coefficient of restitution - avg. ish and for 38 degree ish 0.6
	float e = 0.6f;
	//Contact duration-this should be small but I am not sure quite how small so I guessed
	float tc = 10000.0f;//100;//0.01f;
	float Kd = 0;
	float Kr = 0;
	//Coefficient of friction for sand, number for general use 0.6
	float mu = 0.6f;
	glm::vec3 Fn = glm::vec3(0,0,0);
	glm::vec3 Ft = glm::vec3(0,0,0);
	glm::vec3 F = glm::vec3(0,0,0);
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4) {

			//Meff = particleList.mass(i)*particleList.mass(nbrs[j])/(particleList.mass(i)+particleList.mass(nbrs[j])+0.000001f);
			Kd = 1.0f;//2.0f*Meff*(-log(e)/tc);

			X2 = particleList.predicted_pos(nbrs[j]);
			V2 = particleList.predicted_vel(nbrs[j]);
			N = X1-X2;
			N /= (glm::length(N)+0.000001f);
			V = V1-V2;
			Vn = glm::dot(V,N)*N;
			Vt = V-Vn;
			Fs = Ks*(d-glm::length(X1-X2))*(Vn/(glm::length(Vn)+0.000001f));
			Fd = Kd*(1.0f+particleList.wet(i)+particleList.wet(nbrs[j]))*Vn;
			Ft = Kt*(1.0f+particleList.wet(i)+particleList.wet(nbrs[j]))*(Vt/(glm::length(Vt)+0.000001f));
			Fn = Fs+Fd;
			F += Fn+Ft;
			//Calulate Fn
			/*X2 = particleList.predicted_pos(nbrs[j]);
			V2 = particleList.predicted_vel(nbrs[j]);
			zeta = max(0.0f,glm::length(R1+R2-(X1-X2)));
			N = X1-X2;
			N /= (glm::length(N)+0.000001f);
			V = V1-V2;
			zetaDot = V * N;
			Vt = V-zetaDot*N;
			Meff = particleList.mass(i)*particleList.mass(nbrs[j])/(particleList.mass(i)+particleList.mass(nbrs[j])+0.000001f);
			Kd = 2.0f*Meff*(-log(e)/tc);
			Kr = (Meff/(pow(tc,2.0f)+0.000001f))*(pow(log(e),2.0f)+pow(3.14159f,2.0f));
			fn = -Kd*pow(zeta,0.5f)*zetaDot-Kr*pow(zeta,3.0f/2.0f);
			Fn = fn*(1.0f+particleList.wet(i)+particleList.wet(j))*N;

			//Calculate Ft
			Vt /= (glm::length(Vt)+0.000001f);
			Ft = -mu*fn*(1.0f+particleList.wet(i)+particleList.wet(j))*Vt;

			//Calculate F
			F += Fn + Ft;*/
		}
	}
	//Fn and Ft are undefined
	particleList.disc(i) = F;//glm::vec3(1,1,1);//F;

}

void WaterSim::calculateStrainRate(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	//Check general vs. pressure kernel
	glm::vec3 deltaVi;
	for(int j = 0; j < nbrs.size(); j++) {
		deltaVi += particleList.mass(nbrs[j])/(particleList.rho(nbrs[j]))*kernelGeneralDeriv(R-particleList.predicted_pos(nbrs[j]),h)*particleList.predicted_vel(nbrs[j]);
	}

	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j])/(particleList.rho(nbrs[j])+0.000001f) * kernelGeneralDeriv(temp,h) * particleList.predicted_vel(nbrs[j]) - deltaVi;
		}
	}
	//CHECK THIS MATH WITH ADDING THE BOUNDARY U
	particleList.strain(i) = -As;//-interpolateU(R,h,i,nbrs);//-(As - (interpolateU(R,h,i,nbrs) + interpolateBoundaryU(R,h,i,nbrs)));
}

void WaterSim::calculateCorrectiveStress(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	//You need to check the whole transpose thing-is that just negative in these cases or can I ignore it?
	//Should the kernel calc be squared or not???
	//ISSUE SOMEWHERE HERE
	glm::vec3 D;
	glm::vec3 temp;
	glm::vec3 kernelCalc;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			kernelCalc = kernelGeneralDeriv(temp,h);
			D += 1.0f/(particleList.rho(nbrs[j])+0.000001f)*kernelCalc*kernelCalc;
		}
	}
	D = 2.0f*pow(particleList.mass(i),2.0f)*dt/(pow(particleList.rho(i),2.0f)+0.000001f) * D;
	particleList.corrective_stress(i) = glm::dot(-D,particleList.strain(i));
	//cout << "Kernel " << kernelCalc.x << kernelCalc.y << kernelCalc.z << endl;
	//cout << "D " << D.x << D.y << D.z << endl;
	//cout << "particleList.strain(i) " << particleList.strain(i).x << particleList.strain(i).y << particleList.strain(i).z << endl;
	//cout << "particleList.corrective_stress(i) " << particleList.corrective_stress(i).x << particleList.corrective_stress(i).y << particleList.corrective_stress(i).z << endl;
}

void WaterSim::testYieldandCohesion(int i, float dt) {
	if (particleList.type(i) != 4 && particleList.type(i) != 0) {
		/*I do not understand how any of this actually is meant to be used
		to change the stress value*/
		float alpha = 0.2f;//0.5f;//sqrt(2.0f/3.0f)*sin(angleRepose);
		/*if (glm::length(particleList.stress(i)) <= alpha*glm::length(particleList.press(i))) {
			//particleList.stress(i) = glm::vec3(0,0,0);//alpha*particleList.press(i);
			particleList.rigid(i) = 1;
		}*/
		//PRESSURE SEEMS TO CONSISTENTLY BE 0
		if (particleList.stress(i) <= alpha*particleList.press(i) && particleList.type(i) != 5 && particleList.type(i) != 0 && particleList.type(i) != 4) {
			particleList.rigid(i) = 1;
			//particleList.vel(i) = glm::vec3(0,0,0);
			//particleList.predicted_pos(i) = particleList.pos(i)+dt*particleList.predicted_vel(i);
		}

		//s = s + corr_s;
		//Silty sands-compacted and clayey sandy gravels -> 20
		float beta = 20.0f;
		float C = 100000.0f;//or 5000 for minor tilt or 100000 for cube 
		/*glm::vec3 C = glm::vec3(200,200,200);//UPPER LIMIT 200//glm::vec3(58,58,58); //length = 100
		if (glm::length(particleList.stress(i)) <= pow(beta,2.0f)*glm::length(C)) {
			//particleList.stress(i) = pow(beta,2.0f)*C;
			//particleList.rigid(i) = 1;
		}*/
		if (particleList.stress(i) <= pow(beta,2.0f)*C) {
			//particleList.rigid(i) = 1;
			//particleList.predicted_vel(i) = glm::vec3(0,0,0);
		}
	}

}

void WaterSim::assignRigidClusters() {
	for (int i = 0; i < particleList.size(); i++) {
		if (particleList.type(i) != 4) {
			if (particleList.rigid(i) == 1) {//IS RIGID
				if (particleList.in_cluster(i) == 0) {//NOT IN A RIGID CLUSTER
					vector<int> temp;
					for (int j = 0; j < particleList.nbr(i).size(); j++) {//IF A NEIGHBORING PARTICLE IS PART OF A CLUSTER AND RIGID
						if (particleList.in_cluster(particleList.nbr(i)[j]) == 1 && particleList.rigid(particleList.nbr(i)[j]) == 1 && 
							particleList.cluster_num(particleList.nbr(i)[j]) != -1) {//ADD TO LIST OF POSSIBLE CLUSTERS IF CLUSTER NUM IS NOT -1
							temp.push_back(particleList.cluster_num(particleList.nbr(i)[j]));
						}
					}
					if (temp.size() > 1) {//MULTIPLE NEIGHBORING PARTICLES ARE RIGID AND HAVE DEFINED CLUSTERS, SO COMBINE THEM
						int clusterLocation = -1;
						for (int k = 0; k < rigidClusters.size(); k++) {
							if (rigidClusters[k].size() == 0) {
								clusterLocation = k;
								break;
							}
						}
						vector<int> c1;
						c1.push_back(i);
						if (clusterLocation == -1) {
							particleList.cluster_num(i) = rigidClusters.size();
						} else {
							particleList.cluster_num(i) = clusterLocation;
						}
						particleList.in_cluster(i) = 1;
						for (int j = 0; j < temp.size(); j++) {
							//c1.resize(c1.size()+rigidClusters[temp[j]].size());
							c1.insert(c1.end(), rigidClusters[temp[j]].begin(), rigidClusters[temp[j]].end());
							for (int k = 0; k < rigidClusters[temp[j]].size(); k++) {
								if (clusterLocation == -1) {
									particleList.cluster_num(rigidClusters[temp[j]][k]) = rigidClusters.size();
								} else {
									particleList.cluster_num(rigidClusters[temp[j]][k]) = clusterLocation;
								}
							}
							rigidClusters[temp[j]].clear();
						}
						if (clusterLocation == -1) {
							rigidClusters.push_back(c1);
						} else {
							rigidClusters[clusterLocation].clear();
							rigidClusters[clusterLocation].insert(rigidClusters[clusterLocation].end(), c1.begin(), c1.end());
						}
						//rigidClustersForces.push_back(particleList.fricoh(i)/(rigidClusters[i].size()+0.000001f));
					} else if (temp.size() == 1) { //IF THERE IS ONLY 1 NEIGHBORING RIGID CLUSTER PARTICLE WITH A DEFINED CLUSTER, ASSIGN THIS PARTICLE THAT CLUSTER NUM
						particleList.cluster_num(i) = particleList.cluster_num(temp[0]);
						particleList.in_cluster(i) = 1;
					} else {//if (temp.size() == 0) {//particleList.in_cluster(i) != 1) {
						int clusterLocation = -1;
						for (int k = 0; k < rigidClusters.size(); k++) {
							if (rigidClusters[k].size() == 0) {
								clusterLocation = k;
								break;
							}
						}
						vector<int> c; //IF THERE ARE NO NEIGHBORING RIGID PARTICLES WITH DEFINED CLUSTERS MAKE NEW CLUSTER
						c.push_back(i);
						if (clusterLocation == -1) {
							particleList.cluster_num(i) = rigidClusters.size();
						} else {
							particleList.cluster_num(i) = clusterLocation;
						}
						particleList.in_cluster(i) = 1;
						for (int j = 0; j < particleList.nbr(i).size(); j++) {
							if (particleList.in_cluster(particleList.nbr(i)[j]) == 0 && particleList.rigid(particleList.nbr(i)[j]) == 1) {//IF NEIGHBORING PARTICLE IS RIGID BUT HAS NO CLUSTER
								c.push_back(particleList.nbr(i)[j]);
								if (clusterLocation == -1) {
									particleList.cluster_num(particleList.nbr(i)[j]) = rigidClusters.size();
								} else {
									particleList.cluster_num(particleList.nbr(i)[j]) = clusterLocation;
								}
								particleList.in_cluster(particleList.nbr(i)[j]) = 1;
							}
						}
						if (clusterLocation == -1) {
							rigidClusters.push_back(c);
						} else {
							rigidClusters[clusterLocation].clear();
							rigidClusters[clusterLocation].insert(rigidClusters[clusterLocation].end(), c.begin(), c.end());
						}
						//rigidClustersForces.push_back(particleList.fricoh(i)/(rigidClusters[rigidClusters.size()-1].size()+0.000001f));
					}
				} else {//PARTICLE IS RIGID AND IS IN A RIGID CLUSTER
					vector<int> temp;
					for (int j = 0; j < particleList.nbr(i).size(); j++) {//IF A NEIGHBORING PARTICLE IS PART OF A CLUSTER AND RIGID
						if (particleList.in_cluster(particleList.nbr(i)[j]) == 1 && particleList.rigid(particleList.nbr(i)[j]) == 1 &&
							particleList.cluster_num(particleList.nbr(i)[j]) != -1 && 
							particleList.cluster_num(particleList.nbr(i)[j]) == particleList.cluster_num(i)) {
							temp.push_back(particleList.cluster_num(particleList.nbr(i)[j]));
						}
					}
					if (temp.size() > 0) {//IF AT LEAST 1 CLUSTER NUMBER IS IN TEMP WE NEED TO COMBINE ALL OF THEM
						int clusterLocation = -1;
						for (int k = 0; k < rigidClusters.size(); k++) {
							if (rigidClusters[k].size() == 0) {
								clusterLocation = k;
								break;
							}
						}
						vector<int> c1;
						c1.push_back(i);
						if (clusterLocation == -1) {
							particleList.cluster_num(i) = rigidClusters.size();
						} else {
							particleList.cluster_num(i) = clusterLocation;
						}
						for (int j = 0; j < temp.size(); j++) {
							//c1.resize(c1.size()+rigidClusters[temp[j]].size());
							c1.insert(c1.end(), rigidClusters[temp[j]].begin(), rigidClusters[temp[j]].end());
							for (int k = 0; k < rigidClusters[temp[j]].size(); k++) {
								if (clusterLocation == -1) {
									particleList.cluster_num(rigidClusters[temp[j]][k]) = rigidClusters.size();
								} else {
									particleList.cluster_num(rigidClusters[temp[j]][k]) = clusterLocation;
								}
							}
							rigidClusters[temp[j]].clear();
						}
						if (clusterLocation == -1) {
							rigidClusters.push_back(c1);
						} else {
							rigidClusters[clusterLocation].clear();
							rigidClusters[clusterLocation].insert(rigidClusters[clusterLocation].end(), c1.begin(), c1.end());
						}
					}
				}
			}
		}
	}
}

void WaterSim::assignRigidBodyVelocity(float dt) {
	vector<glm::vec3> f;
	vector<glm::vec3> g;
	vector<glm::vec3> v;
	vector<glm::vec3> d;
	vector<glm::vec3> vi;
	vector<glm::vec3> p;
	vector<float> n;
	for (int i = 0; i < rigidClusters.size(); i++) {
		glm::vec3 fricohsum;
		glm::vec3 discsum;
		glm::vec3 viscsum;
		glm::vec3 gravsum;
		glm::vec3 psum;
		glm::vec3 vel;
		int num = 0;
		for (int j = 0; j < rigidClusters[i].size(); j++) {
			fricohsum += particleList.fricoh(rigidClusters[i][j]);
			discsum += particleList.disc(rigidClusters[i][j]);
			viscsum += particleList.visc(rigidClusters[i][j]);
			psum += particleList.press_force(rigidClusters[i][j]);
			gravsum += glm::vec3(0,-9.81,0);
			vel += particleList.predicted_vel(rigidClusters[i][j]);
			num += 1;
		}
		f.push_back(fricohsum);///(rigidClusters[i].size()+0.000001f));
		d.push_back(discsum);///(rigidClusters[i].size()+0.000001f));
		vi.push_back(viscsum/(rigidClusters[i].size()+0.000001f));
		p.push_back(psum/(rigidClusters[i].size()+0.000001f));
		g.push_back(gravsum);
		v.push_back(vel);
		n.push_back(num);
	}
	for (int i = 0; i < rigidClusters.size(); i++) {
		for (int j = 0; j < rigidClusters[i].size(); j++) {
			particleList.predicted_vel(rigidClusters[i][j]) = v[i]/n[i];//particleList.vel(i) + dt * (d[i] + glm::vec3(0,-9.81,0) + f[i]);// + p[i] + vi[i]);//g[i]);// + d[i]);// + vi[i]);//v[i]/n[i];//1.0f/dt * (f[i]-rigidClustersForces[i]+glm::vec3(0,-9.8f,0));
			particleList.predicted_pos(rigidClusters[i][j]) = particleList.pos(rigidClusters[i][j]) + dt * particleList.predicted_vel(rigidClusters[i][j]);
		}
	}
	//rigidClustersForces.clear();
	//rigidClustersForces.insert(rigidClustersForces.end(), f.begin(), f.end());
}

void WaterSim::interpolateBridge(glm::vec3 R, float h, int i, vector<int> nbrs) {//ONLY FOR WET GRANULAR PARTICLES
	glm::vec3 As = glm::vec3(0,0,0);
	glm::vec3 temp;
	float k_bridge = 0.6f;
	float Wf = 1.0f;
	for (int j = 0; j < nbrs.size(); j++) {
		if ((particleList.type(nbrs[j]) == 1 || particleList.type(nbrs[j]) == 5) && particleList.wet(i) > 0 //low or high res granular
			&& glm::dot(particleList.vel(nbrs[j])-particleList.vel(i), particleList.pos(nbrs[j])-particleList.pos(i)) > 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += k_bridge*max(0.0f,Wf-((particleList.wet(i)+particleList.wet(nbrs[j]))/2.0f))*(particleList.vel(nbrs[j])-particleList.vel(i));
		}	
	}
	particleList.bridge(i) = As;
}

void WaterSim::disperseWaterWetness(glm::vec3 R, float h, int i, vector<int> nbrs) {//ONLY FOR WATER PARTICLES
	float count = 0.0f;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4 && particleList.wet(nbrs[j]) < Wmax) {
			count++;
		}	
	}
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4 && particleList.wet(nbrs[j]) < Wmax) {
			particleList.wet(nbrs[j]) += particleList.wet(i)/(count+0.000001f);
		}	
	}
	//particle should then disappear, so in future if an index is in the following vector it cannot be used and will no longer be updated
	//particleList.type(i) = 0;
	//dispersedWaterParticles.push_back(i);
	/*
	if(std::find(dispersedWaterParticles.begin(), dispersedWaterParticles.end(), i) != dispersedWaterParticles.end()) {
    /* v contains x */
	//} else {
    /* v does not contain x */
	//}
}

void WaterSim::disperseGranularWetness(glm::vec3 R, float h, int i, vector<int> nbrs) {//ONLY FOR GRANULAR PARTICLES
	if (particleList.wet(i)>Wthreshold) {
		float count = 0;
		for (int j = 0; j < nbrs.size(); j++) {
			if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4 && particleList.wet(nbrs[j]) < Wmax) {
				count++;
			}	
		}
		for (int j = 0; j < nbrs.size(); j++) {
			if (particleList.type(nbrs[j]) != 0 && particleList.type(nbrs[j]) != 4 && particleList.wet(nbrs[j]) < Wmax) {
				particleList.wet(nbrs[j]) += (particleList.wet(i)-Wthreshold)/(count+0.000001f);
			}	
		}
		particleList.wet(i) = Wthreshold;
	}
}

//Current one we are using for granular flow
//The k<3 loop is the one that is changing the simulation 
void WaterSim::update(const Scene* const scene, float dt) {
	//For all particles
	populateGrid(); 
	for (int i = 0; i < particleList.size(); i++) {
		//Find neighboring particles
		particleList.nbr(i) = neighborSearch(i);
		//particleList.rigid(i) = 0;
	}
	//For all particles
	for (int i = 0; i < particleList.size(); i++) {
		//if (particleList.rigid(i) == 0) {
			//Calculate gravity, external forces, and viscosity
			particleList.grav(i) = glm::vec3(0,-9.81,0);
			//interpolateViscosity(particleList.predicted_pos(i), nbrRadius,i,particleList.nbr(i));
		
			//Initialize pressure and dissipative stress
			particleList.press(i) = 0.0f;//glm::vec3(0,0,0);
			particleList.stress(i) = 0.0f;//glm::vec3(0,0,0);
			particleList.press_force(i) = glm::vec3(0,0,0);
			particleList.fricoh(i) = glm::vec3(0,0,0);
		//}
	} 

	//float stressTensor = calculateStressTensor(dt);
	//float Ferr = calculateFerr();
	//while (Ferr > tolP || stressTensor > tolS) {
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < particleList.size(); i++) {
			//if (particleList.rigid(i) == 0) {
				//Add pressure and dissipative forces (friction and cohesion)
			if (particleList.type(i) != 4 && particleList.type(i) != 0) {
				//Predict velocity and apply forces Vi <- Vi + dt*fext(Xi) 
				particleList.predicted_vel(i) = particleList.vel(i) + dt * (particleList.grav(i) + particleList.disc(i) + particleList.bridge(i));
			} else {
				particleList.predicted_vel(i) = particleList.vel(i) + dt * (particleList.grav(i));
			}
				//Predict position Xi* <- Xi + dt*Vi
				particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.predicted_vel(i);
			//}
		}

		for (int i = 0; i < particleList.size(); i++) {
			//if (particleList.rigid(i) == 0) {
				//Predict Density
				interpolateDensity(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

				if (particleList.rho(i) > rhoMax) {
					//PRESSURE
					//Compute corrective pressure
					calculateCorrectivePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);
					//Increment pressure
					particleList.press(i) = particleList.press(i) + particleList.corrective_press(i);
					if (particleList.type(i) != 0 && particleList.type(i) != 4) {
						//STRESS
						//Predict strain rate
						calculateStrainRate(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);
						//Compute corrective dissipative stress
						calculateCorrectiveStress(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);
						//Increment stress
						particleList.stress(i) = particleList.stress(i) + particleList.corrective_stress(i);

						testYieldandCohesion(i, dt);
					}
				} else {
					if (particleList.type(i) != 4 && particleList.type(i) != 0) {
						//Add discrete particle forces to close particles
						calculateDiscreteParticleForces(i, particleList.nbr(i));
					}
				}
			//}
		}

		for (int i = 0; i < granularParticleNum; i++) {
			//interpolatePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
			//interpolateFriCoh(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
			interpolateBridge(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
		}
	}
	//PUT BACK
	assignRigidClusters();
	assignRigidBodyVelocity(dt); //updates rigid body velocity and position
	resolve_constraints(scene);
	float Kp = 0.5f;
	for (int i = 0; i < particleList.size(); i++) {
		//Update velocity
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));
		//Apply viscosity
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		collision_detection(i, scene);
		particleList.pos(i) = particleList.predicted_pos(i);
		
		disperseWaterWetness(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
		disperseGranularWetness(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

		//Assign new particle wetness value
		//Kp is speed of propagation
		//Larger Kp will yield faster propagation of wetness while lower Kp will yield lower propagation of wetness
		if (particleList.type(i) != 4 && particleList.type(i) != 0) {
			float deltaW = particleList.wet(i)-Wthreshold;
			if (deltaW < 0) {
				deltaW = 0;
			}
			particleList.wet(i) = particleList.wet(i) + Kp*((deltaW)/(particleList.nbr(i).size()+0.000001f))*dt;
		}
	}

	for (int i = 1000; i < particleList.size()-1000; i++) {
		//HIGH RES PARTICLE CALCULATIONS
		int count = 0;
		for (int j = 0; j < particleList.nbr(i).size(); j++) {
			if (particleList.type(particleList.nbr(i)[j]) == 1) {
				count++;
			}
		}
		if (count > 0) {
			hrCalcVelocity(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);
			hrCalcPosition(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);
		}
	}
	
	updateSpheres();
}

void WaterSim::updateSpheres() {
	for (int i = 0; i < dots.size()/*-3*/; i++) {
		dots[i]->m_center = particleList.pos(i);
		dots[i]->mat_color = glm::vec3(abs(particleList.vel(i).x/glm::length(particleList.vel(i))), abs(particleList.vel(i).y/glm::length(particleList.vel(i))), abs(particleList.vel(i).z/glm::length(particleList.vel(i))));
		if (particleList.type(i) == 4) {
			dots[i]->mat_color = glm::vec3(0,0,1);
		} else {
			dots[i]->mat_color = glm::vec3(sqrt(pow(particleList.wet(i),2.0f)), sqrt(pow(particleList.wet(i),2.0f)), sqrt(pow(particleList.wet(i),2.0f)));
		}
		dots[i]->update();
	}
}

void WaterSim::resolve_constraints(const Scene* const scene) {
	for (unsigned int iter = 0; iter < solverIterations; ++iter) {
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Calculate lambda
			//MOST TIME DRAIN HERE
			if (particleList.type(i) != 5) {
				particleList.lambda(i) = calculateLambda(i); 
			}
		}
		//For all particles
		for (int i = 0; i < particleList.size(); i++) {
			//Calculate deltaP
			if (particleList.type(i) != 5) {
				particleList.deltaP(i) = calculateDeltaP(i); 
			}

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
			if (particleList.type(i) != 5) {
				particleList.predicted_pos(i) = particleList.predicted_pos(i) + particleList.deltaP(i);
			}
			//collision_detection(i, scene);
		}
	}
}

float WaterSim::kernelViscositySecondDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	return (45.0f/(PI*pow(h,6.0f)+0.000001f))*(h-r);
}

glm::vec3 WaterSim::kernelPressureDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	float temp = (45.0f/((PI*pow(h,6.0f)*r)+0.000001f)) * pow((h-r),2.0f);
//	float temp = (45.0f/((PI*pow(h,6.0f)*r)+0.000001f)) * pow((h-r),3.0f);
	return -R * temp; //temp*(R/(glm::length(R)+0.000001f));//
	//-r should be the actual R)
}

glm::vec3 WaterSim::kernelGeneralDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	//Issue here
	float temp = (945.0f/((32.0f*PI*pow(h,9.0f))+0.000001f)) * pow((pow(h,2.0f)-pow(r,2.0f)),2.0f);
	
	temp = ((945.0f/((32.0f*PI*pow(h,8.0f))+0.000001f)) * pow((pow(h,2.0f)-pow(r,2.0f)),2.0f))/* -
		((2835.0f*pow((pow(h,2.0f)-pow(r,2.0f)),3.0f))/(64.0f*PI*pow(h,10.0f)+0.000001f))*/;
	return -R * temp;
}

float WaterSim::kernelGeneralSecondDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	return (945/((8*PI*pow(h,9))+0.000001f)) * (pow(h,2) - pow(r,2)) * (pow(r,2) - ((3.f/4.f) * (pow(h,2) - pow(r,2))));
}

void WaterSim::applyXSPH(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	float norm = 0;
	float c = 0.01f;//0.1f;//0.1f;//0.01;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			norm = glm::length(temp);//sqrt(pow(temp.x,2)+pow(temp.y,2)+pow(temp.z,2));
			//As += (particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscosity(norm, h); 
			if ((particleList.type(nbrs[j]) != 4 && particleList.type(i) != 4) || (particleList.type(nbrs[j]) == 4 && particleList.type(i) == 4)) {
				As += (particleList.mass(j)/(0.000001f + particleList.rho(j)))*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelGeneral(norm,h);//SecondDeriv(temp, h); 
			} else {
				if (particleList.type(i) == 4) {
					As += (particleList.mass(i)/(particleList.rho(i)+0.000001f))*(particleList.mass(nbrs[j])/(particleList.rho(nbrs[j])+0.000001f))/(2.0f)*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelGeneral(norm, h);
				} else {
					As += -1.0f * (particleList.mass(i)/(particleList.rho(i)+0.000001f))*(particleList.mass(nbrs[j])/(particleList.rho(nbrs[j])+0.000001f))/(2.0f)*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelGeneral(norm, h);
				}
			}
		}
	}
	particleList.vel(i) = particleList.vel(i) + c*As;
}

float WaterSim::hrCalcWeight(float r) {
	return max(0.0f, pow((1.0f-(pow(r,2.0f)/pow(3.0f*particleRad,2.0f))),3.0f));
}

glm::vec3 WaterSim::hrCalcAvgVelocity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	float denom = 0.0f;
	glm::vec3 As;
	glm::vec3 temp;
	float w = 0.0f;
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) == 1) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			w = hrCalcWeight(glm::length(temp));
			denom += w;
			As += w*particleList.vel(nbrs[j]);
		}
	}
	return (1.0f/(denom + 0.000001f))*As;
}

void WaterSim::hrCalcVelocity(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	float alpha = hrCalcAlpha(R,h,i,nbrs);
	particleList.vel(i) = (1.0f-alpha)*hrCalcAvgVelocity(R,h,i,nbrs)+alpha*(particleList.vel(i)+dt*(glm::vec3(0,-9.81,0)/particleList.mass(i)));
}

float WaterSim::hrCalcAlpha(glm::vec3 R, float h, int i, vector<int> nbrs) {
	float argmax = FLT_MIN;
	float w = 0.0f;
	float sum = 0.0f;
	glm::vec3 temp = glm::vec3(0,0,0);
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) == 1) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			w = hrCalcWeight(glm::length(temp));
			if (w > argmax) {
				argmax = w;
			}
			sum +=w;
		}
	}

	float wrLR = hrCalcWeight(particleRad);

	if (argmax/(sum+0.000001f) >= 0.6f) {
		return 1.0f - argmax;
	} else if (argmax <= wrLR) {
		return 1.0f - argmax;
	} else {
		return 0;
	}
}

void WaterSim::hrCalcPosition(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	particleList.pos(i) = particleList.pos(i) + dt*particleList.vel(i);
}
