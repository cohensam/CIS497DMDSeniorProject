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

	tolP = 0.000001f;
	tolS = 0.000001f;
	rhoMax = 2000.0f;
	//Degrees, for dry loose beach sand
	angleRepose = 38.0f;

	particleNum = 1000;//8000;//2197;//1728;//1000;//4096;//
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

	if ((particleNum == 8000 && particleRad == 0.1f) || particleNum == 8000+3375)
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
	}

	nbrRadius = 3*particleRad;
	/*gridXmin = xmin;
	gridXmax = xmax;
	gridYmin = ymin;
	gridYmax = ymax;
	gridZmin = zmin;
	gridZmax = zmax; */
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

	xy = (xmax-xmin)*(ymax-ymin);
	xz = (xmax-xmin)*(zmax-zmin);
	yz = (ymax-ymin)*(zmax-zmin);
	boundaryParticleNum = 2*xy+2*xz+2*yz;

	particleList.resize(particleNum+boundaryParticleNum);

	// density of water is 1000 kg/m^3 = 1 g/cm^3
	// particle has volume: (4/3)pi*0.1^3 = 0.004189
	// mass = density * volume = 4.189

	// initialize triangle mesh limits
	//tri_xmax = 0.0f; 
	//tri_xmin = 0.0f; 
	//tri_ymax = 0.0f; 
	//tri_ymin = 0.0f; 
	//tri_zmax = 0.0f; 
	//tri_zmin = 0.0f; 



		 // These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

// a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
/*a2iEdgeConnection[0][0] = 0;
a2iEdgeConnection[0][1] = 1;
a2iEdgeConnection[1][0] = 1;
a2iEdgeConnection[1][1] = 2;
a2iEdgeConnection[2][0] = 2;
a2iEdgeConnection[2][1] = 3;
a2iEdgeConnection[3][0] = 3;
a2iEdgeConnection[3][1] = 0;
a2iEdgeConnection[4][0] = 4;
a2iEdgeConnection[4][1] = 5;
a2iEdgeConnection[5][0] = 5;
a2iEdgeConnection[5][1] = 6;
a2iEdgeConnection[6][0] = 6;
a2iEdgeConnection[6][1] = 7;
a2iEdgeConnection[7][0] = 7;
a2iEdgeConnection[7][1] = 4;
a2iEdgeConnection[8][0] = 0;
a2iEdgeConnection[8][1] = 4;
a2iEdgeConnection[9][0] = 1;
a2iEdgeConnection[9][1] = 5;
a2iEdgeConnection[10][0] = 2;
a2iEdgeConnection[10][1] = 6;
a2iEdgeConnection[11][0] = 3;
a2iEdgeConnection[11][1] = 7;*/
	/*[12][2] = 
{
    {0,1}, {1,2}, {2,3}, {3,0},
    {4,5}, {5,6}, {6,7}, {7,4},
    {0,4}, {1,5}, {2,6}, {3,7}
};*/

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
/*float a2fEdgeDirection[12][3] =
{
    {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
    {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
    {0.0, 0.0, 1.0},{0.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{0.0,  0.0, 1.0}
};

// aiCubeEdgeFlags lists flags for all possible cases
// each flag has 12 bits, corresponding to 12 edges. 0 for no intersection, 1 for intersection.
static int aiCubeEdgeFlags[256]=
{
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460, 
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0, 
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230, 
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

// a2iTriangleConnectionTable lists for each case, how many triangles are there and which edges they resides on.
// -1 for no triangle any more, 0 ~ 11 indicating the edge id.
// for any case, up to 5 triangles are possible, and we need another value in the end to indicate the table item is done.
// so each item has 5 * 3 + 1 = 16 elements.
a2iTriangleConnectionTable[0][0] = -1;
a2iTriangleConnectionTable[0][1] = -1;
a2iTriangleConnectionTable[0][2] = -1;
a2iTriangleConnectionTable[0][3] = -1;
a2iTriangleConnectionTable[0][4] = -1; 
a2iTriangleConnectionTable[0][5] = -1;
a2iTriangleConnectionTable[0][6] = -1;
a2iTriangleConnectionTable[0][7] = -1;
a2iTriangleConnectionTable[0][8] = -1;
a2iTriangleConnectionTable[0][9] = -1;
a2iTriangleConnectionTable[0][10] = -1; 
a2iTriangleConnectionTable[0][11] = -1; 
a2iTriangleConnectionTable[0][12] = -1; 
a2iTriangleConnectionTable[0][13] = -1; 
a2iTriangleConnectionTable[0][14] = -1; 
a2iTriangleConnectionTable[0][15] = -1;

a2iTriangleConnectionTable[1][0] = 0;
a2iTriangleConnectionTable[1][1] = 8;
a2iTriangleConnectionTable[1][2] = 3;
a2iTriangleConnectionTable[1][3] = -1;
a2iTriangleConnectionTable[1][4] = -1;
a2iTriangleConnectionTable[1][5] = -1;
a2iTriangleConnectionTable[1][6] = -1;
a2iTriangleConnectionTable[1][7] = -1;
a2iTriangleConnectionTable[1][8] = -1;
a2iTriangleConnectionTable[1][9] = -1;
a2iTriangleConnectionTable[1][10] = -1;
a2iTriangleConnectionTable[1][11] = -1;
a2iTriangleConnectionTable[1][12] = -1;
a2iTriangleConnectionTable[1][13] = -1;
a2iTriangleConnectionTable[1][14] = -1;
a2iTriangleConnectionTable[1][15] = -1;

a2iTriangleConnectionTable[2][0] = 0;
a2iTriangleConnectionTable[2][1] = 1;
a2iTriangleConnectionTable[2][2] = 9;
a2iTriangleConnectionTable[2][3] = -1;
a2iTriangleConnectionTable[2][4] = -1;
a2iTriangleConnectionTable[2][5] = -1;
a2iTriangleConnectionTable[2][6] = -1;
a2iTriangleConnectionTable[2][7] = -1;
a2iTriangleConnectionTable[2][8] = -1;
a2iTriangleConnectionTable[2][9] = -1;
a2iTriangleConnectionTable[2][10] = -1;
a2iTriangleConnectionTable[2][11] = -1;
a2iTriangleConnectionTable[2][12] = -1;
a2iTriangleConnectionTable[2][13] = -1;
a2iTriangleConnectionTable[2][14] = -1;
a2iTriangleConnectionTable[2][15] = -1;

a2iTriangleConnectionTable[3][0] = 1;
a2iTriangleConnectionTable[3][1] = 8;
a2iTriangleConnectionTable[3][2] = 3;
a2iTriangleConnectionTable[3][3] = 9;
a2iTriangleConnectionTable[3][4] = 8;
a2iTriangleConnectionTable[3][5] = -1;
a2iTriangleConnectionTable[3][6] = -1;
a2iTriangleConnectionTable[3][7] = -1;
a2iTriangleConnectionTable[3][8] = -1;
a2iTriangleConnectionTable[3][9] = -1;
a2iTriangleConnectionTable[3][10] = -1;
a2iTriangleConnectionTable[3][11] = -1;
a2iTriangleConnectionTable[3][12] = -1;
a2iTriangleConnectionTable[3][13] = -1;

a2iTriangleConnectionTable[4][0] = 1;
a2iTriangleConnectionTable[4][1] = 2;
a2iTriangleConnectionTable[4][2] = 10;
a2iTriangleConnectionTable[4][3] = -1;
a2iTriangleConnectionTable[4][4] = -1;
a2iTriangleConnectionTable[4][5] = -1;
a2iTriangleConnectionTable[4][6] = -1;
a2iTriangleConnectionTable[4][7] = -1;
a2iTriangleConnectionTable[4][8] = -1;
a2iTriangleConnectionTable[4][9] = -1;
a2iTriangleConnectionTable[4][10] = -1;
a2iTriangleConnectionTable[4][11] = -1;
a2iTriangleConnectionTable[4][12] = -1;
a2iTriangleConnectionTable[4][13] = -1;
a2iTriangleConnectionTable[4][14] = -1;
a2iTriangleConnectionTable[4][15] = -1;
a2iTriangleConnectionTable[4][14] = -1;
a2iTriangleConnectionTable[4][15] = -1;

a2iTriangleConnectionTable[5][0] = 0;
a2iTriangleConnectionTable[5][1] = 8;
a2iTriangleConnectionTable[5][2] = 3;
a2iTriangleConnectionTable[5][3] = 1;
a2iTriangleConnectionTable[5][4] = 2;
a2iTriangleConnectionTable[5][5] = 10;
a2iTriangleConnectionTable[5][6] = -1;
a2iTriangleConnectionTable[5][7] = -1;
a2iTriangleConnectionTable[5][8] = -1;
a2iTriangleConnectionTable[5][9] = -1;
a2iTriangleConnectionTable[5][10] = -1;
a2iTriangleConnectionTable[5][11] = -1;
a2iTriangleConnectionTable[5][12] = -1;
a2iTriangleConnectionTable[5][13] = -1;
a2iTriangleConnectionTable[5][14] = -1;
a2iTriangleConnectionTable[5][15] = -1;

a2iTriangleConnectionTable[6][0] = 9;
a2iTriangleConnectionTable[6][1] = 2;
a2iTriangleConnectionTable[6][2] = 10;
a2iTriangleConnectionTable[6][3] = 0;
a2iTriangleConnectionTable[6][4] = 2;
a2iTriangleConnectionTable[6][5] = 9;
a2iTriangleConnectionTable[6][6] = -1;
a2iTriangleConnectionTable[6][7] = -1;
a2iTriangleConnectionTable[6][8] = -1;
a2iTriangleConnectionTable[6][9] = -1;
a2iTriangleConnectionTable[6][10] = -1;
a2iTriangleConnectionTable[6][11] = -1;
a2iTriangleConnectionTable[6][12] = -1;
a2iTriangleConnectionTable[6][13] = -1;
a2iTriangleConnectionTable[6][14] = -1;
a2iTriangleConnectionTable[6][15] = -1;

a2iTriangleConnectionTable[7][0] = 2;
a2iTriangleConnectionTable[7][1] = 8;
a2iTriangleConnectionTable[7][2] = 3;
a2iTriangleConnectionTable[7][3] = 2;
a2iTriangleConnectionTable[7][4] = 10;
a2iTriangleConnectionTable[7][5] = 8;
a2iTriangleConnectionTable[7][6] = 10;
a2iTriangleConnectionTable[7][7] = 9;
a2iTriangleConnectionTable[7][8] = 8;
a2iTriangleConnectionTable[7][9] = -1;
a2iTriangleConnectionTable[7][10] = -1;
a2iTriangleConnectionTable[7][11] = -1;
a2iTriangleConnectionTable[7][12] = -1;
a2iTriangleConnectionTable[7][13] = -1;
a2iTriangleConnectionTable[7][14] = -1;
a2iTriangleConnectionTable[7][15] = -1;*/
	/*[256][16] =  
{
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};*/
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
		// density of water is 1000 kg/m^3 = 1 g/cm^3
		// particle has volume: (4/3)pi*0.1^3 = 0.004189
		// mass = density * volume = 4.189
		particleList.set_mass(i, rho_0*(4.0f/3.0f)*PI*pow(particleRad,3));
		//particleList.set_mass(i, (4.0f/3.0f)*PI*pow(particleRad,3));
		particleList.type(i) = 1;
    }

	//Set original position of each particle to be a cube of 20x20x20
	int pos = 0;
	float lim = 2*10*particleRad; 
	if (particleNum == 8000 || particleNum == 8000+3375) 
		lim = 4*10*particleRad; 
	for (float i = 0; i < /*4*//*3.2*/lim; i+=(2*particleRad)) 
	{
		for (float j = 0; j < /*4*//*3.2*/lim; j+=(2*particleRad))
		{
			for (float k = 0; k < /*4*//*3.2*/lim; k+=(2*particleRad))
			{
				particleList.pos(pos) = glm::vec3(i-0.5f,j+4.15f,k-0.5f);
				pos++;
			}
		}
	}
	//lim = 3*10*particleRad;  
	//for (float i = 0; i < /*4*//*3.2*/lim; i+=(2*particleRad)) 
	//{
	//	for (float j = 0; j < /*4*//*3.2*/lim; j+=(2*particleRad))
	//	{
	//		for (float k = 0; k < /*4*//*3.2*/lim; k+=(2*particleRad))
	//		{
	//			particleList.pos(pos) = glm::vec3(i-1.15f,j+0.15f,k-0.15f);
	//			pos++;
	//		}
	//	}
	//}

	//Set up boundary particles, type=0
	int id = particleNum; 
	for (int x = xmin; x < xmax; x++) {
		for (int y = ymin; y < ymax; y++) {
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(x,y,zmin);
			id++;
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(x,y,zmax);
			id++;
		}
	}

	for (int x = xmin; x < xmax; x++) {
		for (int z = zmin; z < zmax; z++) {
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(x,ymin,z);
			id++;
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(x,ymax,z);
			id++;
		}
	}

	for (int y = ymin; y < ymax; y++) {
		for (int z = zmin; z < zmax; z++) {
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(xmin,y,z);
			id++;
			particleList.type(id) = 0;
			particleList.set_mass(id,1);
			particleList.pos(id) = glm::vec3(xmax,y,z);
			id++;
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
    glm::vec3 x, p, q, n, v;
    x = particleList.pos(id);
    p = particleList.predicted_pos(id);
	v = particleList.vel(id); 
	float rest = 1.0f; 
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
		if (n[0] != 0)
			particleList.vel(id) = glm::vec3(-rest*v[0],v[1],v[2]);
		else if (n[1] != 0) 
			particleList.vel(id) = glm::vec3(v[0],-rest*v[1],v[2]);
		else if (n[2] != 0) 
			particleList.vel(id) = glm::vec3(v[0],v[1],-rest*v[2]);
    }
}*/

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

/*void WaterSim::resolve_collisions()
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
}*/

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

bool WaterSim::triangle_collision(/*const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& v, */int tri, int id) {
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
	glm::vec3 tempy = glm::vec3(0.0f);
	dots.resize(particleNum);
	for (int i = 0; i < particleNum; i++) {
		
		tempy = glm::vec3(abs(particleList.vel(i).x/glm::length(particleList.vel(i))), abs(particleList.vel(i).y/glm::length(particleList.vel(i))), abs(particleList.vel(i).z/glm::length(particleList.vel(i))));
		//glm::vec3 tempy = glm::vec3(0.0f, 0.0f, 1.0f); 
		s = new Sphere(particleList.pos(i),particleRad, tempy);
		dots[i] = s;//.push_back(s);
	}
}

void WaterSim::populateGrid() {
	/*grid.resize(133);
	for (int i = 0; i < 133; i++) {
		grid[i].resize(133);
		for(int j = 0; j < 133; j++) {
			grid[i][j].resize(133);
		}
	}*/
	
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
	float rho_i = 0.0; 
	float Wj = 0.0; 
	std::vector<float> W; 
	int nbrid = 0;
	float dist2 = 0.0;
	float mj = 0.0; 
	rho_i = interpolateDensity(particleList.predicted_pos(id), nbrRadius,id,particleList.nbr(id));

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
			As += particleList.mass(nbrs[j])*((particleList.press(i)+particleList.press(nbrs[j]))/(2.0f*particleList.rho(nbrs[j])+0.000001f)) * kernelPressureDeriv(temp,h);
		}	
	}
	//for water only we did:
	particleList.press(i) = -As;
	//particleList.press(i) = -particleList.mass(i)*As;
}

void WaterSim::interpolateViscosity(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As = glm::vec3(0,0,0);
	glm::vec3 temp;
	//mu between 0 and 1
	float mu = 0.6f;//1.0f; 
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j])*((particleList.vel(nbrs[j])-particleList.vel(i))/(particleList.rho(nbrs[j])+0.000001f))*kernelViscositySecondDeriv(temp, h); 
		}
	}
		particleList.visc(i) = mu*As;
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
		if (particleList.type(nbrs[j]) != 0) {
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
	}
	
	particleList.vort(i) = fvort; 
	particleList.vel(i) += particleList.vort(i);
}

void WaterSim::interpolateFriCoh(glm::vec3 R, float h, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j]) * (particleList.stress(i)/(pow(particleList.rho(i),2.0f)+0.000001f) + particleList.stress(j)/(pow(particleList.rho(j),2.0f)+0.000001f)) * kernelGeneralDeriv(temp,h);
		}
	}
	particleList.fricoh(i) = As;
}

void WaterSim::calculateCorrectivePressure(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	//You need to check to see if we are using the general or the pressure kernel
	glm::vec3 As1;
	glm::vec3 As2;
	glm::vec3 temp;
	glm::vec3 kernelCalc;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			kernelCalc = kernelPressureDeriv(temp,h);
			As1 += kernelCalc;
			As2 += kernelCalc*kernelCalc;
		}
	}
	particleList.corrective_press(i) = pow(particleList.rho(i),2.0f)/((2.0f*pow(particleList.mass(i),2.0f)+0.000001f)*pow(dt,2.0f)*-As1*As1-As2);
}

void WaterSim::calculateDiscreteParticleForces(int i, vector<int> nbrs) {
	float R1 = particleRad;
	float R2 = particleRad;
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
	float tc = 0.01f;
	float Kd = 0;
	float Kr = 0;
	//Coefficient of friction for sand, number for general use 0.6
	float mu = 0.6f;
	glm::vec3 Fn = glm::vec3(0,0,0);
	glm::vec3 Ft = glm::vec3(0,0,0);
	glm::vec3 F = glm::vec3(0,0,0);
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			//Calulate Fn
			X2 = particleList.predicted_pos(nbrs[j]);
			V2 = particleList.predicted_vel(nbrs[j]);
			zeta = max(0.0f,R1+R2-(X1-X2).length());
			N = X1-X2;
			if (N.length() >= 0) {
				N /= N.length();
			} else {
				N /= 0.000001f;
			}
			V = V1-V2;
			zetaDot = V * N;
			Vt = V-zetaDot*N;
			Meff = particleList.mass(i)*particleList.mass(nbrs[j])/(particleList.mass(i)+particleList.mass(nbrs[j])+0.000001f);
			Kd = 2.0f*Meff*(-log(e)/tc);
			Kr = (Meff/pow(tc,2.0f))*(pow(log(e),2.0f)+pow(3.14159f,2.0f));
			fn = -Kd*pow(zeta,0.5f)*zetaDot-Kr*pow(zeta,3.0f/2.0f);
			Fn = fn*N;

			//Calculate Ft
			if (Vt.length() > 0) {
				Vt /= Vt.length();
			} else {
				Vt /= 0.000001f;
			}
			Ft = -mu*fn*Vt;

			//Calculate F
			F += Fn + Ft;
		}
	}
	
	particleList.disc(i) = F;
}

glm::vec3 WaterSim::interpolateU(glm::vec3 R, float h, float dt, int i, vector<int> nbrs) {
	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j])/(particleList.rho(j)+0.000001f) * kernelGeneralDeriv(temp,h) * particleList.predicted_vel(j);
		}
	}
	return As;
}

void WaterSim::calculateStrainRate(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	//Check general vs. pressure kernel
	glm::vec3 As;
	glm::vec3 temp;
	for(int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			As += particleList.mass(nbrs[j])/(particleList.rho(j)+0.000001f) * kernelGeneralDeriv(temp,h) * particleList.predicted_vel(j);
		}
	}
	particleList.strain(i) = -(As - interpolateU(R,h, dt,i,nbrs));
}

void WaterSim::calculateCorrectiveStress(glm::vec3 R, float h, int i, vector<int> nbrs, float dt) {
	//You need to check the whole transpose thing-is that just negative in these cases or can I ignore it?
	glm::vec3 D;
	glm::vec3 temp;
	glm::vec3 kernelCalc;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			kernelCalc = kernelGeneralDeriv(temp,h);
			D += 1.0f/(particleList.rho(j)+0.000001f)*kernelCalc*kernelCalc;
		}
	}
	D = 2.0f*pow(particleList.mass(i),2.0f)*dt/(pow(particleList.rho(i),2.0f)+0.000001f) * D;
	particleList.corrective_stress(i) = -D*particleList.strain(i);
}

void WaterSim::testYieldandCohesion(int i) {
	/*I do not understand how any of this actually is meant to be used
	to change the stress value*/
	float alpha = 0.5f;//sqrt(2.0f/3.0f)*sin(angleRepose);
	float s = particleList.stress(i).length();
	float corr_s = particleList.corrective_stress(i).length();
	//CHECK IF PRESSURE IS ONLY IN THE Y DIRECTION
	if (corr_s > sqrt(3.0f)*alpha*particleList.press(i).length()) {

		//calculateCorrectiveStress(R,h,i,nbrs,dt);
		corr_s = particleList.corrective_stress(i).length();
	}
	//s = s + corr_s;
	//Silty sands-compacted and clayey sandy gravels -> 20
	float beta = 20.0f;
	float C = 100.0f;//or 5000 for minor tilt or 100000 for cube 
	if (s > pow(beta,2.0f)*C) {
		s = pow(beta,2.0f)*C;
		//s = s + corr_s;
	}

}

float WaterSim::calculateFerr() {
	float maxi = -(1.1754943351e-38);
	float density = 0;
	for (int i = 0; i < particleList.size(); i++) {
		density = max(particleList.rho(i)-rhoMax,0.0f);
		if (density > maxi) {
			maxi = density;
		}
	}
	return maxi;
}

float WaterSim::calculateStressTensor(float dt) {
	float maxi = -(1.1754943351e-38);
	for (int i = 0; i < particleList.size(); i++) {
		glm::vec3 D;
		glm::vec3 temp;
		glm::vec3 kernelCalc;
		glm::vec3 sc;
		glm::vec3 s;
		float As = 0;
		float stress = 0;
		for (int j = 0; j < particleList.nbr(i).size(); j++) {
			temp = particleList.predicted_pos(i)-particleList.predicted_pos(particleList.nbr(i)[j]);
			kernelCalc = kernelGeneralDeriv(temp,nbrRadius);
			D = 2.0f*pow(particleList.mass(i),2.0f)*dt/(pow(particleList.rho(i),2.0f)+0.000001f) * 1.0f/(particleList.rho(j)+0.000001f)*kernelCalc*kernelCalc;
			
			sc = -D*particleList.strain(i);
			s = particleList.stress(i) + sc;
			As += s.x*s.x+s.y*s.y+s.z*s.z;
		}
		//Should this be done in this loop or the outer loop???
		stress = sqrt(As);
		if (stress > maxi) {
			maxi = stress;
		}
	}
	return maxi;
}

//Current one we are using for granular flow
//The k<3 loop is the one that is changing the simulation 
void WaterSim::update(const Scene* const scene, float dt) {
	//For all particles
	populateGrid(); 
	for (int i = 0; i < particleList.size(); i++) {
		//Find neighboring particles
		particleList.nbr(i) = neighborSearch(i);
	}
	//For all particles
	for (int i = 0; i < particleNum; i++) {
		//Calculate gravity, external forces, and viscosity
		particleList.grav(i) = glm::vec3(0,-9.81,0);
	//	interpolateViscosity(particleList.predicted_pos(i), nbrRadius,i,particleList.nbr(i));
		
		//Initialize pressure and dissipative stress
		particleList.press(i) = glm::vec3(0,0,0);
		particleList.stress(i) = glm::vec3(0,0,0);
	} 

	//float stressTensor = calculateStressTensor(dt);
	//float Ferr = calculateFerr();
	//while (Ferr > tolP || stressTensor > tolS) {
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < particleNum; i++) {
			//Add pressure and dissipative forces (friction and cohesion)
			//PUT BACK
		//	interpolatePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
		//	interpolateFriCoh(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

			//Predict velocity and apply forces Vi <- Vi + dt*fext(Xi) 
			particleList.vel(i) = particleList.vel(i) + dt * (particleList.grav(i) + particleList.fricoh(i));// + particleList.press(i));// + particleList.visc(i) );
			
			//Predict position Xi* <- Xi + dt*Vi
			particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);

			//Apply boundary conditions
			//collision_detection();
		}

		for (int i = 0; i < particleNum; i++) {
			//Predict Density
			interpolateDensity(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

			if (particleList.rho(i) > rhoMax) {
				//Compute corrective pressure
				calculateCorrectivePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

				//Increment pressure
				particleList.press(i) = particleList.press(i) + particleList.corrective_press(i);
			} else {
				//Add discrete particle forces to close particles
				calculateDiscreteParticleForces(i, particleList.nbr(i));
				/*particleList.predicted_vel(i) = particleList.vel(i) + dt * (particleList.grav(i) + particleList.visc(i) + 
					particleList.press(i) + particleList.fricoh(i) + particleList.disc(i));
				particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);*/
			}
		}

		for (int i = 0; i < particleNum; i++) {
			//Predict strain rate
			calculateStrainRate(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

			//Compute corrective dissipative stress
			calculateCorrectiveStress(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

			//Increment stress
			particleList.stress(i) = particleList.stress(i) + particleList.corrective_stress(i);

			//Test yield and cohesion
			//testYieldandCohesion(i);
		}	
		//stressTensor = calculateStressTensor(dt);
		//Ferr = calculateFerr();
	}
	//PUT BACK
	resolve_constraints(scene);
	for (int i = 0; i < particleNum; i++) {
		//Update velocity
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));//particleList.predicted_vel(i);
		//visc temp
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		//Update position (and check for collision - not in original algorithm but I want to do it here)
		//particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);
		collision_detection(i, scene);
		particleList.pos(i) = particleList.predicted_pos(i);
	}
	updateSpheres();
}

void WaterSim::updateGranular(const Scene* const scene, float dt) {
	//For all particles
	populateGrid(); 
	for (int i = 0; i < particleList.size(); i++) {
		//Find neighboring particles
		particleList.nbr(i) = neighborSearch(i);
	}
	//For all particles
	for (int i = 0; i < particleList.size(); i++) {
		//Calculate gravity, external forces, and viscosity
		particleList.grav(i) = glm::vec3(0,-9.81,0);
	//	interpolateViscosity(particleList.predicted_pos(i), nbrRadius,i,particleList.nbr(i));
		
		//Initialize pressure and dissipative stress
	//	particleList.press(i) = glm::vec3(0,0,0);
		particleList.stress(i) = glm::vec3(0,0,0);
	} 

	//float stressTensor = calculateStressTensor(dt);
	//float Ferr = calculateFerr();
	//while (Ferr > tolP || stressTensor > tolS) {
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < particleList.size(); i++) {
			//Add pressure and dissipative forces (friction and cohesion)
			//PUT BACK
		//	interpolatePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
			interpolateFriCoh(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

			//Predict velocity and apply forces Vi <- Vi + dt*fext(Xi) 
			particleList.vel(i) = particleList.vel(i) + dt * (particleList.grav(i) + particleList.fricoh(i));// + particleList.press(i));// + particleList.visc(i) );
			
			//Predict position Xi* <- Xi + dt*Vi
			particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);

			//Apply boundary conditions
			//collision_detection();
		}

		/*for (int i = 0; i < particleList.size(); i++) {
			//Predict Density
			interpolateDensity(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));

			if (particleList.rho(i) > rhoMax) {
				//Compute corrective pressure
				calculateCorrectivePressure(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

				//Increment pressure
				particleList.press(i) = particleList.press(i) + particleList.corrective_press(i);
			} else {
				//Add discrete particle forces to close particles
				calculateDiscreteParticleForces(i, particleList.nbr(i));
				/*particleList.predicted_vel(i) = particleList.vel(i) + dt * (particleList.grav(i) + particleList.visc(i) + 
					particleList.press(i) + particleList.fricoh(i) + particleList.disc(i));
				particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);*/
			/*}
		}*/

		for (int i = 0; i < particleList.size(); i++) {
			//Predict strain rate
			calculateStrainRate(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

			//Compute corrective dissipative stress
			calculateCorrectiveStress(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i),dt);

			//Increment stress
			particleList.stress(i) = particleList.stress(i) + particleList.corrective_stress(i);

			//Test yield and cohesion
			testYieldandCohesion(i);
		}	
		//stressTensor = calculateStressTensor(dt);
		//Ferr = calculateFerr();
	}
	//PUT BACK
	resolve_constraints(scene);
	for (int i = 0; i < particleList.size(); i++) {
		//Update velocity
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));//particleList.predicted_vel(i);
		//visc temp
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		//Update position (and check for collision - not in original algorithm but I want to do it here)
		//particleList.predicted_pos(i) = particleList.pos(i) + dt * particleList.vel(i);
		collision_detection(i, scene);
		particleList.pos(i) = particleList.predicted_pos(i);
	}
	updateSpheres();
}

void WaterSim::updateWater(const Scene* const scene, float dt) {
	//For all particles
	for (int i = 0; i < particleNum; i++) {
		//Apply forces Vi <- Vi + dt*fext(Xi) 
		//TODO: Are there other external forces??? Does mass come into play???
		particleList.vel(i) = particleList.vel(i) + dt * /*particleList.mass(i) **/ glm::vec3(0,-9.81,0);

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
	for (int i = 0; i < particleNum; i++) {
		//Update velocity v = (1/dt) * (Xi* - Xi)
		particleList.vel(i) = (1/dt) * (particleList.predicted_pos(i) - particleList.pos(i));

		//Apply vorticity confinement and XSPH viscosity
		//With just this nothing happens due to vorticity, particles hit floor and stay there
		//interpolateVorticity(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		//With just this column of bubbles happens
		applyXSPH(particleList.predicted_pos(i), nbrRadius, i, particleList.nbr(i));
		//for (int j = 0; j < particleList.nbr(i).size(); j++) {
			//interpolateViscosity(particleList.predicted_pos(j) - particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
		//}
			collision_detection(i, scene);

		//Update position Xi <- Xi*
		particleList.pos(i) = particleList.predicted_pos(i);
	}
	//particleToSphere();
	updateSpheres();
    //m_constraints_ext.clear();
}

void WaterSim::updateSpheres() {
	for (int i = 0; i < dots.size(); i++) {
		dots[i]->m_center = particleList.pos(i);
		dots[i]->mat_color = glm::vec3(abs(particleList.vel(i).x/glm::length(particleList.vel(i))), abs(particleList.vel(i).y/glm::length(particleList.vel(i))), abs(particleList.vel(i).z/glm::length(particleList.vel(i))));
		dots[i]->update();
	}
}

void WaterSim::resolve_constraints(const Scene* const scene) {
	for (unsigned int iter = 0; iter < solverIterations; ++iter) {
		//For all particles
		for (int i = 0; i < particleNum; i++) {
			//Calculate lambda
			//MOST TIME DRAIN HERE
			particleList.lambda(i) = calculateLambda(i); 
		}
		//For all particles
		for (int i = 0; i < particleNum; i++) {
			//Calculate deltaP
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
		for (int i = 0; i < particleNum; i++) {
			//Update position Xi* = Xi* + deltaP
			particleList.predicted_pos(i) = particleList.predicted_pos(i) + particleList.deltaP(i);
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
	return -R * temp;
	//-r should be the actual R)
}

glm::vec3 WaterSim::kernelVorticityDeriv(glm::vec3 R, float h) {
	return glm::vec3(0.0f);
}

glm::vec3 WaterSim::kernelGeneralDeriv(glm::vec3 R, float h) {
	float r = glm::length(R);
	float temp = (945.0f/((32.0f*PI*pow(h,9.0f))+0.000001f)) * pow((pow(h,2.0f)-pow(r,2.0f)),2.0f);
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
	float c = 0.01;
	for (int j = 0; j < nbrs.size(); j++) {
		if (particleList.type(nbrs[j]) != 0) {
			temp = R-particleList.predicted_pos(nbrs[j]);
			norm = glm::length(temp);//sqrt(pow(temp.x,2)+pow(temp.y,2)+pow(temp.z,2));
			//As += (particleList.vel(nbrs[j])-particleList.vel(i))*kernelViscosity(norm, h); 
			As += (particleList.mass(j)/(0.00001f + particleList.rho(j)))*(particleList.vel(nbrs[j])-particleList.vel(i))*kernelGeneral(norm,h);//SecondDeriv(temp, h); 
		}
	}
	particleList.vel(i) = particleList.vel(i) + c*As;
}

vector<int> WaterSim::surfaceSampling(vector<int> levelSet, float r, float t, float e) {
	vector<int> hola;
	return hola;
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
/*
//Modified from https://raw.githubusercontent.com/christopherbatty/Fluid3D/master/main.cpp
void WaterSim::export_particles(string path, int frame, const std::vector<glm::vec3>& particles, float radius) {
   //Write the output
   
   std::stringstream strout;
   strout << path << "particles_" << frame << ".txt";
   string filepath = strout.str();
   
   ofstream outfile(filepath.c_str());
   //write vertex count and particle radius
   outfile << particles.size() << " " << radius << std::endl;
   //write vertices
   for(unsigned int i = 0; i < particles.size(); ++i)
      outfile << particles[i][0] << " " << particles[i][1] << " " << particles[i][2] << std::endl;
   outfile.close();
}

//Modified from https://raw.githubusercontent.com/christopherbatty/Fluid3D/master/fluidsim.cpp
void WaterSim::compute_phi() {
   
   //grab from particles
   liquid_phi.assign(3*phi_dx);
   glm::vec3 tempv = glm::vec3(0.f);
   for(unsigned int p = 0; p < particleList.size(); ++p) {
      glm::vec3 cell_ind(particleList.predicted_pos(p) / phi_dx);
      for(int k = max(0.f, cell_ind[2]-1); k <= min(cell_ind[2]+1.f,nk-1.f); ++k) {
         for(int j = max(0.f,cell_ind[1] - 1); j <= min(cell_ind[1]+1.f,nj-1.f); ++j) {
            for(int i = max(0.f,cell_ind[0] - 1); i <= min(cell_ind[0]+1.f,ni-1.f); ++i) {
               glm::vec3 sample_pos((i+0.5f)*phi_dx, (j+0.5f)*phi_dx,(k+0.5f)*phi_dx);
               float test_val = glm::length(sample_pos - particleList.predicted_pos(p)) - particleRad;
               if(test_val < liquid_phi(i,j,k))
                  liquid_phi(i,j,k) = test_val;
            }
         }
      }
   }
   
   //extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
   Array3f phi_temp = liquid_phi;
   for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) {
         for(int i = 0; i < ni; ++i) {
            if(liquid_phi(i,j,k) < 0.5*phi_dx) {
               float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
                  + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
               if(solid_phi_val < 0)
                  phi_temp(i,j,k) = -0.5f*phi_dx;
            }
         }
      }
   }
   liquid_phi = phi_temp;
}*/



float WaterSim::interpolatePhi(glm::vec3 R, float h, int i, vector<int> nbrs) {
	float phi = 0;
	glm::vec3 r = glm::vec3(0.f);
	for (int j = 0; j < nbrs.size(); j++) {
		r = R - particleList.predicted_pos(nbrs[j]);
		phi += (particleList.mass(nbrs[j])/particleList.rho(nbrs[j])) * kernelAnisotropic(r,h);
	}
	return phi;
}

float WaterSim::kernelAnisotropic(glm::vec3 R, float h) {
	//Scaling factor
	float s = 1.f;
	//Dimension of simulation
	float d = 1.f;
	return (s/pow(h,d))*kernelAnisotropicSpline(R,h);
}

float WaterSim::kernelAnisotropicSpline(glm::vec3 R, float h) {
	float k = glm::length(R)/h;
	float a_D = 1.f;
	if (0 <= k && k < 1) {
		return a_D*(1 - (3.f/2.f)*pow(k,2.f) + (3.f/4.f)*pow(k,3.f));
	} else if (1 <= k && k < 2) {
		return a_D*((1.f/4.f)*pow(2-k,3.f));
	} else {
		return 0;
	}
}

void WaterSim::marchingCube(const VBO& vbos)
{
	phis.resize(133);
	for (int i = 0; i < 133; i++) {
		phis[i].resize(133);
		for (int j = 0; j < 133; j++) {
			phis[i][j].resize(133);
		}
	}
	for (int i = 0; i < particleList.size(); i++) {
		interpolatePhi(particleList.predicted_pos(i),nbrRadius,i,particleList.nbr(i));
	}

	for(auto kv = table.begin(); kv != table.end(); kv++) {
		keys.push_back(kv->first);
		vals.push_back(kv->second);  
	} 
	phiCell.resize(vals.size());

	// vbo is used for displaying triangle mesh generated by marching cube with OpenGL.
	std::vector<float> position;// position for all triangle vertices
	std::vector<float> normal;// normal for all triangle vertices
	std::vector<float> color;// color for all triangle vertices
	std::vector<unsigned int> indices;// triangle indices

	vector<glm::vec3> m_positions;
	vector<glm::vec3> m_normals;
	vector<glm::vec3> m_colors;
	vector<unsigned short> m_indices;

	glm::vec3 mat_color = glm::vec3(0.0,0.0,1.0);// color for the triangle mesh generated
	unsigned int triangleCount(0);

	float iso_value(0.0f);// the value for surface, by default is 0.0f
	float key = 0.0;
	float avg = 0.0;

	// Average the phi values in every cell
	for (int i = 0; i < 133; i++) {
		for (int j = 0; j < 133; j++) {
			for (int k = 0; k < 133; k++) {
				avg = 0.0;
				vector<int> neighbors = table[grid[i][j][k]];
				for (int h = 0; h < neighbors.size(); h++) {
					avg += particleList.phi(neighbors[h]);
				}
				avg /= (vals[i].size() + 0.000001);
				phis[i][j][k] = avg;
			}
		}
	}
	// Loop over all cells in the grid.
	for (int i = 0; i < 133; i++) {
		for (int j = 0; j < 133; j++) {
			for (int k = 0; k < 133; k++) {
	// for a current cell ijk:
			
				//1. sample values of 8 corners for current cell.
				float cornerValue[8];
				unsigned int valueMask = 0;//8 bit for indicating polarity for 8 corner.
				// Loop over 8 corner
				for (int current_corner = 0; current_corner < 8; current_corner++)
				{
					float corner_value = 0.0f;// get sample from the grid.
					if (i != 0 && j != 0 && k != 0 && i != 132 && j != 132 && k != 132) {
						if (current_corner == 0) {
							corner_value = (phis[i][j][k] + phis[i][j-1][k] + phis[i-1][j-1][k] + phis[i-1][j][k] + phis[i][j][k+1] 
							+ phis[i-1][j][k+1] + phis[i][j-1][k+1] + phis[i-1][j-1][k+1])/8.0f;
						} else if (current_corner == 1) {
							corner_value = (phis[i][j][k] + phis[i][j-1][k] + phis[i+1][j-1][k] + phis[i+1][j][k] + phis[i][j][k+1] 
							+ phis[i+1][j][k+1] + phis[i][j-1][k+1] + phis[i+1][j-1][k+1])/8.0f;
						} else if (current_corner == 2) {
							corner_value = (phis[i][j][k] + phis[i][j-1][k] + phis[i+1][j-1][k] + phis[i+1][j][k] + phis[i][j][k-1] 
							+ phis[i+1][j][k-1] + phis[i][j-1][k-1] + phis[i+1][j-1][k-1])/8.0f;
						} else if (current_corner == 3) {
							corner_value = (phis[i][j][k] + phis[i][j-1][k] + phis[i-1][j-1][k] + phis[i][j-1][k-1] + phis[i-1][j-1][k-1] 
							+ phis[i-1][j][k] + phis[i][j][k-1] + phis[i-1][j][k-1])/8.0f;
						} else if (current_corner == 4) {
							corner_value = (phis[i][j][k] + phis[i][j+1][k] + phis[i-1][j+1][k] + phis[i-1][j][k] + phis[i][j][k+1] 
							+ phis[i-1][j][k+1] + phis[i][j+1][k+1] + phis[i-1][j+1][k+1])/8.0f;
						} else if (current_corner == 5) {
							corner_value = (phis[i][j][k] + phis[i][j+1][k] + phis[i+1][j+1][k] + phis[i+1][j][k] + phis[i][j][k+1] 
							+ phis[i+1][j][k+1] + phis[i][j+1][k+1] + phis[i+1][j+1][k+1])/8.0f;
						} else if (current_corner == 6) {
							corner_value = (phis[i][j][k] + phis[i][j+1][k] + phis[i+1][j+1][k] + phis[i+1][j][k] + phis[i][j][k-1] 
							+ phis[i+1][j][k-1] + phis[i][j+1][k-1] + phis[i+1][j+1][k-1])/8.0f;
						} else {
							corner_value = (phis[i][j][k] + phis[i][j+1][k] + phis[i][j+1][k-1] + phis[i-1][j+1][k] + phis[i][j][k-1] 
							+ phis[i-1][j][k] + phis[i-1][j+1][k-1] + phis[i-1][j][k-1])/8.0f;
						}
					}
					iso_value = phis[i][j][k];
					cornerValue[current_corner] = corner_value - iso_value;

					if(cornerValue[current_corner] <= 0)
						valueMask |= (1 << current_corner);// bit operation to mark corresponding corner.
				}

				//2. get vertex position and normal for intersection point.
				int edgeFlag = aiCubeEdgeFlags[valueMask];// flag for indicating which edges is intersecting
				// not intersecting at all.
				if(edgeFlag == 0)
					continue;
				// 12 correspond to 12 edges. e.g. a intersection vertex on edge 6 would be stored in edgeVertex[6] and edgeNormal[6]
				glm::vec3 edgeVertex[12];
				glm::vec3 edgeNormal[12];
				// Loop over 12 edges.
				for (int edge_id = 0; edge_id < 12; edge_id++)
				{
					// if edge is intersecting with the surface.
					//1 for intersection
					if ((aiCubeEdgeFlags[edge_id] & (1 << edge_id)) != 0) {//HELP WITH THIS!!!!!!!!!!!!!
					// use bit-wise operation to find out. intersection information is stored in edgeFlag.
					{
						// get the two end point ids for current edge.
						unsigned int v1 = a2iEdgeConnection[edge_id][0];
						unsigned int v2 = a2iEdgeConnection[edge_id][1];
						// fraction is the intersection point between v1 and v2.
						float fraction = cornerValue[v1] / ((cornerValue[v1] - cornerValue[v2])+0.000001);
						// edgeDirection is the unit vector of the direction of current edge in world space.
						glm::vec3 edgeDirection = glm::vec3(a2fEdgeDirection[edge_id][0],a2fEdgeDirection[edge_id][1],a2fEdgeDirection[edge_id][2]);
						edgeDirection /= glm::length(edgeDirection);
						// get edge direction from a2fEdgeDirection[edge_id]

						glm::vec3 pos;
						unsigned int i1, j1, k1;
						//TO DO: INDEXING
						// i, j, k is the index for the corner with id 0;
						// i1, j1, k1 is the index for the corner with id v1;
						
						// computer i1, j1, k1;
						if (v1 == 0) {
							i1 = i;
							j1 = j;
							k1 = k;
						} else if (v1 == 1) {
							i1 = i+1;
							j1 = j;
							k1 = k;
						} else if (v1 == 2) {
							i1 = i+1;
							j1 = j;
							k1 = k-1;
						} else if (v1 == 3) {
							i1 = i;
							j1 = j;
							k1 = k-1;
						} else if (v1 == 4) {
							i1 = i;
							j1 = j+1;
							k1 = k;
						} else if (v1 == 5) {
							i1 = i+1;
							j1 = j+1;
							k1 = k;
						} else if (v1 == 6) {
							i1 = i+1;
							j1 = j+1;
							k1 = k-1;
						} else {
							i1 = i;
							j1 = j+1;
							k1 = k-1;
						}

                        // compute position of corner i1, j1, k1 in world space;

                        float cellsize = 0.3;//CHECK THIS

						pos += fraction * cellsize * edgeDirection;

						edgeVertex[edge_id] = pos;
                        
						//TO DO:
                        // compute normal for current vertex.
                        // works only for level set field, since gradient of level set is the normal.

                        glm::vec3 norm;
                        // compute normal using the gradient.
						glm::vec3 a1, b1, c1, a2, b2, c2, a3, b3, c3;
						if (v1 == 0) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1+1,k1);
							c1 = glm::vec3(i1+1,j1,k1);

							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1,k1-1);
							c2 = glm::vec3(i1,j1+1,k1);

							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1+1,j1,k1);
							c3 = glm::vec3(i1,j1,k1-1);
						} else if (v1 == 1) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1+1,k1);
							c1 = glm::vec3(i1,j1,k1-1);

							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1-1,j1,k1);
							c2 = glm::vec3(i1,j1+1,k1);

							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1,j1,k1-1);
							c3 = glm::vec3(i1-1,j1,k1);
						} else if (v1 == 2) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1+1,k1);
							c1 = glm::vec3(i1,j1,k1+1);

							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1,k1+1);
							c2 = glm::vec3(i1,j1+1,k1);

							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1-1,j1,k1);
							c3 = glm::vec3(i1,j1,k1+1);
						} else if (v1 == 3) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1+1,j1,k1);
							c1 = glm::vec3(i1,j1+1,k1);

							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1+1,k1);
							c2 = glm::vec3(i1,j1,k1+1);

							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1,j1,k1+1);
							c3 = glm::vec3(i1+1,j1,k1);
						} else if (v1 == 4) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1,k1-1);
							c1 = glm::vec3(i1+1,j1,k1);
											
							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1-1,k1);
							c2 = glm::vec3(i1,j1,k1-1);

							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1+1,j1,k1);
							c3 = glm::vec3(i1,j1-1,k1);
						} else if (v1 == 5) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1-1,j1,k1);
							c1 = glm::vec3(i1,j1,k1-1);
											
							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1,k1-1);
							c2 = glm::vec3(i1,j1-1,k1);
											
							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1,j1-1,k1);
							c3 = glm::vec3(i1-1,j1,k1);
						} else if (v1 == 6) {
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1,k1+1);
							c1 = glm::vec3(i1-1,j1,k1);

							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1,j1-1,k1);
							c2 = glm::vec3(i1,j1,k1+1);
											
							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1-1,j1,k1);
							c3 = glm::vec3(i1,j1-1,k1);
						} else {			
							a1 = glm::vec3(i1,j1,k1);
							b1 = glm::vec3(i1,j1-1,k1);
							c1 = glm::vec3(i1+1,j1,k1);
											
							a2 = glm::vec3(i1,j1,k1);
							b2 = glm::vec3(i1+1,j1,k1);
							c2 = glm::vec3(i1,j1,k1+1);
											
							a3 = glm::vec3(i1,j1,k1);
							b3 = glm::vec3(i1,j1,k1+1);
							c3 = glm::vec3(i1,j1-1,k1);
						}

						glm::vec3 d1,e1,d2,e2,d3,e3;
						d1 = b1-a1;
						e1 = c1-a1;
						glm::vec3 n1 = glm::cross(d1,e1);
						n1 /= (n1.x+n1.y+n1.z)+0.000001;

						d2 = b2-a2;
						e2 = c2-a2;
						glm::vec3 n2 = glm::cross(d2,e2);
						n2 /= (n2.x+n2.y+n2.z)+0.000001;

						d3 = b3-a3;
						e3 = c3-a3;
						glm::vec3 n3 = glm::cross(d3,e3);
						n3 /= (n3.x+n3.y+n3.z)+0.000001;
                        // take samples from the gradient.
						
						norm = n1+n2+n3;

                        norm = glm::normalize(norm);
                        edgeNormal[edge_id] = norm;
					}
				}

				//3. generate triangles from the vertices and normal.
				// triangle information is stored in a2iTriangleConnectionTable[valueMask]
				int * triangleTablePtr = a2iTriangleConnectionTable[valueMask];
				for(unsigned int num = 0; num < 5; ++num)
				{// up to 5 triangles
					if(*(triangleTablePtr + 3 * num) < 0)
						break;

					for(unsigned int idx = 0; idx < 3; ++idx)
					{
						// vertex id is used for extracting position and normal from edgeVertex and edgeNormal array.
						int vertex_idx = *(triangleTablePtr + 3 * num + idx);

						glm::vec3& p = edgeVertex[vertex_idx];

						position.push_back(p.x);
						position.push_back(p.y);
						position.push_back(p.z);

						m_positions.push_back(p);

                        glm::vec3& n = edgeNormal[vertex_idx];

                        normal.push_back(n.x);
                        normal.push_back(n.y);
                        normal.push_back(n.z);

						m_normals.push_back(n);

						color.push_back(mat_color.r);
						color.push_back(mat_color.g);
						color.push_back(mat_color.b);

						m_colors.push_back(mat_color);

						indices.push_back(3 * triangleCount + idx);

						m_indices.push_back(3 * triangleCount + idx);
					}

					triangleCount++;
				}
	}
			}
		}
	}
	// post work: clean up and upload mesh to somewhere.
	////
	triangle* t = new triangle();
	t->draw(vbos,m_positions,m_colors,m_normals,m_indices);
}