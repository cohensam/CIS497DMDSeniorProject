#include "point.h"
#include "math_headers.h"
#include "openGL_headers.h"
#include "glm.hpp"
#include <vector>

class triangle {
public:
	triangle();
	~triangle();
	glm::vec3 p0;
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 color;
	glm::vec3 normal();
	glm::vec3 a;
	glm::vec3 b;
	glm::vec3 c;
	glm::vec3 norm;
	void draw(const VBO& vbos, std::vector<glm::vec3> m_positions, std::vector<glm::vec3> m_colors, std::vector<glm::vec3> m_normals, std::vector<unsigned short> m_indices) const;
};