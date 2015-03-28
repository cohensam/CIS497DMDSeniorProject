#include "triangle.h"

triangle::triangle() {
	p0 = glm::vec3(0,0,0);
	p1 = glm::vec3(0,0,0);
	p2 = glm::vec3(0,0,0);
	a = glm::vec3(0,0,0);
	b = glm::vec3(0,0,0);
	c = glm::vec3(0,0,0);
	color = glm::vec3(0,0,0);
	norm = glm::vec3(0,0,0);
}

triangle::~triangle(){

}

glm::vec3 triangle::normal() {
	//a = p1 - p0;
	//b = p2 - p0;
	glm::vec3 V = p1 - p0;
	glm::vec3 W = p2 - p0;
	glm::vec3 N;
	N.x = (glm::dot(V.y,W.z)-glm::dot(V.z,W.y));
	N.y = (glm::dot(V.z,W.x)-glm::dot(V.x,W.z));
	N.z = (glm::dot(V.x,W.y)-glm::dot(V.y,W.x));
	float denom = N.x + N.y + N.z;
	N.x /= denom;
	N.y /= denom;
	N.z /= denom;
	return N;
	/*if (glm::length(glm::cross(a,b)) != 0) {
		c = glm::cross(a,b) / glm::length(glm::cross(a,b));
	} else {
		c = glm::cross(a,b);
	}
	return c;*/
}

void triangle::draw(const VBO& vbos, std::vector<glm::vec3> m_positions, std::vector<glm::vec3> m_colors, std::vector<glm::vec3> m_normals, std::vector<unsigned short> m_indices) const {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    // position
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_STREAM_DRAW);

    // color
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_STREAM_DRAW);

    // normal
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_SHORT, 0);//GL_UNSIGNED_INT

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}