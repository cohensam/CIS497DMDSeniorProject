#include "sphere.h"
#include <cassert>

//----------Sphere Class----------//
void Sphere::init_visualization()
{
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    //mat_color = glm::vec3(0.0f, 0.0f, 1.0f);//0.6f);
    unsigned int slice = 24, stack = 10;

    glm::vec3 tnormal(0.0f, 1.0f, 0.0f), tpos;
	tpos = m_center + m_radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

	float theta_z, theta_y, sin_z;
    float delta_y = 360.0f / slice, delta_z = 180.0f / stack;
	//loop over the sphere
	for(theta_z = delta_z; theta_z < 179.99f; theta_z += delta_z)
	{
		for(theta_y = 0.0f; theta_y < 359.99f; theta_y += delta_y)
		{
			sin_z = sin(glm::radians(theta_z));
			
            tnormal.x = sin_z * cos(glm::radians(theta_y));
			tnormal.y = cos(glm::radians(theta_z));
			tnormal.z = -sin_z * sin(glm::radians(theta_y));

			tpos = m_center + m_radius * tnormal;

            m_positions.push_back(tpos);
            m_normals.push_back(tnormal);
            m_colors.push_back(mat_color);
		}
	}
	tnormal = glm::vec3(0.0f, -1.0f, 0.0f);
    tpos = m_center + m_radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

	//indices
	unsigned int j = 0, k = 0;
	for(j = 0; j < slice - 1; ++j)
	{
		m_indices.push_back(0);
		m_indices.push_back(j + 1);
		m_indices.push_back(j + 2);
	}
	m_indices.push_back(0);
	m_indices.push_back(slice);
	m_indices.push_back(1);

	for(j = 0; j < stack - 2; ++j)
	{
		for(k = 1 + slice * j; k < slice * (j + 1); ++k)
		{
			m_indices.push_back(k);
			m_indices.push_back(k + slice);
			m_indices.push_back(k + slice + 1);

			m_indices.push_back(k);
			m_indices.push_back(k + slice + 1);
			m_indices.push_back(k + 1);
		}
		m_indices.push_back(k);
		m_indices.push_back(k + slice);
		m_indices.push_back(k + 1);

		m_indices.push_back(k);
		m_indices.push_back(k + 1);
		m_indices.push_back(k + 1 - slice);
	}

    unsigned int bottom_id = (stack - 1) * slice + 1;
    unsigned int offset = bottom_id - slice;
	for(j = 0; j < slice - 1; ++j)
	{
		m_indices.push_back(j + offset);
		m_indices.push_back(bottom_id);
		m_indices.push_back(j + offset + 1);
	}
	m_indices.push_back(bottom_id - 1);
	m_indices.push_back(bottom_id);
	m_indices.push_back(offset);

	if(m_indices.size() != 6 * (stack - 1) * slice)
		printf("indices number not correct!\n");
}

void Sphere::draw(const VBO& vbos) const
{
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
	
	/*
	//Makes points instead
	glHint(GL_POINT_SMOOTH, GL_NICEST);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(100);
	glColor4f(0,0,1,1);
	glBegin(GL_POINTS);
		glVertex3f(m_center.x,m_center.y,m_center.z);
	glEnd();
	*/
}

bool Sphere::line_intersection(const glm::vec3& p1, const glm::vec3& p2, float threshold, glm::vec3& intersect, glm::vec3& normal) const
{// TODO: implement line-sphere intersection. you can refer to line-plane intersection.
    
	//JUST USE NORMAL SPHERE INTERSECTION
	
	float v1, v2;// v1 v2 are distance to sphere for p1 and p2.
    v1 = abs(glm::length(p1 - m_center) - (m_radius+threshold));//glm::dot(p1, normal) - m_radius;  // m_value;
    v2 = abs(glm::length(p2 - m_center) - (m_radius+threshold));//glm::dot(p2, normal) - m_radius;  // m_value;
    if(v2 < threshold)
    {
        if(v1 >= threshold)
        {// continuous collision handling.
			//Try line-sphere intersection use math from online to calculate the intersection
			float t = 0;
			float a = pow((p1.x-m_center.x),2.f) + pow((p1.y-m_center.y),2.f) + pow((p1.z-m_center.z),2.f) - pow((m_radius+threshold),2.f);//A = (x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2 - R^2
			float c = pow((p1.x-p2.x),2.f) + pow((p1.y-p2.y),2) + pow((p1.z-p2.z),2);					//C = (x0-x1)^2 + (y0-y1)^2 + (z0-z1)^2
			float b = pow((p2.x-m_center.x),2) + pow((p2.y-m_center.y),2) + pow((p2.z-m_center.z),2) - a - c - pow((m_radius+threshold),2.f);//B = (x1-xc)^2 + (y1-yc)^2 + (z1-zc)^2 - A - C - R^2
			if (a != 0 && b * b - 4 * a * c >= 0) {
				t = (-b - sqrt(b * b - 4 * a * c))/(2 * a);
			}
			intersect = glm::vec3(p1.x*(1.0-t)+t*p2.x,p1.y*(1.0-t)+t*p2.y,p1.z*(1.0-t)+t*p2.z);
			//x(t) = x0*(1-t) + t*x1
			//y(t) = y0*(1-t) + t*y1
			//z(t) = z0*(1-t) + t*z1
            normal = glm::normalize(intersect - m_center);//glm::vec3(0.0f);
        }
        else
        {// static collision handling.
			//Need to find closest point on sphere to p2, and this is the intersection point
            glm::vec3 dir = glm::normalize(p2 - m_center);//p1-p2);//
			//intersect = p2 - (v2 - threshold) * normal;
            intersect = p2 + (dir*v2);//m_center + (dir * v2);//dir * (m_radius + threshold);//m_center + dir * v2; //CHECK THIS MATH ONE OF THESE IS WRONG
			normal = glm::normalize(intersect - m_center);//glm::vec3(0.0f);
			//Remembr to NORMALIZE THE NORMAL
			
        }
        return true;
    }
    else
        return false;
}
void Sphere::update() {

    glm::vec3 tnormal(0.0f, 1.0f, 0.0f), tpos;

	//loop over the sphere
	for(int i = 0; i <  m_positions.size(); i++)
	{
		tpos = m_center + m_radius * m_normals[i];

		m_positions[i] = tpos;
        m_colors[i] = mat_color;
	}

	tnormal = glm::vec3(0.0f, -1.0f, 0.0f);
    tpos = m_center + m_radius * tnormal;

}