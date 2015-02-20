// Cloth Simulation using Position Based Dynamics
// Courtesy of Aline Normoyle
// Copyright 2013 Xing Du

#ifndef SPHERE_H_INCLUDED
#define SPHERE_H_INCLUDED

#include "openGL_headers.h"
#include "math_headers.h"
#include "glm.hpp"
#include <tinyxml.h>
#include <vector>

class Sphere
    {
    public:
        Sphere() :  m_center(glm::vec3(0.0f)), m_radius(0.0f) {init_visualization();};
        Sphere(const glm::vec3 pos, float radius, glm::vec3 color) : m_center(pos), m_radius(radius) {mat_color = color; init_visualization();};
        Sphere(const Sphere& other) : 
			m_positions(other.m_positions), m_colors(other.m_colors), m_normals(other.m_normals),
            m_indices(other.m_indices), m_center(other.m_center), m_radius(other.m_radius) {};
        virtual ~Sphere()
        {
            m_positions.clear();
            m_colors.clear();
            m_normals.clear();
            m_indices.clear();
        }

	protected:
        virtual void init_visualization();
    protected:
        std::vector<glm::vec3> m_positions, m_colors, m_normals;
        std::vector<unsigned short> m_indices;
	public:
		glm::vec3 mat_color;
        virtual void draw(const VBO& vbos) const;
        virtual bool line_intersection(const glm::vec3& p1, const glm::vec3& p2, float threshold, glm::vec3& intersect, glm::vec3& normal) const;
		virtual void update();
	//protected:
	public:
        glm::vec3 m_center;
        float m_radius;
    };

#endif