#ifndef CONSTRAINT_H_INCLUDED
#define CONSTRAINT_H_INCLUDED

#include "math_headers.h"
#include "particlelist.h"
#include <vector>

class Constraint
{
public:
    Constraint();
    Constraint(ParticleList *verts, float stiff);
    Constraint(const Constraint& other);
    virtual ~Constraint();
    
    void set_stiffness(float k)
    {
        m_stiffness = k;
    }
    virtual bool project_constraint();
protected:
    // a pointer to all the vertices.
    ParticleList *m_vertices;
    // stiffness decide how much the vertex would move towards the constrained position.
    float m_stiffness;
};

class CollisionConstraint : public Constraint
{
public:
    CollisionConstraint();
    CollisionConstraint(ParticleList *verts, unsigned int p0, const glm::vec3& q, const glm::vec3& n);
    CollisionConstraint(const CollisionConstraint& other);
    virtual ~CollisionConstraint();

    virtual bool project_constraint();
    
    const glm::vec3& normal() const
    {
        return m_normal;
    }
    unsigned int index() const
    {
        return m_p0;
    }
protected:
    // cardinality for bend constraint is 1.
    unsigned int m_p0;
    glm::vec3 m_ref_point, m_normal;
};
#endif