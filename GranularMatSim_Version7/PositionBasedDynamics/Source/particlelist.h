// Cloth Simulation using Position Based Dynamics
// Copyright 2013 Xing Du

#ifndef PARTICLE_LIST_H_INCLUDED
#define PARTICLE_LIST_H_INCLUDED

#include "math_headers.h"
#include <vector>
#include <cassert>

// Structure of Array
struct ParticleList
{
    ParticleList() {};
    virtual ~ParticleList()
    {
        m_pos.clear();
        m_predicted_pos.clear();
        m_vel.clear();
		m_predicted_vel.clear();
        m_lock_pos.clear();
        m_mass.clear();
        m_rho.clear();
		m_visc.clear();
		m_press.clear();
		m_vort.clear();
		m_nbr.clear();
		m_deltaP.clear();
		m_lambda.clear();
		m_phi.clear();
		m_size = 0;

		m_grav.clear();
		m_stress.clear();
		m_fricoh.clear();
		m_corrective_press.clear();
		m_corrective_stress.clear();
		m_strain.clear();
		m_disc.clear();
		m_press_force.clear();

		m_type.clear();
		m_bound.clear();

		m_coh.clear();

		m_parent_id.clear();
    }

    void resize(unsigned int size)
    {
        m_pos.clear();
        m_predicted_pos.clear();
        m_vel.clear();
		m_predicted_vel.clear();
        m_lock_pos.clear();
        m_mass.clear();
		m_phi.clear();

		m_nbr.clear();
		m_deltaP.clear();
		m_lambda.clear();
		m_rho.clear();
		m_visc.clear();
		m_press.clear();
		m_vort.clear();
		m_vort.resize(size);
		m_rho.resize(size);
		m_visc.resize(size);
		m_press.resize(size);
		m_nbr.resize(size);
		m_deltaP.resize(size);
		m_lambda.resize(size);
		m_phi.resize(size);

        m_pos.resize(size);
        m_predicted_pos.resize(size);
        m_vel.resize(size);
		m_predicted_vel.resize(size);
        m_lock_pos.resize(size, false);
        m_mass.resize(size, 1.0f);

        m_size = size;

		m_grav.clear();
		m_grav.resize(size);
		m_stress.clear();
		m_stress.resize(size);
		m_fricoh.clear();
		m_fricoh.resize(size);
		m_corrective_press.clear();
		m_corrective_press.resize(size);
		m_corrective_stress.clear();
		m_corrective_stress.resize(size);
		m_strain.clear();
		m_strain.resize(size);
		m_disc.clear();
		m_disc.resize(size);
		m_press_force.clear();
		m_press_force.resize(size);

		m_type.clear();
		m_type.resize(size);
		m_bound.clear();
		m_bound.resize(size);

		m_coh.clear();
		m_coh.resize(size);

		m_parent_id.clear();
		m_parent_id.resize(size);
    }
    // these two functions are used for attach a vertex to a certain point.
    void lock_pos(unsigned int n)
    {
        assert(n < m_size);
        m_lock_pos[n] = true;
    }
    void unlock_pos(unsigned int n)
    {
        assert(n < m_size);
        m_lock_pos[n] = false;
    }
    void unlock_pos_all()
    {
        for(std::vector<bool>::iterator i = m_lock_pos.begin(); i != m_lock_pos.end(); ++i)
        {
            *i = false;
        }
    }
    // inverse mass accessor.
    float mass(unsigned int n) const
    {
        assert(n < m_size);
        if(m_lock_pos[n])
            return 0.0f;
        else
            return m_mass[n];
    }
    void set_mass(unsigned int n, float mass)
    {
        assert(n < m_size);
        m_mass[n] = mass;
    }
    // position accessor.
    const glm::vec3& pos(unsigned int n) const
    {
        assert(n < m_size);
        return m_pos[n];
    }
    glm::vec3& pos(unsigned int n)
    {
        assert(n < m_size);
        return m_pos[n];
    }
    // predicted position accessor.
    const glm::vec3& predicted_pos(unsigned int n) const
    {
        assert(n < m_size);
        return m_predicted_pos[n];
    }
    glm::vec3& predicted_pos(unsigned int n)
    {
        assert(n < m_size);
        return m_predicted_pos[n];
    }
    // velocity accessor.
    const glm::vec3& vel(unsigned int n) const
    {
        assert(n < m_size);
        return m_vel[n];
    }
    glm::vec3& vel(unsigned int n)
    {
        assert(n < m_size);
        return m_vel[n];
    }
	// predicted velocity accessor.
    const glm::vec3& predicted_vel(unsigned int n) const
    {
        assert(n < m_size);
        return m_predicted_vel[n];
    }
    glm::vec3& predicted_vel(unsigned int n)
    {
        assert(n < m_size);
        return m_predicted_vel[n];
    }
    // size accessor
    unsigned int size() const
    {
        return m_size;
    }

	// grid accessor.
    glm::vec3 grid(unsigned int n) const
    {
        assert(n < m_size);
        return m_grid[n];
    }

	// density accessor.
    const float& rho(unsigned int n) const
    {
        assert(n < m_size);
        return m_rho[n];
    }
    float& rho(unsigned int n)
    {
        assert(n < m_size);
        return m_rho[n];
    }

	// viscosity accessor
    const glm::vec3& visc(unsigned int n) const
    {
        assert(n < m_size);
        return m_visc[n];
    }
    glm::vec3& visc(unsigned int n)
    {
        assert(n < m_size);
        return m_visc[n];
    }

	// pressure accessor
	const glm::vec3& press(unsigned int n) const
    {
        assert(n < m_size);
        return m_press[n];
    }
    glm::vec3& press(unsigned int n)
    {
        assert(n < m_size);
        return m_press[n];
    }

	// vorticity accessor
	const glm::vec3& vort(unsigned int n) const
    {
        assert(n < m_size);
        return m_vort[n];
    }
    glm::vec3& vort(unsigned int n)
    {
        assert(n < m_size);
        return m_vort[n];
    }

	// lambda accessor.
	const float& lambda(unsigned int n) const
    {
        assert(n < m_size);
        return m_lambda[n];
    }
    float& lambda(unsigned int n)
    {
        assert(n < m_size);
        return m_lambda[n];
    }

	// deltaP accessor.
    const glm::vec3& deltaP(unsigned int n) const
    {
        assert(n < m_size);
        return m_deltaP[n];
    }
    glm::vec3& deltaP(unsigned int n)
    {
        assert(n < m_size);
        return m_deltaP[n];
    }

	//neighbor accessor
	const std::vector<int>& nbr(unsigned int n) const
    {
        assert(n < m_size);
        return m_nbr[n];
    }
    std::vector<int>& nbr(unsigned int n)
    {
        assert(n < m_size);
        return m_nbr[n];
    }

	//phi accessor
	const float& phi(unsigned int n) const
    {
        assert(n < m_size);
        return m_phi[n];
    }
    float& phi(unsigned int n)
    {
        assert(n < m_size);
        return m_phi[n];
    }

	//gravity accessor
	 const glm::vec3& grav(unsigned int n) const
    {
        assert(n < m_size);
        return m_grav[n];
    }
    glm::vec3& grav(unsigned int n)
    {
        assert(n < m_size);
        return m_grav[n];
    }

	//stress accessor
	 const glm::vec3& stress(unsigned int n) const
    {
        assert(n < m_size);
        return m_stress[n];
    }
    glm::vec3& stress(unsigned int n)
    {
        assert(n < m_size);
        return m_stress[n];
    }

	//friction and cohesion accessor
	 const glm::vec3& fricoh(unsigned int n) const
    {
        assert(n < m_size);
        return m_fricoh[n];
    }
    glm::vec3& fricoh(unsigned int n)
    {
        assert(n < m_size);
        return m_fricoh[n];
    }

	//corrective pressure accessor
	 const glm::vec3& corrective_press(unsigned int n) const
    {
        assert(n < m_size);
        return m_corrective_press[n];
    }
    glm::vec3& corrective_press(unsigned int n)
    {
        assert(n < m_size);
        return m_corrective_press[n];
    }

	//corrective stress accessor
	 const glm::vec3& corrective_stress(unsigned int n) const
    {
        assert(n < m_size);
        return m_corrective_stress[n];
    }
    glm::vec3& corrective_stress(unsigned int n)
    {
        assert(n < m_size);
        return m_corrective_stress[n];
    }

	//strain accessor
	 const glm::vec3& strain(unsigned int n) const
    {
        assert(n < m_size);
        return m_strain[n];
    }
    glm::vec3& strain(unsigned int n)
    {
        assert(n < m_size);
        return m_strain[n];
    }

	//discrete particle forces accessor
	 const glm::vec3& disc(unsigned int n) const
    {
        assert(n < m_size);
        return m_disc[n];
    }
    glm::vec3& disc(unsigned int n)
    {
        assert(n < m_size);
        return m_disc[n];
    }

	// particle type accessor.
    const float& type(unsigned int n) const
    {
        assert(n < m_size);
        return m_type[n];
    }
    float& type(unsigned int n)
    {
        assert(n < m_size);
        return m_type[n];
    }

	// pressure force accessor
	const glm::vec3& press_force(unsigned int n) const
    {
        assert(n < m_size);
        return m_press_force[n];
    }
    glm::vec3& press_force(unsigned int n)
    {
        assert(n < m_size);
        return m_press_force[n];
    }

	// boundary particle force accessor
	const glm::vec3& bound(unsigned int n) const
    {
        assert(n < m_size);
        return m_bound[n];
    }
    glm::vec3& bound(unsigned int n)
    {
        assert(n < m_size);
        return m_bound[n];
    }

	// cohesion force accessor
	const glm::vec3& coh(unsigned int n) const
    {
        assert(n < m_size);
        return m_coh[n];
    }
    glm::vec3& coh(unsigned int n)
    {
        assert(n < m_size);
        return m_coh[n];
    }

	// parent particle accessor.
    const float& parent_id(unsigned int n) const
    {
        assert(n < m_size);
        return m_parent_id[n];
    }
    float& parent_id(unsigned int n)
    {
        assert(n < m_size);
        return m_parent_id[n];
    }

protected:
    // number of vertices.
    unsigned int m_size;
    // position of all vertices.
    std::vector<glm::vec3> m_pos;
	// grid indices of all vertices.
    std::vector<glm::vec3> m_grid;
    // predicted position.
    std::vector<glm::vec3> m_predicted_pos;
    // velocity of a vertex
    std::vector<glm::vec3> m_vel;
    // predicted velocity
    std::vector<glm::vec3> m_predicted_vel;
    // used for fixing vertex to a certain point
    std::vector<bool> m_lock_pos;
    // weight for resolving the constraints.
    std::vector<float> m_mass;
	// Phi values for level set calculations
    std::vector<float> m_phi;

	// density of all vertices.
    std::vector<float> m_rho;
	// viscosity of all vertices.
    std::vector<glm::vec3> m_visc;
	// pressure of all vertices.
    std::vector<glm::vec3> m_press;
	// vorticity of all vertices.
    std::vector<glm::vec3> m_vort;

	// neighbors
    std::vector<std::vector<int>> m_nbr;

	// constraint projection 
	std::vector<float> m_lambda; 
	std::vector<glm::vec3> m_deltaP; 

	// gravity
    std::vector<glm::vec3> m_grav;
	// dissipative stress
    std::vector<glm::vec3> m_stress;
	//friction and cohesion
	std::vector<glm::vec3> m_fricoh;
	//corrective pressure
	std::vector<glm::vec3> m_corrective_press;
	//corrective dissipative stress
	std::vector<glm::vec3> m_corrective_stress;
	//strain rate
	std::vector<glm::vec3> m_strain;
	//discrete particle forces
	std::vector<glm::vec3> m_disc;
	//pressure force
	std::vector<glm::vec3> m_press_force;

	/*particle type: 
	-1->surrounding particles to change low-res shape
	0->boundary 
	1->granular low res 
	2->partially wet (not implemented yet) 
	3->wet (not implemented yet) 
	4->water (implemented separately need to combine)
	5->granular high res*/
	std::vector<float> m_type; 
	//boundary particle force
	std::vector<glm::vec3> m_bound;
	//cohesion accessor
	std::vector<glm::vec3> m_coh;
	//id of parent particle for low res particles' surrounding particles to make them oblong
	std::vector<float> m_parent_id; 
};
#endif