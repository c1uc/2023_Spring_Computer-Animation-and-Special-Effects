#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    int particleNum = jelly.getParticleNum();
    for (int i = 0; i < particleNum; i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f position = particle.getPosition();
        Eigen::Vector3f velocity = particle.getVelocity();
        Eigen::Vector3f force = particle.getForce();

        Eigen::Vector3f velocityN = this->normal.dot(velocity) * this->normal;
        Eigen::Vector3f velocityT = velocity - velocityN;

        Eigen::Vector3f dist = position - this->hole_position;
        dist[1] = 0.0f;
        if (dist.norm() <= this->hole_radius) continue;

        if (this->normal.dot(position - this->position) < eEPSILON && this->normal.dot(velocity) < 0) {
            Eigen::Vector3f velocity_ = velocityT - coefResist * velocityN;
            particle.setVelocity(velocity_);
        }

        if (this->normal.dot(position - this->position) < eEPSILON && this->normal.dot(velocity) < eEPSILON && this->normal.dot(force) < 0) {
            Eigen::Vector3f forceC = -this->normal.dot(force) * this->normal;
            Eigen::Vector3f forceF = -coefFriction * -this->normal.dot(force) * velocityT;
            particle.addForce(forceC + forceF);
        }
    }
    
}
// BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data. 
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    int particleNum = jelly.getParticleNum();
    for (int i = 0; i < particleNum; i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f position = particle.getPosition();
        Eigen::Vector3f velocity = particle.getVelocity();
        Eigen::Vector3f force = particle.getForce();

        Eigen::Vector3f dist = position - this->position;
        if (dist.norm() - this->radius > eEPSILON || this->radius - dist.norm() > eEPSILON || position[1] > eEPSILON)
            continue;

        Eigen::Vector3f normal = (this->position - position).normalized();
        Eigen::Vector3f velocityN = normal.dot(velocity) * normal;
        Eigen::Vector3f velocityT = velocity - velocityN;


        if (normal.dot(position - this->position) < eEPSILON && normal.dot(velocity) < 0) {
            Eigen::Vector3f velocity_ = velocityT - coefResist * velocityN;
            particle.setVelocity(velocity_);
        }

        if (normal.dot(position - this->position) < eEPSILON && normal.dot(velocity) < eEPSILON && normal.dot(force) < 0) {
            Eigen::Vector3f forceC = -normal.dot(force) * normal;
            Eigen::Vector3f forceF = -coefFriction * -normal.dot(force) * velocityT;
            particle.addForce(forceC + forceF);
        }
    }
   
}
}  // namespace simulation
