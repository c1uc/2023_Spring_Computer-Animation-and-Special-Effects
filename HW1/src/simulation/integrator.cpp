#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.15 - p.16

    float timeDelta = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();
    for (int i = 0; i < jellyCount; i++) {
        Jelly *jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);
            
            Eigen::Vector3f pos = particle.getPosition();
            Eigen::Vector3f vel = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + posDelta);
            particle.setVelocity(vel + velDelta);
            particle.setForce(force);
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review ¡§ODE_basics.pptx¡¨ from p.18 - p.19

    std::vector<Jelly> jelliesBefore(particleSystem.jellies);

    float timeDelta = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            particle.setForce(force);
        }
    }

    particleSystem.computeAllForce();
    std::vector<Jelly> jelliesAfter(particleSystem.jellies);

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();

            Eigen::Vector3f vel_ = jelliesAfter[i].getParticle(j).getVelocity();
            Eigen::Vector3f acc = jelliesAfter[i].getParticle(j).getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel_;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + posDelta);
            particle.setVelocity(vel + velDelta);
            particle.setForce(force);
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review ¡§ODE_basics.pptx¡¨ from p .18 - p .20

    std::vector<Jelly> jelliesBefore(particleSystem.jellies);

    float timeDelta = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();
    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = particle.getPosition();
            Eigen::Vector3f vel = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = 0.5f * timeDelta * vel;
            Eigen::Vector3f velDelta = 0.5f * timeDelta * acc;

            particle.setPosition(pos + posDelta);
            particle.setVelocity(vel + velDelta);
            particle.setForce(force);
        }
    }

    particleSystem.computeAllForce();
    std::vector<Jelly> jelliesAfter(particleSystem.jellies);

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();

            Eigen::Vector3f vel_ = jelliesAfter[i].getParticle(j).getVelocity();
            Eigen::Vector3f acc = jelliesAfter[i].getParticle(j).getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel_;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + posDelta);
            particle.setVelocity(vel + velDelta);
            particle.setForce(force);
        }
    }

}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.21
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };

    std::vector<Jelly> k1, k2, k3, k4;

    k1 = particleSystem.jellies;

    std::vector<Jelly> jelliesBefore(particleSystem.jellies);

    float timeDelta = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();
            Eigen::Vector3f acc = jelliesBefore[i].getParticle(j).getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + 0.5 * posDelta);
            particle.setVelocity(vel + 0.5 * velDelta);
            particle.setForce(force);
        }
    }

    particleSystem.computeAllForce();
    k2 = particleSystem.jellies;

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();

            Eigen::Vector3f vel_ = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel_;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + 0.5 * posDelta);
            particle.setVelocity(vel + 0.5 * velDelta);
            particle.setForce(force);
        }
    }

    particleSystem.computeAllForce();
    k3 = particleSystem.jellies;

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();

            Eigen::Vector3f vel_ = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * vel_;
            Eigen::Vector3f velDelta = timeDelta * acc;

            particle.setPosition(pos + 0.5 * posDelta);
            particle.setVelocity(vel + 0.5 * velDelta);
            particle.setForce(force);
        }
    }

    particleSystem.computeAllForce();
    k4 = particleSystem.jellies;

    for (int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleCount = jelly->getParticleNum();
        for (int j = 0; j < particleCount; j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f pos = jelliesBefore[i].getParticle(j).getPosition();
            Eigen::Vector3f vel = jelliesBefore[i].getParticle(j).getVelocity();
            Eigen::Vector3f force = Eigen::Vector3f::Zero();

            Eigen::Vector3f posDelta = timeDelta * (k1[i].getParticle(j).getVelocity() + 2 * k2[i].getParticle(j).getVelocity() +
                                        2 * k3[i].getParticle(j).getVelocity() + k4[i].getParticle(j).getVelocity()) /
                                       6;
            Eigen::Vector3f velDelta =
                timeDelta * 
                (k1[i].getParticle(j).getAcceleration() + 2 * k2[i].getParticle(j).getAcceleration() +
                 2 * k3[i].getParticle(j).getAcceleration() + k4[i].getParticle(j).getAcceleration()) /
                                       6;

            particle.setPosition(pos + posDelta);
            particle.setVelocity(vel + velDelta);
            particle.setForce(force);
        }
    }

}
}  // namespace simulation
