#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO (FK)
    // Same as HW2
    // Hint:
    //   1. If you don't use `axis` in this function, you can copy-paste your code
    if (bone->idx == 0)
        bone->start_position = posture.bone_translations[0], bone->end_position = posture.bone_translations[0],
        bone->rotation = util::rotateDegreeZYX(posture.bone_rotations[0]).toRotationMatrix();

    auto chi = bone->child;
    while (chi != nullptr) {
        chi->start_position = bone->end_position;
        chi->rotation = bone->rotation * chi->rot_parent_current *
                        util::rotateDegreeZYX(posture.bone_rotations[chi->idx]).toRotationMatrix();
        chi->end_position =
            chi->start_position + (chi->rotation * (chi->length * chi->dir) + posture.bone_translations[chi->idx]);
        forwardSolver(posture, chi);
        chi = chi->sibling;
    }
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::VectorXd deltatheta(Jacobian.cols());
    deltatheta.setZero();

    deltatheta = Jacobian.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(target);

    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    std::vector<acclaim::Bone*> boneList;
    // TODO
    // Calculate number of bones need to move to perform IK, store in `bone_num` 
    // a.k.a. how may bones from end_bone to its parent then to start_bone (include both start_bone and end_bone)
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start_bone and end_bone should be in the list
    acclaim::Bone* current = end_bone;
    
    while (current != start_bone) {
        boneList.emplace_back(current);
        bone_num++;
        current = current->parent;

        if (current == nullptr) {
            current = start_bone;
            while (current != root_bone) {
                boneList.emplace_back(current);
                bone_num++;
                current = current->parent;
            }
            break;
        }
    }
    if (current == start_bone)
        boneList.emplace_back(current);
    
    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        if (desiredVector.norm() < epsilon) {
            break;
        }
        // TODO (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */
        for (long long i = 0; i < bone_num; i++) {
            auto b = boneList[i];
            auto rotation = b->rotation.matrix();
            auto rot_ = rotation.block<3, 3>(0, 0);
            Eigen::Vector4d r = end_bone->end_position - b->start_position;
            Eigen::Vector3d r_;
            r_ << r[0], r[1], r[2];
            if (b->dofrx) {
                Eigen::Vector3d v = rot_.col(0);
                auto w = v.cross(r_);
                Jacobian.col(3 * i) << w, 0;
            }
            if (b->dofry) {
                Eigen::Vector3d v = rot_.col(1);
                auto w = v.cross(r_);
                Jacobian.col(3 * i + 1) << w, 0;
            }
            if (b->dofrz) {
                Eigen::Vector3d v = rot_.col(2);
                auto w = v.cross(r_);
                Jacobian.col(3 * i + 2) << w, 0;
            }
        }
        
        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        
        // TODO (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.
        
        for (long long i = 0; i < bone_num; i++) {
            for (int j = 0; j < 3; j++) {
                posture.bone_rotations[boneList[i]->idx][j] += deltatheta[3 * i + j] * 180 / util::PI;
            }
            posture.bone_rotations[boneList[i]->idx][0] = std::clamp(
                posture.bone_rotations[boneList[i]->idx][0], (double)boneList[i]->rxmin, (double)boneList[i]->rxmax);
            posture.bone_rotations[boneList[i]->idx][1] = std::clamp(
                posture.bone_rotations[boneList[i]->idx][1], (double)boneList[i]->rymin, (double)boneList[i]->rymax);
            posture.bone_rotations[boneList[i]->idx][2] = std::clamp(
                posture.bone_rotations[boneList[i]->idx][2], (double)boneList[i]->rzmin, (double)boneList[i]->rzmax);
        }

    }
    // TODO (Bonus)
    // Return whether IK is stable (i.e. whether the ball is reachable) and let the skeleton not swing its hand in the air
    if ((target_pos - end_bone->end_position).norm() < epsilon) {
        return true;
    } else {
        posture = original_posture;
        return false;
    }
}
}  // namespace kinematics
