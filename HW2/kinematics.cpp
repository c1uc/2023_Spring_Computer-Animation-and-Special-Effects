#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
	// TODO#1 (FK)
	// You should set these variables:
	//     bone->start_position = Eigen::Vector4d::Zero();
	//     bone->end_position = Eigen::Vector4d::Zero();
	//     bone->rotation = Eigen::Matrix4d::Zero();
	// The sample above just set everything to zero
	// Hint:
	//   1. posture.bone_translations, posture.bone_rotations
	// Note:
	//   1. This function will be called with bone == root bone of the skeleton
	//   2. we use 4D vector to represent 3D vector, so keep the last dimension as "0"
	//   3. util::rotate{Degree | Radian} {XYZ | ZYX}
	//      e.g. rotateDegreeXYZ(x, y, z) means:
	//      x, y and z are presented in degree rotate z degrees along z - axis first, then y degrees along y - axis, and then x degrees along x - axis 
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

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int allframe_old, int allframe_new) {

	int total_frames = static_cast<int>(postures.size());
	int total_bones = static_cast<int>(postures[0].bone_rotations.size());
	std::vector<acclaim::Posture> new_postures;
	for (int i = 0; i <= allframe_new; ++i) {
		acclaim::Posture new_poseture(total_bones);
		for (int j = 0; j < total_bones; ++j) {

			// TODO#2 (Time warping)
			// original: |--------------|
			// new     : |----------------------|
			// OR
			// original: |--------------|
			// new     : |-------|
			// You should set these variables:
			//     new_postures[i].bone_translations[j] = Eigen::Vector4d::Zero();
			//     new_postures[i].bone_rotations[j] = Eigen::Vector4d::Zero();
			// The sample above just set everything to zero
			// Hint:
			//   1. Scale the frames.
			//   2. You can use linear interpolation with translations.
			//   3. You should use spherical linear interpolation for rotations.

			float ratio = (float)allframe_old / (float)allframe_new;
            int idx_a = std::floor(ratio * (float)i), idx_b = std::ceil(ratio * (float)i);
            float t = ratio * (float)i - idx_a;

            Eigen::Quaterniond rot_a(postures[idx_a].bone_rotations[j]);
            Eigen::Quaterniond rot_b(postures[idx_b].bone_rotations[j]);
            Eigen::Quaterniond rot = rot_a.slerp(t, rot_b);

            Eigen::Vector4d new_rot(rot.x(), rot.y(), rot.z(), 0);

            new_poseture.bone_rotations[j] = new_rot;
            new_poseture.bone_translations[j] =
                (1 - t) * postures[idx_a].bone_translations[j] + t * postures[idx_b].bone_translations[j];
		}

		new_postures.push_back(new_poseture);
	}
	return new_postures;
}
}  // namespace kinematics
