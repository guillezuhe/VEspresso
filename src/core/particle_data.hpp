/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CORE_PARTICLE_DATA_HPP
#define CORE_PARTICLE_DATA_HPP
/** \file
 *  Particles property access.
 *
 *  This file contains everything related to particle properties. If you want to
 *  add a new property to the particles, it is probably a good idea to modify
 *  Particle to give scripts access to that property. You always have to modify
 *  two positions: first the print section, where you should add your new
 *  data at the end, and second the read section where you have to find a nice
 *  and short name for your property to appear in the Python code.
 */

#include "config.hpp"

#include "Particle.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <vector>

/** Call only on the head node: set particle velocity.
 *  @param part the particle.
 *  @param v its new velocity.
 */
void set_particle_v(int part, Utils::Vector3d const &v);

/** Call only on the head node: set particle Lees-Edwards offset.
 *  @param part the particle.
 *  @param v new value for Lees-Edwards offset
 */
void set_particle_lees_edwards_offset(int part, double v);

#ifdef ENGINE
/** Call only on the head node: set particle velocity.
 *  @param part the particle.
 *  @param swim struct containing swimming parameters
 */
void set_particle_swimming(int part, ParticleParametersSwimming swim);
#endif

/** Call only on the head node: set particle force.
 *  @param part the particle.
 *  @param F its new force.
 */
void set_particle_f(int part, const Utils::Vector3d &F);

/** Call only on the head node: set particle mass.
 *  @param part the particle.
 *  @param mass its new mass.
 */
void set_particle_mass(int part, double mass);

#ifdef ROTATIONAL_INERTIA
/** Call only on the head node: set particle rotational inertia.
 *  @param part the particle.
 *  @param rinertia its new inertia.
 */
void set_particle_rotational_inertia(int part, Utils::Vector3d const &rinertia);
#endif

/** Call only on the head node: Specifies whether a particle's rotational
 *  degrees of freedom are integrated or not. If set to zero, the content of
 *  the torque and omega variables are meaningless
 *  @param part the particle.
 *  @param flag the degrees of freedom flag.
 */
void set_particle_rotation(int part, Utils::Vector3i const &flag);

/** @brief rotate a particle around an axis
 *
 *  @param part particle id
 *  @param axis rotation axis
 *  @param angle rotation angle
 */
void rotate_particle(int part, const Utils::Vector3d &axis, double angle);

/** Call only on the head node: set particle charge.
 *  @param part the particle.
 *  @param q its new charge.
 */
void set_particle_q(int part, double q);

#ifdef LB_ELECTROHYDRODYNAMICS
/** Call only on the head node: set particle electrophoretic mobility.
 *  @param part the particle.
 *  @param mu_E its new mobility.
 */
void set_particle_mu_E(int part, Utils::Vector3d const &mu_E);
#endif

/** Call only on the head node: set particle type.
 *  @param p_id the particle.
 *  @param type its new type.
 */
void set_particle_type(int p_id, int type);

/** Call only on the head node: set particle's molecule id.
 *  @param part the particle.
 *  @param mid  its new mol id.
 */
void set_particle_mol_id(int part, int mid);

#ifdef ROTATION
/** Call only on the head node: set particle orientation using quaternions.
 *  @param part the particle.
 *  @param quat its new value for quaternions.
 */
void set_particle_quat(int part, Utils::Quaternion<double> const &quat);

/** Call only on the head node: set particle orientation using director.
 *  The particle director defines the z-axis in the body-fixed frame.
 *  @param part the particle.
 *  @param director its new director vector (will be normalized if necessary)
 */
void set_particle_director(int part, const Utils::Vector3d &director);

/** Call only on the head node: set particle angular velocity from lab frame.
 *  @param part the particle.
 *  @param omega_lab its new angular velocity.
 */
void set_particle_omega_lab(int part, const Utils::Vector3d &omega_lab);

/** Call only on the head node: set particle angular velocity in body frame.
 *  @param part the particle.
 *  @param omega its new angular velocity.
 */
void set_particle_omega_body(int part, const Utils::Vector3d &omega);

/** Call only on the head node: set particle torque from lab frame.
 *  @param part the particle.
 *  @param torque_lab its new torque.
 */
void set_particle_torque_lab(int part, const Utils::Vector3d &torque_lab);
#endif

#ifdef DIPOLES
/** Call only on the head node: set particle dipole orientation.
 *  @param part the particle.
 *  @param dip its new dipole orientation.
 */
void set_particle_dip(int part, Utils::Vector3d const &dip);

/** Call only on the head node: set particle dipole moment (absolute value).
 *  @param part the particle.
 *  @param dipm its new dipole moment.
 */
void set_particle_dipm(int part, double dipm);
#endif

#ifdef VIRTUAL_SITES
/** Call only on the head node: set particle virtual flag.
 *  @param part the particle.
 *  @param is_virtual new @ref ParticleProperties::is_virtual "is_virtual" flag.
 */
void set_particle_virtual(int part, bool is_virtual);
#endif
#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part,
                          Utils::Quaternion<double> const &vs_relative_quat);
void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                              Utils::Quaternion<double> const &rel_ori);
#endif

#ifdef THERMOSTAT_PER_PARTICLE
/** Call only on the head node: set particle frictional coefficient.
 *  @param part the particle.
 *  @param gamma its new frictional coefficient.
 */
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma(int part, double gamma);
#else
void set_particle_gamma(int part, Utils::Vector3d const &gamma);
#endif
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma_rot(int part, double gamma);
#else
void set_particle_gamma_rot(int part, Utils::Vector3d const &gamma_rot);
#endif
#endif
#endif // THERMOSTAT_PER_PARTICLE

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
/** Call only on the head node: set particle external torque.
 *  @param part  the particle.
 *  @param torque new value for ext_torque.
 */
void set_particle_ext_torque(int part, const Utils::Vector3d &torque);
/** Call only on the head node: set particle viscoelastic torque.
 *  @param part  the particle.
 *  @param torque new value for visc_torque.
 */
void set_particle_visc_torque(int part, const Utils::Vector3d &torque);
/** Call only on the head node: set particle viscoelastic torque matrix.
 *  @param part  the particle.
 *  @param torque_mat new value for visc_torque_mat.
 */
void set_particle_visc_torque_mat(int part, const Utils::Matrix<double,20,3> &torque_mat);

void set_particle_omegacrit(int part, const std::vector<double> &omegacrit);
#endif
/** Call only on the head node: set particle external force.
 *  @param part  the particle.
 *  @param force new value for ext_force.
 */
void set_particle_ext_force(int part, const Utils::Vector3d &force);
/** Call only on the head node: set particle viscoelastic force.
 *  @param part  the particle.
 *  @param force new value for visc_force.
 */
void set_particle_visc_force(int part, const Utils::Vector3d &force);
/** Call only on the head node: set particle viscoelastic force matrix.
 *  @param part  the particle.
 *  @param force_mat new value for visc_force_mat.
 */
void set_particle_visc_force_mat(int part, const Utils::Matrix<double,20,3> &force_mat);

void set_particle_visc_gamma_vec(int part, const Utils::Vector3d &visc_gamma_vec);

void set_particle_visc_gamma(int part, const std::vector<double> &visc_gamma);

void set_particle_taum(int part, const std::vector<double> &taum);

void set_particle_vcrit(int part, const std::vector<double> &vcrit);

void set_particle_aexp(int part, const std::vector<double> &aexp);

void set_particle_bexp(int part, const std::vector<double> &bexp);

void set_particle_Nm(int part, int Nm);

/** Call only on the head node: set coordinate axes for which the particles
 *  motion is fixed.
 *  @param part  the particle.
 *  @param flag  coordinates to be fixed.
 */
/** Call only on the head node: set coordinate axes for which the particles
 *  motion is fixed.
 *  @param part  the particle.
 *  @param flag  coordinates to be fixed.
 */
void set_particle_fix(int part, Utils::Vector3i const &flag);
#endif // EXTERNAL_FORCES

/** Call only on the head node: remove bond from particle.
 *  @param part     identity of principal atom of the bond.
 *  @param bond     field containing the bond type number and the identity
 *                  of all bond partners (secondary atoms of the bond).
 */
void delete_particle_bond(int part, Utils::Span<const int> bond);

/** Call only on the head node: remove all bonds from particle.
 *  @param part     identity of principal atom of the bond.
 */
void delete_particle_bonds(int part);

/** @brief Remove the specified bond from the particle
 *  @param p        The particle to update
 *  @param bond     The bond in the form
 *                  <tt>{bond_id, partner_1, partner_2, ...}</tt>
 */
void local_remove_bond(Particle &p, std::vector<int> const &bond);

/** @brief Remove all pair bonds on the particle which have the specified
 *  particle id as partner.
 *  @param p         The particle to update
 *  @param other_pid The particle id to filter for
 */
void local_remove_pair_bonds_to(Particle &p, int other_pid);

/** Call only on the head node: Add bond to particle.
 *  @param part     identity of principal atom of the bond.
 *  @param bond     field containing the bond type number and the identity
 *                  of all bond partners (secondary atoms of the bond).
 */
void add_particle_bond(int part, Utils::Span<const int> bond);

const std::vector<BondView> &get_particle_bonds(int part);

#ifdef EXCLUSIONS
/** @brief Remove particle exclusion.
 *  Call only on the head node.
 *  @param part1    identity of particle for which the exclusion is set.
 *  @param part2    identity of particle for which the exclusion is set.
 */
void remove_particle_exclusion(int part1, int part2);

/** @brief Add particle exclusion.
 *  Call only on the head node.
 *  @param part1    identity of particle for which the exclusion is set.
 *  @param part2    identity of particle for which the exclusion is set.
 */
void add_particle_exclusion(int part1, int part2);

/** Automatically add the next \<distance\> neighbors in each molecule to the
 *  exclusion list.
 *  This uses the bond topology obtained directly from the particles.
 *  To easily setup the bonds, all data should be on a single node,
 *  therefore the \ref partCfg array is used. With large amounts of
 *  particles, you should avoid this function and setup exclusions manually.
 */
void auto_exclusions(int distance);
#endif

/** Rescale all particle positions in direction @p dir by a factor @p scale. */
void mpi_rescale_particles(int dir, double scale);

// The following functions are used by the python interface to obtain
// properties of a particle, which are only compiled in in some configurations
// This is needed, because cython does not support conditional compilation
// within a ctypedef definition

#ifdef VIRTUAL_SITES_RELATIVE
inline Utils::Quaternion<double> get_particle_vs_quat(Particle const *p) {
  return p->vs_relative().quat;
}
inline Utils::Quaternion<double> get_particle_vs_relative(Particle const *p,
                                                          int &vs_relative_to,
                                                          double &vs_distance) {
  vs_relative_to = p->vs_relative().to_particle_id;
  vs_distance = p->vs_relative().distance;
  return p->vs_relative().rel_orientation;
}
#endif

#ifdef EXTERNAL_FORCES
inline Utils::Vector3i get_particle_fix(Particle const *p) {
  return Utils::Vector3i{
      {p->is_fixed_along(0), p->is_fixed_along(1), p->is_fixed_along(2)}};
}
#endif // EXTERNAL_FORCES

#ifdef THERMOSTAT_PER_PARTICLE
#ifdef PARTICLE_ANISOTROPY
inline Utils::Vector3d get_particle_gamma(Particle const *p) {
  return p->gamma();
}
#else
inline double get_particle_gamma(Particle const *p) { return p->gamma(); }
#endif // PARTICLE_ANISOTROPY
#ifdef ROTATION
#ifdef PARTICLE_ANISOTROPY
inline Utils::Vector3d get_particle_gamma_rot(Particle const *p) {
  return p->gamma_rot();
}
#else
inline double get_particle_gamma_rot(Particle const *p) {
  return p->gamma_rot();
}
#endif // PARTICLE_ANISOTROPY
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE

#ifdef ROTATION
inline Utils::Vector3i get_particle_rotation(Particle const *p) {
  return Utils::Vector3i{{p->can_rotate_around(0), p->can_rotate_around(1),
                          p->can_rotate_around(2)}};
}
#endif // ROTATION

#endif
