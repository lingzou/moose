//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFEFluidMomentumBC.h"

registerMooseObject("NavierStokesApp", INSFEFluidMomentumBC);
registerMooseObjectRenamed("NavierStokesApp",
                           MDFluidMomentumBC,
                           "02/01/2024 00:00",
                           INSFEFluidMomentumBC);

InputParameters
INSFEFluidMomentumBC::validParams()
{
  InputParameters params = INSFEFluidIntegratedBCBase::validParams();
  params.addClassDescription("Specifies flow of momentum through a boundary");
  
  MooseEnum bc_type("VelocityInlet Pressure CoupledPressure Wall Invalid", "Invalid");
  params.addParam<MooseEnum>("bc_type", bc_type, "The type of physical boundary condition");

  params.addRequiredParam<unsigned>("component", "0,1,or 2 for x-, y-, or z- direction");
  params.addParam<FunctionName>("p_fn", "Pressure function with time at the boundary");
  params.addParam<FunctionName>("v_fn", "Velocity function with time at the boundary");

  // coupled with branch pressure and density
  // The 'branch_center' is a little bit tricky, because SAM 1D and multi-D could be in
  // different mesh system.
  //   * The volume branch center is always defined in physical 3D XYZ coordinate system,
  //   * but multi-D flow could be simulated in 2D XY coordinate system,
  //   * the volume brance center needs be mapped to the 2D/3D flow mesh system
  // The pressure at the multi-D boundary and the branch pressure is related by:
  //   p_boundary = p_branch + rho_branch * (point_boundary - branch_center) * gravity
  params.addCoupledVar("p_branch", "Coupled scalar branch pressure");
  params.addCoupledVar("rho_branch", "Coupled scalar branch density for gravity head calculation");
  params.addParam<VectorValue<Real>>("gravity", "Gravity vector in 2D/3D flow mesh system");
  params.addParam<Point>("branch_center", "Position of branch center in 2D/3D flow mesh system");
  return params;
}

INSFEFluidMomentumBC::INSFEFluidMomentumBC(const InputParameters & parameters)
  : INSFEFluidIntegratedBCBase(parameters),
    _bc_type(getParam<MooseEnum>("bc_type")),
    _component(getParam<unsigned>("component")),
    _mu(getMaterialProperty<Real>("dynamic_viscosity")),
    _mu_t(getMaterialProperty<Real>("turbulence_viscosity")),
    _has_pbc(parameters.isParamValid("p_fn")),
    _has_vbc(parameters.isParamValid("v_fn")),
    _p_fn(_has_pbc ? &getFunction("p_fn") : nullptr),
    _v_fn(_has_vbc ? &getFunction("v_fn") : nullptr),
    _has_pbranch(parameters.isParamValid("p_branch")),
    _p_branch(_has_pbranch ? coupledScalarValue("p_branch") : _zero),
    _p_branch_var_number(_has_pbranch ? coupledScalar("p_branch") : libMesh::invalid_uint),
    _rho_branch(_has_pbranch ? coupledScalarValue("rho_branch") : _zero)
{
  if (_bc_type == "Invalid")
    mooseError("A valid physical boundary type, bc_type, must be specified for INSFEFluidMomentumBC"
              "Candidates are: 'VelocityInlet', 'Pressure', 'CoupledPressure', or 'Wall'.");
  else if (_bc_type == "VelocityInlet")
  {
    if (!isParamValid("v_fn"))
      mooseError("A velocity function, v_fn, is required for 'VelocityInlet' boundary.");
    if (isParamValid("p_fn") || isParamValid("p_branch"))
      mooseError("A pressure value/function cannot be specified for 'VelocityInlet' boundary.");
  }
  else if (_bc_type == "Pressure")
  {
    if (isParamValid("v_fn"))
      mooseError("A velocity value/function cannot be specified for 'Pressure' boundary.");
    if (!isParamValid("p_fn"))
      mooseError("A pressure function, p_fn, is required for 'Pressure' boundary.");
  }
  else if (_bc_type == "CoupledPressure")
  {
    if (isParamValid("v_fn"))
      mooseError("A velocity value/function cannot be specified for 'CoupledPressure' boundary.");
    if (!(isParamValid("p_branch") && isParamValid("rho_branch") && isParamValid("gravity")))
      mooseError(
          name(),
          ": this boundary is coupled to a volume branch, ",
          "please provide 'p_branch', for pressure at the volume branch,"
          "'gravity' vector and 'branch_center' for gravity head calculation.");
    _vec_g = getParam<VectorValue<Real>>("gravity");
    _branch_center = getParam<Point>("branch_center");
  }
  else if (_bc_type == "Wall")
  {
    if ((isParamValid("v_fn") || isParamValid("p_fn") || isParamValid("p_branch")))
      mooseError("None of 'v_fn', 'p_fn', or 'p_branch' is expected for the 'Wall' boundary.");
  }
  else
  {
    // should never arrive here because of the default "Invalid" value
    mooseError("Code should never reach here."
                  "If you see this line, please contact the development team.");
  }
}

Real
INSFEFluidMomentumBC::computeQpResidual()
{
  Real porosity = _has_porosity ? _porosity[_qp] : 1;
  RealVectorValue vec_vel(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  Real p_bc = 0;
  Real v_dot_n = 0;
  if (_bc_type == "VelocityInlet")
  {
     p_bc = _pressure[_qp];
    v_dot_n = -_v_fn->value(_t, _q_point[_qp]);
  }
  else if (_bc_type == "Pressure")
  {
    p_bc = _p_fn->value(_t, _q_point[_qp]);
      v_dot_n = vec_vel * _normals[_qp];
  }
  else if (_bc_type == "CoupledPressure")
  {
    Real dH = (_q_point[_qp] - _branch_center) * _vec_g;
    p_bc = _p_branch[0] + _rho_branch[0] * dH;
    v_dot_n = vec_vel * _normals[_qp];
  }
  else if (_bc_type == "Wall")
  {
    p_bc = _pressure[_qp];
    v_dot_n = 0;
  }
  else
    mooseError("Unknown types of bc_type in 'INSFEFluidMomentumBC'.");

  Real viscous_part = (porosity > 0.99)
                          ? -(_mu[_qp] + _mu_t[_qp]) * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp]
                          : 0;

  // BC contribution is: eps * p_bc * psi * n_i                   (pressure part)
  //                      + psi * rho * v_i * v_dot_n / eps       (convection part)
  //                      - psi * mu * (grad(u_i) dot n)          (viscous part/conditional)
  return (porosity * p_bc * _normals[_qp](_component) + _rho[_qp] * _u[_qp] * v_dot_n / porosity) *
             _test[_i][_qp] +
         viscous_part;
}

Real
INSFEFluidMomentumBC::computeQpJacobian()
{
  Real porosity = _has_porosity ? _porosity[_qp] : 1;
  RealVectorValue vec_vel(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  Real v_dot_n = 0;
  Real convection_part = 0;
  if (_bc_type == "VelocityInlet")
  {
    v_dot_n = -_v_fn->value(_t, _q_point[_qp]);
    convection_part = _rho[_qp] * _phi[_j][_qp] * v_dot_n / porosity * _test[_i][_qp];
  }
  else if ((_bc_type == "Pressure") || (_bc_type == "CoupledPressure"))
  {
    v_dot_n = vec_vel * _normals[_qp];
    convection_part =
      _rho[_qp] * _phi[_j][_qp] * v_dot_n / porosity * _test[_i][_qp] +
      _rho[_qp] * _u[_qp] * _phi[_j][_qp] * _normals[_qp](_component) / porosity * _test[_i][_qp];
  }
  else if (_bc_type == "Wall")
  {
    // do nothing
  }
  else
    mooseError("Unknown types of bc_type in 'INSFEFluidMomentumBC'.");

  Real viscous_part = (porosity > 0.99)
                          ? -(_mu[_qp] + _mu_t[_qp]) * _grad_phi[_j][_qp](_component) *
                                _normals[_qp](_component) * _test[_i][_qp]
                          : 0;

  return convection_part + viscous_part;
}

Real
INSFEFluidMomentumBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned m = this->mapVarNumber(jvar);
  RealVectorValue vec_vel(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real porosity = _has_porosity ? _porosity[_qp] : 1;

  // this is the jacobian w.r.t branch pressure
  if (jvar == _p_branch_var_number)
  {
    return porosity * _normals[_qp](_component) * _test[_i][_qp];
  }

  // this is the jacobian w.r.t multi-dimensional flow variables
  Real jac = 0;
  switch (m)
  {
    case 0:
      if ((_bc_type == "VelocityInlet") || (_bc_type == "Wall"))
        jac = porosity * _phi[_j][_qp] * _normals[_qp](_component) * _test[_i][_qp];
      break;
    case 1:
    case 2:
    case 3:
    {
      if (m != (_component + 1))
        jac =
            _rho[_qp] / porosity * _phi[_j][_qp] * _test[_i][_qp] * _u[_qp] * _normals[_qp](m - 1);
      break;
    }
    case 4:
    default:
      jac = 0;
  }
  return jac;
}
