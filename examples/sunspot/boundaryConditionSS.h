#ifndef _BOUNDARY_CONDITION_SS_H
#define _BOUNDARY_CONDITION_SS_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryCondition.h"
#include "equationsMhd.h"
#include "initialConditionSS.h"
#include "parametersSS.h"

template <int dim>
class BoundaryConditionSSWithVortices : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionSSWithVortices(Parameters<dim>&, SSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  SSParameters& cs_parameters;
  double eps;
  double y_1, y_2;
  double r_1_bar(double x, double y) const;
  double r_2_bar(double x, double y) const;
  double omega_1(double x, double y) const;
  double omega_2(double x, double y) const;
  double omega(double time) const;
};

template <int dim>
class BoundaryConditionSSFree : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionSSFree(Parameters<dim>&, SSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
};

template <int dim>
class BoundaryConditionSSTest : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionSSTest(Parameters<dim>&, SSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
private:
  SSParameters& cs_parameters;
};

template <int dim>
class BoundaryConditionSSInitialState : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionSSInitialState(Parameters<dim>&, SSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result,
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  SSParameters& cs_parameters;
  InitialConditionSS<EquationsTypeMhd, dim> ic;
  double invL_G;
  double iSgn;
  double d2R;
  double H;
  double L2R;
  double R2L;
  double q_mag;
};
#endif
