#include "postProcessor.h"

using namespace dealii;

Postprocessor::Postprocessor() : DataPostprocessor<DIM>()
{}

void Postprocessor::compute_derived_quantities_scalar(const std::vector<d>         &uh0,
    const std::vector<Tensor<1, DIM> > &duh0,
    const std::vector<Tensor<2, DIM> > &dduh0,
    const std::vector<Point<DIM> >    &normals0,
    const std::vector<Point<DIM> >    &evaluation_points0,
    std::vector<Vector<double> >      &computed_quantities) const
{
    const ui n_quadrature_points = uh0.size();

    for (ui q = 0; q < n_quadrature_points; ++q)
    {
       computed_quantities[q](0) = uh0[q];
    }
}

void Postprocessor::compute_derived_quantities_vector(const std::vector<Vector<d> > &uh,
  const std::vector<vecDimVec> &duh,
  const std::vector<std::vector<Tensor<2, DIM> > > &dduh,
  const std::vector<Point<DIM> > &normals,
  const std::vector<Point<DIM> > &evaluation_points, 
  std::vector<Vector<d> > &computed_quantities) const
{
  const ui n_quadrature_points = uh.size();

  for (ui q = 0; q < n_quadrature_points; ++q)
  {
    for (ui d = 0; d < COMPONENT_COUNT; ++d)
      computed_quantities[q](d) = uh[q](d);
  }
}

std::vector<std::string> Postprocessor::get_names() const
{
  std::vector<std::string> names;
  for (ui d = 0; d < COMPONENT_COUNT; ++d)
  {
    std::stringstream ss;
    switch(d)
    {
      case 0:
      ss << "Rho";
      break;
    case 1:
      ss << "Pi_x";
      break;
    case 2:
      ss << "Pi_Y";
      break;
    case 3:
      ss << "Pi_z";
      break;
    case 4:
      ss << "B_x";
      break;
    case 5:
      ss << "B_y";
      break;
    case 6:
      ss << "B_z";
      break;
    case 7:
      ss << "U";
      break;
    case 8:
      ss << "J_x";
      break;
    case 9:
      ss << "J_y";
      break;
    case 10:
      ss << "J_z";
      break;
    }

    names.push_back(ss.str());
  }
  return names;
}

std::vector<DataComponentInterpretation::DataComponentInterpretation> Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
  for (ui d = 0; d < COMPONENT_COUNT; ++d)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  return interpretation;
}

UpdateFlags Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients | update_quadrature_points;
}
