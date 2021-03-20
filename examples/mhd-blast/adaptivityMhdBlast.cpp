#include "adaptivityMhdBlast.h"

template <int dim>
AdaptivityMhdBlast<dim>::AdaptivityMhdBlast(Parameters<dim>& parameters, MPI_Comm& mpi_communicator) :
  Adaptivity<dim>(parameters, mpi_communicator),
  last_time_step(0),
  adaptivity_step(0)
{
}

template <int dim>
void AdaptivityMhdBlast<dim>::calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator)
{
  FEValuesExtractors::Scalar scalars[2];
  scalars[0].component = 0;
  scalars[1].component = 4;
  const QGauss<dim - 1> face_quadrature(1);
  UpdateFlags face_update_flags = UpdateFlags(update_values | update_JxW_values);
  FEFaceValues<dim> fe_v_face(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, dof_handler.get_fe(), face_quadrature, update_values);
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;
    Point<dim> jump;
    Point<dim> area;
    for (unsigned int face_no = 0; face_no < 4; ++face_no)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
      if (!face->at_boundary())
      {
        Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
        typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
        std::vector<double> u(fe_v_face.n_quadrature_points);
        std::vector<double> u_neighbor(fe_v_face.n_quadrature_points);
        if (face->has_children())
        {
          unsigned int neighbor2 = cell->neighbor_face_no(face_no);
          for (unsigned int subface_no = 0; subface_no < face->number_of_children(); ++subface_no)
          {
            typename DoFHandler<dim>::cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);
            Assert(!neighbor_child->has_children(), ExcInternalError());
            fe_v_subface.reinit(cell, face_no, subface_no);
            fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
            const std::vector<double> &JxW = fe_v_subface.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_subface.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_subface[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_subface.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
        }
        else
        {
          if (!cell->neighbor_is_coarser(face_no))
          {
            unsigned int neighbor2 = cell->neighbor_of_neighbor(face_no);
            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, neighbor2);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
          else //i.e. neighbor is coarser than cell
          {
            std::pair<unsigned int, unsigned int> neighbor_face_subface
              = cell->neighbor_of_coarser_neighbor(face_no);
            Assert(neighbor_face_subface.first < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
            Assert(neighbor_face_subface.second < neighbor->face(neighbor_face_subface.first)->number_of_children(),
              ExcInternalError());
            Assert(neighbor->neighbor_child_on_subface(neighbor_face_subface.first, neighbor_face_subface.second)
              == cell, ExcInternalError());
            fe_v_face.reinit(cell, face_no);
            fe_v_subface.reinit(neighbor, neighbor_face_subface.first,
              neighbor_face_subface.second);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_subface[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
        }
      }
    }
    double average_jumps[dim];
    double sum_of_average_jumps = 0.;
    for (unsigned int i = 0; i < 2; ++i)
    {
      average_jumps[i] = jump(i) / area(i);
      sum_of_average_jumps += average_jumps[i];
    }
    for (int i = 0; i < this->parameters.volume_factor; i++)
      sum_of_average_jumps *= cell->diameter();
    gradient_indicator(cell->active_cell_index()) = sum_of_average_jumps;
  }
  for (int i = 0; i < gradient_indicator.size(); i++)
    if (gradient_indicator[i] < SMALL)
      gradient_indicator[i] = 0.;
}

template <int dim>
bool AdaptivityMhdBlast<dim>::refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
  , const Mapping<dim>& mapping)
{
	//std::cout << "parameter1: " << time_step % this->parameters.refine_every_nth_time_step;
	//std::cout << " parameter2: " << (this->parameters.perform_n_initial_refinements : 1);
  if (time_step % this->parameters.refine_every_nth_time_step)
    return false;
  if (++adaptivity_step > (time_step == 0 ? this->parameters.perform_n_initial_refinements : 1))
  {
    adaptivity_step = 0;
    return false;
  }
  Vector<double> gradient_indicator(triangulation.n_active_cells());

  calculate_jumps(solution, dof_handler, mapping, gradient_indicator);
  

  int max_calls_ = this->parameters.max_cells + (int)std::floor(time * this->parameters.max_cells * this->parameters.time_interval_max_cells_multiplicator / this->parameters.final_time);
  GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, gradient_indicator, this->parameters.refine_threshold, this->parameters.coarsen_threshold, max_calls_);
  
//The following is tricky but usefull for anisotropic mesh refinement 
//CANT WORK WITH distributed MESH !! (fundamental dealii-p4est limitation) 
 //aditionally, mesh smoothing in declaration of Triangulation class must be properly set
#ifdef HAVE_MPI
 // if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)  std::cout << "smoothing: " << triangulation.get_mesh_smoothing() << "/n";
 // triangulation.prepare_coarsening_and_refinement();//!!!
#else
  triangulation.prepare_coarsening_and_refinement();
  
#endif
  



  return true;






}

template class AdaptivityMhdBlast<2>;
