#ifndef _ADAPTIVITY_SS_H
#define _ADAPTIVITY_SS_H

#include "adaptivity.h"
#include "parametersSS.h"

// Initial condition
template <int dim>
class AdaptivitySS : public Adaptivity<dim>
{
public:
  AdaptivitySS(Parameters<dim>&, MPI_Comm& mpi_communicator);
  bool refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
    , const Mapping<dim>& mapping);
  
  void calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator);
  int last_time_step;
  int adaptivity_step;

  int max_cells;
  int refine_every_nth_time_step;
  int perform_n_initial_refinements;
  double refine_threshold;
  double coarsen_threshold;
  std::array <unsigned short, BASIS_FN_COUNT> component_ii;
  std::array <bool, BASIS_FN_COUNT> is_primitive;
  FEValuesExtractors::Vector mag;
  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_neighbor;
};

#endif