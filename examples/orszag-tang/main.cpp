#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialConditionOT.h"

// Dimension of the problem - passed as  template parameter to pretty much every class.
#define DIMENSION 3
// Type of equations, must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd

#ifdef HAVE_MPI
void set_triangulation(parallel::distributed::Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#else
void set_triangulation(Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#endif
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, parameters.refinements, parameters.corner_a, parameters.corner_b, true);
  
  std::vector<dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<DIMENSION> > > > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
}

void set_parameters(Parameters<DIMENSION>& parameters)  
{
  parameters.slope_limiter = parameters.vertexBased;
  parameters.corner_a = Point<DIMENSION>(0., 0., 0.);
  parameters.corner_b = Point<DIMENSION>(1., 1., 0.001);
  parameters.refinements = { 200, 200, 1 };
  parameters.limit = false;
  parameters.use_div_free_space_for_B = true;
  parameters.periodic_boundaries = { { 0, 1, 0 },{ 2, 3, 1 } };
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .05;
  parameters.start_limiting_at = .05;
  parameters.quadrature_order = 1;
  parameters.polynomial_order_dg = 1;
  parameters.patches = 0;
  parameters.output_step = 2.e-2;
  parameters.final_time = .5;
  parameters.debug = parameters.Adaptivity;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    set_parameters(parameters);
    parameters.delete_old_outputs(mpi_communicator);
    std::cout << std::endl <<
        parameters.limit << std::endl <<
        parameters.limitB << std::endl <<
        parameters.limit_edges_and_vertices << std::endl <<
        parameters.output_file_prefix << std::endl <<
        parameters.use_div_free_space_for_B << std::endl <<

        parameters.lax_friedrich_stabilization_value << std::endl <<
        parameters.cfl_coefficient << std::endl <<
        parameters.quadrature_order << std::endl <<
        parameters.polynomial_order_dg << std::endl <<
        parameters.patches << std::endl <<
        parameters.output_step << std::endl <<
        parameters.final_time << std::endl;

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::none), parallel::distributed::Triangulation<DIMENSION>::no_automatic_repartitioning);
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    set_triangulation(triangulation, parameters);

    InitialConditionOT<EQUATIONS, DIMENSION> initial_condition(parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryCondition<EQUATIONS, DIMENSION> boundary_conditions(parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Run the problem - entire transient problem.
    problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
      << exc.what() << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  };

  return 0;
}
