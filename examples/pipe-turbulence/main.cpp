#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "parametersCS.h"
#include "initialConditionCS.h"
#include "boundaryConditionCS.h"
#include "adaptivityCS.h"

// Dimension of the problem - passed as a template parameter to pretty much every class.
#define DIMENSION 2
// Type of equations, must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd

#ifdef HAVE_MPI
void set_triangulation(parallel::distributed::Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#else
void set_triangulation(Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#endif
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, parameters.refinements, parameters.corner_a, parameters.corner_b, true);
}

void set_parameters(Parameters<DIMENSION>& parameters, CSParameters& cs_parameters)
{
  parameters.slope_limiter = parameters.vertexBased;
  //parameters.corner_a = Point<DIMENSION>(-5, -10., 0.);
  //parameters.corner_b = Point<DIMENSION>(5., 10., 0.5);
  //parameters.refinements = { 50, 100 , 1 };//ok je na oase napr. 300:9000-15cpu na nod pri 8 st. vol. na nod.
  parameters.corner_a = Point<DIMENSION>(-1., -10.);
  parameters.corner_b = Point<DIMENSION>(1., 10.);
  parameters.refinements = { 1000, 10000 };//ok je na oase napr. 300:9000-15cpu na nod pri 8 st. vol. na nod.

  parameters.limit = true;
  parameters.limitB = true;
  parameters.use_div_free_space_for_B = false;
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .01;
  parameters.start_limiting_at = -1e-6;//e-6
  parameters.quadrature_order = 1;
  parameters.polynomial_order_dg = 0;  
  parameters.patches = 0;
  parameters.output_step = 0.05;
  parameters.final_time = 10.;
  parameters.output_file_prefix = "solution";

  parameters.max_cells = 100000;
  parameters.refine_every_nth_time_step = 100;
  parameters.perform_n_initial_refinements = 2;//15
  parameters.refine_threshold = 0.5;
  parameters.coarsen_threshold = 0.2;
  parameters.volume_factor = 4;
  parameters.time_interval_max_cells_multiplicator = 0.;
  parameters.gas_gamma = 1.25;
  // plasma beta
  cs_parameters.beta = 0.10;
  

  // coronal height scale
  cs_parameters.L_G = 0.;

  // Gravity acceleration
  // g = ( l_0 / v^2_0 ) g
  if (cs_parameters.L_G > NEGLIGIBLE)
  {
    double l_0 = 1.2e8 / cs_parameters.L_G;
    double t_0 = 10.;
    double v_0 = l_0 / t_0;
    parameters.g = 0.;
  }

  // Density
  cs_parameters.rho_0 = 1.;

  // Torus winding number
  cs_parameters.N_t = 5.;

  // Torus major radius
  cs_parameters.R = 3.0;

  // Torus minor radius
  cs_parameters.r = 1.0;

  // Magnetic charge separation distance
  cs_parameters.L = 1.5;

  // Geometrical factor
  cs_parameters.d = 1.5;

  // The coronal/prominence temperature ratio
  cs_parameters.Tc2Tp = 1.;

  cs_parameters.omega_0 = 0.3;

  cs_parameters.t_drive = 2.0;

  cs_parameters.t_ramp = 1.0;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    CSParameters cs_parameters;
    set_parameters(parameters, cs_parameters);
    parameters.delete_old_outputs(mpi_communicator);

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::limit_level_difference_at_vertices));
#else
    Triangulation<DIMENSION> triangulation(Triangulation<DIMENSION>::limit_level_difference_at_vertices);
#endif    
    set_triangulation(triangulation, parameters);

    InitialConditionCS<EQUATIONS, DIMENSION> initial_condition(parameters, cs_parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditionSlab<DIMENSION> boundary_conditions(parameters, cs_parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    //AdaptivityCS<DIMENSION> adaptivity(parameters, mpi_communicator);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Set adaptivity
    // problem.set_adaptivity(&adaptivity);
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
