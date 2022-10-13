#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "parametersTD.h"
#include "initialConditionTD.h"
#include "boundaryConditionTD.h"
#include "adaptivityTD.h"


// Dimension of the problem - passed as a template parameter to pretty much every class.
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
}
std::vector<std::vector<std::string>> inputvector;
void set_parameters(Parameters<DIMENSION>& parameters, TitovDemoulinParameters& td_parameters)
{
  parameters.slope_limiter = parameters.vertexBased;
  parameters.corner_a = Point<DIMENSION>(-5., -10., 0.);
  parameters.corner_b = Point<DIMENSION>(5., 10., 10.);
  parameters.refinements = { 3, 6, 3 };//60 120 60
  parameters.limit = false;
  parameters.limitB = false;
  parameters.use_div_free_space_for_B = false;
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .01;
  parameters.start_limiting_at = -1e-6;
  parameters.quadrature_order = 2;
  parameters.polynomial_order_dg = 1;
  parameters.patches = 0;
  parameters.output_step = 5.e-1;
  parameters.final_time = 0.1;//12.
  parameters.output_file_prefix = "vtx_limit";

  parameters.max_cells = 2500;
  parameters.refine_every_nth_time_step = 25;
  parameters.perform_n_initial_refinements = 0;//25
  parameters.refine_threshold = 0.5;
  parameters.coarsen_threshold = 0.2;
  parameters.volume_factor = 4;
  parameters.time_interval_max_cells_multiplicator = 0.;

  // plasma beta
  td_parameters.beta = 0.05;

  // coronal height scale
  td_parameters.L_G = 20.;

  // Gravity acceleration
  // g = ( l_0 / v^2_0 ) g
  if (td_parameters.L_G > NEGLIGIBLE)
  {
    double l_0 = 1.2e8 / td_parameters.L_G;
    double t_0 = 10.;
    double v_0 = l_0 / t_0;
    parameters.g = (l_0 / (v_0 * v_0)) * 274.;
  }

  // Density
  td_parameters.rho_0 = 1.;

  // Torus winding number
  td_parameters.N_t = 5.;

  // Torus major radius
  td_parameters.R = 3.0;

  // Torus minor radius
  td_parameters.r = 1.0;

  // Magnetic charge separation distance
  td_parameters.L = 1.5;

  // Geometrical factor
  td_parameters.d = 1.5;

  // The coronal/prominence temperature ratio
  td_parameters.Tc2Tp = 1.;

  td_parameters.omega_0 = 0.3;

  td_parameters.t_drive = 2.0;

  td_parameters.t_ramp = 1.0;
}
 
 
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {  
      
      std::vector<std::string> row;
      std::string line, word;

      std::fstream input("/r/2D.csv", std::ios::in);
          if (input.is_open())
          {
              while (getline(input, line))
              {
                  row.clear();
                  std::stringstream str(line);
                  while (getline(str, word, ';'))
                      row.push_back(word);
                  inputvector.push_back(row);
              }
          }
          else std::cout << "Could not open input file";

    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    TitovDemoulinParameters td_parameters;
    set_parameters(parameters, td_parameters);
    parameters.delete_old_outputs(mpi_communicator);
    
    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::limit_level_difference_at_vertices));
#else
    Triangulation<DIMENSION> triangulation(Triangulation<DIMENSION>::limit_level_difference_at_vertices);
#endif    
    set_triangulation(triangulation, parameters);



    InitialConditionTitovDemoulin<EQUATIONS, DIMENSION> initial_condition(parameters, td_parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditionTDTest<DIMENSION> boundary_conditions(parameters, td_parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    AdaptivityTD<DIMENSION> adaptivity(parameters, mpi_communicator);
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
