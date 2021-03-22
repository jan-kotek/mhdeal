#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "initialConditionMhdBlast.h"
#include "adaptivityMhdBlast.h"
#include "parameters.h"

// Dimension of the problem - passed as a template parameter to pretty much every class.
#define DIMENSION 2
// Type of equations, must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd

#ifdef HAVE_MPI
void set_triangulation(parallel::distributed::Triangulation<DIMENSION>& triangulation,  Parameters<DIMENSION>& parameters)
#else
void set_triangulation(Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#endif
{
  //GridGenerator::subdivided_hyper_rectangle(triangulation, parameters.refinements, parameters.corner_a, parameters.corner_b, true);
    GridGenerator::hyper_ball(triangulation);
    static const SphericalManifold<DIMENSION> boundary;
    triangulation.set_all_manifold_ids_on_boundary(0);
    triangulation.set_manifold(0, boundary);
    const Point<DIMENSION> mesh_center;
    const double core_radius = 1.0 / 5.0,
        inner_radius = 1.0 / 3.0;


    // Step 1: Shrink the inner cell
    //
    // We cannot get a circle out of the inner cell because of
    // the degeneration problem mentioned above. Rather, shrink
    // the inner cell to a core radius of 1/5 that stays
    // sufficiently far away from the place where the
    // coefficient will have a discontinuity and where we want
    // to have cell interfaces that actually lie on a circle.
    // We do this shrinking by just scaling the location of each
    // of the vertices, given that the center of the circle is
    // simply the origin of the coordinate system.
    for (typename Triangulation<DIMENSION>::active_cell_iterator cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
        if (mesh_center.distance(cell->center()) < 1e-5)
        {
            for (unsigned int v = 0;
                v < GeometryInfo<DIMENSION>::vertices_per_cell;
                ++v)
                cell->vertex(v) *= core_radius / mesh_center.distance(cell->vertex(v));
        }
    // Step 2: Refine all cells except the central one
    for (typename Triangulation<DIMENSION>::active_cell_iterator cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
        if (mesh_center.distance(cell->center()) >= 1e-5)
            cell->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
    // Step 3: Resize the inner children of the outer cells
    //
    // The previous step replaced each of the four outer cells
    // by its four children, but the radial distance at which we
    // have intersected is not what we want to later refinement
    // steps. Consequently, move the vertices that were just
    // created in radial direction to a place where we need
    // them.
    for (typename Triangulation<DIMENSION>::active_cell_iterator cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
        for (unsigned int v = 0; v < GeometryInfo<DIMENSION>::vertices_per_cell; ++v)
        {
            const double dist = mesh_center.distance(cell->vertex(v));
            if (dist > core_radius * 1.0001 && dist < 0.9999)
                cell->vertex(v) *= inner_radius / dist;
        }
    // Step 4: Apply curved manifold description
    //
    // As discussed above, we can not expect to subdivide the
    // inner four cells (or their faces) onto concentric rings,
    // but we can do so for all other cells that are located
    // outside the inner radius. To this end, we loop over all
    // cells and determine whether it is in this zone. If it
    // isn't, then we set the manifold description of the cell
    // and all of its bounding faces to the one that describes
    // the spherical manifold already introduced above and that
    // will be used for all further mesh refinement.
    for (typename Triangulation<DIMENSION>::active_cell_iterator cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
    {
        bool is_in_inner_circle = false;
        for (unsigned int v = 0; v < GeometryInfo<DIMENSION>::vertices_per_cell; ++v)
            if (mesh_center.distance(cell->vertex(v)) < inner_radius)
            {
                is_in_inner_circle = true;
                break;
            }
        if (is_in_inner_circle == false)
            cell->set_all_manifold_ids(0);
    }

    //triangulation.refine_global(4);

 /*std::vector<dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<DIMENSION> > > > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
  */
}

void set_parameters(Parameters<DIMENSION>& parameters)
{
    parameters.corner_a = Point<DIMENSION>(-0.5, -0.75);
    parameters.corner_b = Point<DIMENSION>(0.5, 0.75);
    //parameters.refinements = { 25, 32 };
    parameters.limit = true;
    parameters.limitB = true;
    parameters.limit_edges_and_vertices = true;
    parameters.output_file_prefix = "RD-mhd-blast";
    parameters.slope_limiter = parameters.vertexBased;
    parameters.use_div_free_space_for_B = false;
    //parameters.periodic_boundaries = {{ 0, 1, 0 }};
    parameters.num_flux_type = Parameters<DIMENSION>::lax_friedrich;
    parameters.lax_friedrich_stabilization_value = 0.5;
    parameters.cfl_coefficient = .05;
    parameters.quadrature_order = 5;
    parameters.polynomial_order_dg = 1;
    parameters.patches = 0;
    parameters.output_step = 1.e-3;
    parameters.final_time = .5;
    //parameters.debug = parameters.Adaptivity; // | parameters.BasicSteps | parameters.PeriodicBoundaries;// | parameters.Assembling;
    //
    /*
    parameters.output_matrix = true;
    parameters.output_rhs = true;
    parameters.output_solution = true;
    */

    parameters.max_cells = 900;
    parameters.refine_every_nth_time_step = 5;
    parameters.perform_n_initial_refinements = 10;
    parameters.refine_threshold = 0.5;
    parameters.coarsen_threshold = 0.2;
    parameters.volume_factor = 2;
    parameters.time_interval_max_cells_multiplicator = 1.;
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

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
  
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::limit_level_difference_at_vertices));
    
#else
    Triangulation<DIMENSION> triangulation(Triangulation<DIMENSION>::MeshSmoothing::limit_level_difference_at_vertices);
    
#endif
   
    set_triangulation(triangulation, parameters);

    InitialConditionMhdBlast<EQUATIONS, DIMENSION> initial_condition(parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryCondition<EQUATIONS, DIMENSION> boundary_conditions(parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    AdaptivityMhdBlast<DIMENSION> adaptivity(parameters, mpi_communicator);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Set adaptivity
    problem.set_adaptivity(&adaptivity);
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
