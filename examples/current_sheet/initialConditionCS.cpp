#include "equationsMhd.h"
#include "completeEllipticIntegrals.h"
#include "initialConditionCS.h"
#include <random>
#include <cstdlib>

template <EquationsType equationsType, int dim>
InitialConditionCS<equationsType, dim>::InitialConditionCS(Parameters<dim>& parameters, CSParameters& cs_pars) :
      InitialCondition<equationsType, dim>(parameters), cs_parameters(cs_pars)
{
}

template <EquationsType equationsType, int dim>
void InitialConditionCS<equationsType, dim>::vector_value(const std::vector<Point<dim> >& points, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >& value_list) const
{
	Point<dim>& cb = this->parameters.corner_b;

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(-1, 1);


	/*random
	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
		value_list[i][0] = 1.; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;
		//momentum density perturbation //if (p[i] < 0.5)
		//value_list[i][1] = 0.005 * std::sin(std::rand());
		//value_list[i][2] = 0.005 * std::sin(std::rand());
		if (((abs(points[i][0]) < (cb[0] - 1.)) and (abs(points[i][1]) < (cb[1] - 1.))))
		{
			value_list[i][1] = 0.005 * dis(gen);
			value_list[i][2] = 0.005 * dis(gen);
		}
		value_list[i][3] = 0.;

			


		value_list[i][5] = 0.0; //magnetic field
		value_list[i][6] = std::tanh(p[0]);
		value_list[i][7] = 0.;

	}
	*/
	
	//start outflow
	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
		value_list[i][0] = 1.; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;

		if (((abs(points[i][0]) < (cb[0] - 0.5)) and (abs(points[i][1]) < (cb[1] - 1.))))
		{
			value_list[i][1] = 0.;
			value_list[i][2] = 0.08 * (0 < points[i][1]) - (points[i][1] < 0);
		}
		value_list[i][3] = 0.;




		value_list[i][5] = 0.0; //magnetic field
		value_list[i][6] = std::tanh(p[0]);
		value_list[i][7] = 0.;

	}

	//energy density

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
		value_list[i][4] = (1.15 + cs_parameters.beta - Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i])) / (2. / 3) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]);
	}
}
template class InitialConditionCS<EquationsTypeMhd, 3>;
