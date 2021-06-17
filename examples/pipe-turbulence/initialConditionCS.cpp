#include "equationsMhd.h"
#include "completeEllipticIntegrals.h"
#include "initialConditionCS.h"
#include <random>
#include <cstdlib>
#include <math.h>

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

	double mydensity = 0.1;//kg/m3
	double mypressure = 0.02;//Pa
	double myB = 5.0-4;//T
	double mysize = 6.0e+6;//m
	double myspeed = 63078;//m/s
	double totalpressure = 1.;

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];


		value_list[i][0] = 0.1;

		value_list[i][1] = 0.;
		value_list[i][2] = 0;// 1. - p[0] * p[0];
		value_list[i][3] = 0.;

		



		value_list[i][5] = 0.;
		value_list[i][6] = p[0];
		value_list[i][7] = 0.;
	}

	//energy density

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		//if (i == 0) std::cout << this->parameters.gas_gamma;
		const Point<dim>& p = points[i];
		//mypressure =  0*std::cos(points[i][0] * M_PI / 300)* std::cos(points[i][0] * M_PI / 300);
		//value_list[i][4] = (1.15 + cs_parameters.beta - Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i])) / (2. / 3) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]);
		mypressure = totalpressure - value_list[i][6] * value_list[i][6];

		value_list[i][4] = mypressure / (this->parameters.gas_gamma - 1.0)+
			Equations<EquationsTypeMhd, dim >::compute_kinetic_energy(value_list[i]) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]);
	} 
}
template class InitialConditionCS<EquationsTypeMhd, 2>;
