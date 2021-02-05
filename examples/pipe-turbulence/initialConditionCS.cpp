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

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
/*		value_list[i][0] = mydensity;//1.+std::sin(points[i][0]* M_PI / 600.); //density

		value_list[i][1] = 2000*std::cos(2.68*points[i][0]* mysize*1.0e-7)+45.249*std::cosh(1.93*points[i][0] * mysize *1.0e-7);//1.*std::cosh(points[i][0] * M_PI / 600.); //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;// 1.*std::sinh(points[i][0] * M_PI / 600.);

		//momentum density perturbation //if (p[i] < 0.5)

		//value_list[i][1] = 0.005 * std::sin(std::rand());
		//value_list[i][2] = 0.005 * std::sin(std::rand());
		  //if (((abs(points[i][0]) < (cb[0] - 1.)) and (abs(points[i][1]) < (cb[1] - 1.))))
		//if (abs(points[i][0])<4.)
		 //{
		 //	value_list[i][1] += 0.0000001 * dis(gen);
		 //	value_list[i][2] += 0.0000001 * dis(gen);
		 //}
	//	value_list[i][3] = 0.;
    
			


		value_list[i][5] =myB;// 1.*std::cosh(points[i][0] * M_PI / 600.); //magnetic field
		value_list[i][6] =0.;//std::tanh(p[0]);
		value_list[i][7] = -7.4 * std::cos(2.68 * points[i][0] * mysize * 1.0e-7) + 3.703 * std::cosh(1.93 * points[i][0] * mysize * 1.0e-7);//1*std::cosh(points[i][0] * M_PI / 600.);
*/

		value_list[i][0] = 1.;

		value_list[i][1] = 0.;
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;

		



		value_list[i][5] = 0.;
		value_list[i][6] = 0.;
		value_list[i][7] = 0.;
	}

	//energy density

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		//if (i == 0) std::cout << this->parameters.gas_gamma;
		const Point<dim>& p = points[i];
		//mypressure =  0*std::cos(points[i][0] * M_PI / 300)* std::cos(points[i][0] * M_PI / 300);
		//value_list[i][4] = (1.15 + cs_parameters.beta - Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i])) / (2. / 3) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]);
		value_list[i][4] = mypressure / (this->parameters.gas_gamma - 1.0)+
			Equations<EquationsTypeMhd, dim >::compute_kinetic_energy(value_list[i]) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]);
	} 
}
template class InitialConditionCS<EquationsTypeMhd, 3>;
