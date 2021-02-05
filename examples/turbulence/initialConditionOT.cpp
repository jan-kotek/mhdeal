#include "completeEllipticIntegrals.h"
#include "initialConditionOT.h"
#include "equationsMhd.h"
#include <random>
#include <cstdlib>

template <EquationsType equationsType, int dim>
InitialConditionOT<equationsType, dim>::InitialConditionOT(Parameters<dim>& parameters) :
  InitialCondition<equationsType, dim>(parameters)
{
}


template <EquationsType equationsType, int dim>
void InitialConditionOT<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&result) const
{
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(-1, 1);

  
  int pressure = 5. / (12. * 3.14);
  
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;

    result[i][1] = 0.;
    result[i][2] = 2.;
    result[i][3] = 0.;

	if ((abs(points[i][0]) < 4.) and (abs(points[i][1]) < 4.))
	{
		result[i][1] += 0.001 * dis(gen);
		result[i][2] += 0.001 * dis(gen);
	}

    result[i][5] = 0;
    result[i][6] = 1.;
    result[i][7] = 0.0;
    result[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (result[i][5] * result[i][5] + result[i][6] * result[i][6] + result[i][7] * result[i][7]) +
      0.5 * (result[i][1] * result[i][1] + result[i][2] * result[i][2] + result[i][3] * result[i][3]) / result[i][0];
	//result[i][4] = (1.15 + beta - Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(result[i])) / (2. / 3) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(result[i]);

  }
}

template class InitialConditionOT<EquationsTypeMhd, 3>;
