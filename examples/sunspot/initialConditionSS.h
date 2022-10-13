#ifndef _INITIAL_CONDITION_SS_H
#define _INITIAL_CONDITION_SS_H

#include "util.h"
#include "parameters.h"
#include "equationsMhd.h"
#include "initialCondition.h"
#include "parametersSS.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

template <EquationsType equationsType, int dim>
class InitialConditionSS : public InitialCondition<equationsType, dim>
{
public:
    InitialConditionSS(Parameters<dim>&, SSParameters&);
    void vector_value(const std::vector<Point<dim> >&, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const;

private:
  SSParameters& cs_parameters;
  double invL_G;
  double d2R;
  double L2R;
  double R2L;
  double q_mag;

};

#endif
