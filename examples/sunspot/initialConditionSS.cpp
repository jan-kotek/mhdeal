#include "equationsMhd.h"
#include "completeEllipticIntegrals.h"
#include "initialConditionSS.h"
#include <random>
#include <cstdlib>

extern std::vector<std::vector<std::string>> inputvector;
template <EquationsType equationsType, int dim>
InitialConditionSS<equationsType, dim>::InitialConditionSS(Parameters<dim>& parameters, SSParameters& cs_pars) :
      InitialCondition<equationsType, dim>(parameters), cs_parameters(cs_pars)
{
}

template <EquationsType equationsType, int dim>
void InitialConditionSS<equationsType, dim>::vector_value(const std::vector<Point<dim> >& points, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >& value_list) const
{
	Point<dim>& cb = this->parameters.corner_b;

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(-1, 1);



	//open input binary file and fill the grid.
	//static background
	
	/*
	double point_array[5];
	double pressure;
	int size = sizeof(point_array);
	
	std::ifstream skvrna("mhs_data.out", std::ios::in | std::ios::binary);
	skvrna.seekg(0, std::ios::beg);
	std::cout << "points size " << points.size() << "  ";
	//for (unsigned int i = 0; i < points.size(); ++i) //input from external file
	{   
		
		skvrna.read(reinterpret_cast<char*>(&point_array), size);
		const Point<dim>& p = points[i];
		//std::cout << " " << p << " "<<point_array;
		value_list[i][0] = point_array[0]; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;
		

		
		pressure = point_array[1];
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (B_loc[0] * B_loc[0] + B_loc[1] * B_loc[1] + B_loc[2] * B_loc[2]); //energy density
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]); //energy density

	


		value_list[i][5] = point_array[2]; //magnetic field
		value_list[i][6] = point_array[3];
		value_list[i][7] = point_array[4];

		value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (point_array[2] * point_array[2] + point_array[3] * point_array[3] + point_array[4] * point_array[4]); //energy density

		
		
	}
	for (unsigned int i = 0; i < points.size(); ++i)
	{

		skvrna.read(reinterpret_cast<char*>(&point_array), size);
		const Point<dim>& p = points[i];
		//std::cout << " " << p << " "<<point_array;
		value_list[i][0] = point_array[0]; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;



		pressure = point_array[1];
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (B_loc[0] * B_loc[0] + B_loc[1] * B_loc[1] + B_loc[2] * B_loc[2]); //energy density
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]); //energy density




		value_list[i][5] = point_array[2]; //magnetic field
		value_list[i][6] = point_array[3];
		value_list[i][7] = point_array[4];

		value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (point_array[2] * point_array[2] + point_array[3] * point_array[3] + point_array[4] * point_array[4]); //energy density



	}
	std::cout << "input done \n";
	skvrna.close();
	*/

	std::cout << inputvector[1][1] << " ";
	double r; //radial coordinate
	double pressure;
	double sdx = 1./0.01;//input file dx
	double sdy = 1./0.01;
	double x;
	double y;
	
	/*std::cout << "open file ";
	std::ifstream skvrna("mhs_data.csv", std::ios::in);
	int rows = 0;
	std::string line;
	std::cout << "while points size " << points.size();
	while (getline(skvrna, line))
	rows++;
	std::cout << " points size " << points.size() << " rows  "<< rows<<" \n";
	*/
	for (unsigned int i = 0; i < points.size(); ++i)
	{

		
		const Point<dim>& p = points[i];
		r = (p[0] * p[0] + p[1] * p[1]);
		x=r;
		y=p[2];
		x *= sdx;
		y *= sdy;

		double xLow = floor(x);
		double yLow = floor(y);

		double xCell = x - xLow;
		double yCell = y - yLow;

		int ix = (int)xLow;
		int iy = (int)yLow;
		/*
		for (int yy = -1; yy < 3; yy++) // loop over adjacent lines
		{
			for (int xx = -1; xx < 3; xx++) // fill the coarse-valued array of neighbours
				cxVals[xx + 1] = current(ix + xx, iy + yy);
			inputvector[1][6]
			cyVals[yy + 1] = hcInterp(cxVals, xCell); // interpolate one line
		}
		*/
		//std::cout << " " << p << " "<<point_array;

		3Dinputvektor[ix][iy][6];=current(ix,iy)

		value_list[i][0] = 1.; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;



		pressure = 1.;
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (B_loc[0] * B_loc[0] + B_loc[1] * B_loc[1] + B_loc[2] * B_loc[2]); //energy density
		//value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]); //energy density




		value_list[i][5] = 1.; //magnetic field
		value_list[i][6] = 1.;
		value_list[i][7] = 1.;

		value_list[i][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (value_list[i][5] * value_list[i][5] + value_list[i][6] * value_list[i][6] + value_list[i][7] * value_list[i][7]); //energy density



	}
	std::cout << "input done \n";
	skvrna.close();
}
template class InitialConditionSS<EquationsTypeMhd, 3>;
