// Must be _before_ #include "field.h" 

//sompile: g++ -I/home/jfk/programy/MHD3d-recent/include -lgsl -lgslcblas -lm -O1 -o mhs_vtk main.cc
//stand alone compile: g++ main.cc -I/home/jfk/prg/MHD3d-recent/include -lgslcblas -lm -lgsl -O1 -o mhs_vtkls


#define real_t   double //float




#include "field.h" //mhd3D c++ library mailto:jankotek@asu.cas.cz
//#include "ffield.h"f

#include "gsl/gsl_sf_ellint.h"


#include <fstream>
#include <algorithm>
#include <cstdio>
#include <ctime>
using namespace std;


//============================================================================
//--- Geometry
bool output3D=0;
bool output_binary=0;
bool output_csv=0;
bool output_csv_2D=1;

int boxHsize=100; // make it even in order to avoid on-axis singularities
int boxVsize=200;

real_t hCellSize=0.01;
real_t vCellSize=0.01;//0.02


//--- B-field params

// reference point
real_t r0=0.1;//1
real_t z0=0.1;//0

// B_phi at the reference point
real_t B_phi_ref=0.;//0.1

// B_phi (torsion) parameter 'k' in Eq. (2.1)
real_t Kt=0.;//1.0 <- 

// free parameter for A_phi in Eq. (5.1)
real_t Ka=0.00000;//100        


//-- Current-loop generated B-field (photo-/chromosphere)

// ring current (5.2)
real_t I=3.;//10,3
// normalised to 2PI (just a shortcut)
real_t I_norm=I/(2.0*M_PI);

// current-loop radius
real_t r_I=0.5;
// current-loop z-coordinate (negative means submerging under photosphere)
real_t z_I=-0.2; // avoid singularities around the 'wire'

//-- Vertical homogeneous B-field (corona)

// vertical field Bz for z -> infty 
real_t B0=2.;//2


     

//--- Temperature/pressure controlling params

// gravity height-scale
real_t L_g=0.5;//3,0.2
// specific heat ratio gamma=c_p/c_V (=5/3 for atomic gas), kappa=gamma-1
real_t cpcv=1.67;
real_t kappa=cpcv-1.0;
// on-axis typical plasma beta
real_t beta=20.;//0.2 
// on-axis typical temperature
real_t w_0=1.05;//0.2            


// on-axis pressure at the base (derived)
//real_t p_0=beta*B0*B0;//*100
real_t p_0=beta*20.;//*100



/*
//--- Geometry

int boxHsize=100; // make it even in order to avoid on-axis singularities
int boxVsize=200;

real_t hCellSize=0.01;
real_t vCellSize=0.01;


//--- B-field params

// reference point
real_t r0=1.0;
real_t z0=0.1;

// B_phi at the reference point
real_t B_phi_ref=0.5;

// B_phi (torsion) parameter 'k' in Eq. (2.1)
real_t Kt=5.0;


//-- Current-loop generated B-field (photo-/chromosphere)

// ring current (5.2)
real_t I=10.0;
// normalised to 2PI (just a shortcut)
real_t I_norm=I/(2.0*M_PI);

// current-loop radius
real_t r_I=0.5;//0.2
// current-loop z-coordinate (negative means submerging under photosphere)
real_t z_I=-0.02; //0.01 avoid singularities around the 'wire'

//-- Vertical homogeneous B-field (corona)

// vertical field Bz for z -> infty 
real_t B0=0.1;//1.

// free parameter for A_phi in Eq. (5.1)
real_t Ka=0.;//0.

//--- Temperature/pressure controlling params

// gravity height-scale
real_t L_g=0.05;//10
// specific heat ratio gamma=c_p/c_V (=5/3 for atomic gas), kappa=gamma-1
real_t cpcv=1.67;
real_t kappa=cpcv-1.0;
// on-axis typical plasma beta
real_t beta=0.2;
// on-axis typical temperature
real_t w_0=0.2;//0.2


// on-axis pressure at the base (derived)
real_t p_0=beta*B0*B0;
*/


//============================================================================


/// Big/little endian encoding reverse function 



template <class T2>
T2 ReverseFloat(T2 in)
{
    char* const p = reinterpret_cast<char*>(&in);
    for (size_t i = 0; i < sizeof(T2) / 2; ++i)
        std::swap(p[i], p[sizeof(T2) - i - 1]);
    return in;
}

/*float ReverseFloat2( const float inFloat )
{  
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}
*/



template <typename T>
T bswap(T val) {
    T retVal;
    char *pVal = (char*) &val;
    char *pRetVal = (char*)&retVal;
    int size = sizeof(T);
    for(int i=0; i<size; i++) {
        pRetVal[size-1-i] = pVal[i];
    }

    return retVal;
}
//============================================================================

int main()
{

  //--- 3D box for display

  TGeometry* box=new TGeometry(boxHsize,boxHsize,boxVsize);
  box->setMeshSizes(hCellSize,hCellSize,vCellSize);


   
  // Place the axis in the center of the box
  real_t x0=0.5*box->dx()*(box->xSize()-1);
  real_t y0=0.5*box->dy()*(box->ySize()-1);


  //--- 2D radial-vertical (r,z,phi=0) half-slice for fields
  //    Make it as the largest (diagonal) slice to the box

  // dr: alias for dx in the 2D slice. NB: dx must be the same as dy for 3D box 
  real_t dr=box->dx();
 
  size_t diagSize=(size_t) (ceil(sqrt(x0*x0+y0*y0)/dr));
  std::cout<<"diagSize "<<diagSize<<"\n";



  // basic field geometry
  TGeometry* fg=new TGeometry(diagSize,box->zSize());
  fg->setMeshSizes(dr,box->dz());

  // extended geometries for B and A - need boundaries for num. derivative calc.
  // 2 boundary points in r for j~A" and even 3 in z as dL_r/dz~A'''
  TGeometry* ag=new TGeometry(*fg);
  ag->setBoundarySizes(2,3);

  // bg: boundaries of 1 less than ag
  TGeometry* bg=new TGeometry(*fg);
  bg->setBoundarySizes(1,2);

  // fg: boundaries of 1 less then bg - just 1 in z-direction for dL_r/dz calc.
  fg->setBoundarySizes(0,1);


  //--- Fill the 2D cylindrical radial slice phi=0


  // A equiv. A_phi - to be selected
  RScalarField A(ag);


  // B components (calculated from A)
  RScalarField B_r(bg);
  RScalarField B_phi(bg);
  RScalarField B_z(bg);


  // j components: calculated as j = rot B
  RScalarField j_r(fg);
  RScalarField j_phi(fg);
  RScalarField j_z(fg);

  // Lorentz force components: calculated as L = j x B
  RScalarField L_r(fg);
  RScalarField L_z(fg);


  // presssure, density and temperature distributions: Calculated from MHS 
  RScalarField rho(fg);
  RScalarField tau(fg);
  RScalarField w(fg);


  RScalarField p(fg);
  RScalarField RBeta(fg);

  //--------- Init free-selectable quantities

  // Init the vector potential

  for(int iz=ag->yBegin();iz<ag->yEnd();iz++)
    {
      real_t z=iz*ag->dy();

      for(int ir=ag->xBegin();ir<ag->xEnd();ir++)
	{
	  // half-cell off axis, symmetric
    real_t signed_r=ir*dr+0.5*dr; //a je treba vytvorit signum r kvuli vektorove operaci.
	  real_t r=fabs(signed_r);//tohle by mel bejt unsigned real_t

     
	  
	  // Argument of elliptical integral (eq. 5.2a)
	  real_t kr=2.0*sqrt(r*r_I/((r+r_I)*(r+r_I)+(z-z_I)*(z-z_I)));
    //   if (iz<0) kr*=-1.;
	  real_t Ek=gsl_sf_ellint_Ecomp(kr,GSL_PREC_DOUBLE);
	  real_t Kk=gsl_sf_ellint_Kcomp(kr,GSL_PREC_DOUBLE);

	  // I-generated mg. field (Eq. 5.2b)
	  real_t Ak=((2.0-kr*kr)*Kk-2.0*Ek)/kr*I_norm*sqrt(r_I/r);
         if (ir<0) Ak*=-1.;
	  // Total A_phi: I-generated + that corresponding to Bz=B0
	  A(ir,iz)=Ak+0.5*B0*signed_r+Ka/signed_r;//zmena
    /*
    if ((iz<2)&&(ir<2)) std::cout<<"ir "<<ir<<" iz "<<iz<<" Ak "<<Ak
       <<" kr "<<kr<<" Ek "<<Ek<<" Kk "<<Kk<<" A(r) "<<A(ir,iz)<<" r_s "<<signed_r<<"\n";
	  */
  }
    }

  // Init the vertical temperature profile at the axis
  TGeometry tp(box->zSize());
  tp.setBoundarySizes(1);
  RScalarField w0(&tp);


  ifstream w0data ("w0i.dat");
  if (w0data.is_open())
  { 
    std::cout<<"\n"<<"C7"<<"\n"; 
    std::string line;
    real_t line_value; 
    int iz=0;
    while ( getline (w0data,line))
    {
     line_value=std::stod(line); 
     if ((iz%2)==0) w0(iz/2)=line_value;
     //std::cout<<iz<<"iz "<<line<<"line "<<line_value<<"\n";
     iz++;
     //std::cout<<"iz"<<iz;
    }
    w0data.close();

  }
  else for(int iz=1;iz<box->zSize();iz++)
  {  
    std::cout<<"\n"<<"analytical_w0"<<"\n"; 
    w0(0)=w_0;

   // w0(iz)=w_0; // isothermal height profile
   //if   (iz<=100)   w0(iz)=w0(iz-1)*(0.99);  else if (iz>100) w0(iz)=w0(iz-1)*(1.02); //v1
   if   (iz<=100)   w0(iz)=w0(0);  else if (iz<=120) w0(iz)=w0(iz-1)+w0(0)*6; else w0(iz)=w0(iz-1);                                                                                   //corona_v2
   //else if (iz<=200) w0(iz)=w0(iz)*(0.5+((iz-100)*0.005));
   
  }//--------- Calculate derived quantities


  //--- On-axis plasma thermo


  // On-axis inverse-height profile tau0(z) (eq. '*')
  RScalarField tau0(&tp);

  // On-axis pressure profile p0(z) (eq. 3a0)
  RScalarField p0(&tp);

  // at the base
  tau0(0)=1.0/L_g;
  p0(0)=p_0;

  // integrated tau0
  real_t tau0_int=0.0;

  for(int iz=1;iz<box->zSize();iz++)
    { 
      tau0(iz)=w_0/(L_g*w0(iz));

      // integrate tau0(z) by trapezoidal rule
      tau0_int+=0.5*(tau0(iz-1)+tau0(iz))*box->dz();

      p0(iz)=p_0*exp(-tau0_int);
      if (p0(iz)<0.) std::cout<<"ERR p0<0: iz,tau,p0: "<<iz<<" , "<<tau0(iz)<<" , "<<p0(iz)<<"\n";
    }


  //--- Magnetic field components


  real_t invDz=1.0/box->dz();
  real_t invDr=1.0/dr;

  // flux function in the reference point
  real_t A_ref=A.atc(r0,z0);


  for(int iz=bg->yBegin();iz<bg->yEnd();iz++)
    for(int ir=bg->xBegin();ir<bg->xEnd();ir++)
      {
  real_t signed_r=ir*dr+0.5*dr;
	real_t r=fabs(signed_r);
  int    signum_r;
  if (signed_r==r) signum_r=1; 
  else signum_r=-1;

	B_r(ir,iz)=-0.5*invDz*(A(ir,iz+1)-A(ir,iz-1));
	B_phi(ir,iz)=B_phi_ref*pow(((A(ir,iz)*signed_r)/(A_ref*r0)),Kt)*r0/r;//orig
  //B_phi(ir,iz)=B_phi_ref*pow(((A(ir,iz)*signed_r)/(A_ref*r0)),Kt)*r0/(1+r); //better nit unstable
  // B_phi(ir,iz)=B_phi_ref*signum_r*powf(((A(ir,iz)*r)/(A_ref*r0)),Kt)*r0/r;//signed?
  
  B_z(ir,iz)=0.5*invDr*(A(ir+1,iz)-A(ir-1,iz))+A(ir,iz)/signed_r;

  if ((iz==100)&&(ir<4)) 
  { std::cout<<"ir "<<ir<<"dr "<< dr <<" invDr "<< invDr<<" iz "<<iz<<" B_z "<<B_z(ir,iz)<<" B_phi "<<B_phi(ir,iz)<<" B_r "<<B_r(ir,iz)
      <<" A(r+1) "<<A(ir+1,iz)<<" A(r) "<<A(ir,iz)<<" A(r-1) "<<A(ir-1,iz)<<" signed_r "<<signed_r<<" idr "<<invDr<<"\n";
    
  }
      }

  


/*   ofstream deb;
  deb.open("2d.data",ios::trunc);
  deb<<"z;x;Br;Bphi;Bz;"<<"\n";
  for(int i=0;i<boxVsize;i++)
  {
      for(int j=0;j<boxHsize;j++) 
    {
      deb<<i*vCellSize<<";"<<j*hCellSize<<";"<<B_r(i,j)<<";"<<B_phi(i,j)<<";"<<B_z(i,j)<<"\n";
      std::cout<<i*vCellSize<<";"<<j*hCellSize<<";"<<B_r(i,j)<<";"<<B_phi(i,j)<<";"<<B_z(i,j)<<"\n";
    }
    
  }
 deb.close(); */

  //--- Current density

  for(int iz=fg->yBegin();iz<fg->yEnd();iz++)
    for(int ir=fg->xBegin();ir<fg->xEnd();ir++)
      {
	 real_t signed_r=ir*dr+0.5*dr;
	 real_t r=fabs(signed_r);

//central differecne
	j_r(ir,iz)=-0.5*invDz*(B_phi(ir,iz+1)-B_phi(ir,iz-1));
	j_phi(ir,iz)=0.5*invDz*(B_r(ir,iz+1)-B_r(ir,iz-1))-0.5*invDr*(B_z(ir+1,iz)-B_z(ir-1,iz));
  
	j_z(ir,iz)=0.5*invDr*(B_phi(ir+1,iz)-B_phi(ir-1,iz))+B_phi(ir,iz)/r; //zmena1
 
  //j_z(ir,iz)=0.;
  //if (ir==0) j_z(ir,iz)=-1.;
//  j_z(ir,iz)=invDr*(B_phi(ir,iz)-B_phi(ir-1,iz))+B_phi(ir,iz)/r; //backward
  
//forward difference
  /*(ir,iz)=-invDz*(B_phi(ir,iz+1)-B_phi(ir,iz));
	j_phi(ir,iz)=invDz*(B_r(ir,iz+1)-B_r(ir,iz))-invDr*(B_z(ir+1,iz)-B_z(ir,iz));
  j_z(ir,iz)=invDr*(B_phi(ir+1,iz)-B_phi(ir,iz))+B_phi(ir,iz)/signed_r; 
  */
  //if ((iz==100)&&(ir<5)) std::cout<<"ir "<<ir<<" signed_r "<<signed_r<<" j_z "<<j_z(ir,iz)<<" j_phi "<<j_phi(ir,iz)<<"\n";
  //  if ((iz<2)&&(ir<2)) std::cout<<"ir "<<ir<<" iz "<<iz<<" j_phi "<<j_phi(ir,iz)<<"\n";
     }
  
  //--- Lorentz force

  for(int iz=fg->yBegin();iz<fg->yEnd();iz++)
    for(int ir=fg->xBegin();ir<fg->xEnd();ir++)
      {
	//L_r(ir,iz)=j_phi(ir,iz)*B_z(ir,iz)-j_z(ir,iz)*B_phi(ir,iz);
  L_r(ir,iz)=j_phi(ir,iz)*B_z(ir,iz)-j_z(ir,iz)*B_phi(ir,iz);
   
	L_z(ir,iz)=j_r(ir,iz)*B_phi(ir,iz)-j_phi(ir,iz)*B_r(ir,iz);
      }




  //=== Plasma thermo: No boundaries, just internal region

  //--- Pressure 

  for(int iz=0;iz<box->zSize();iz++)
    {
      // Integral_0^r L_r(r',z) dr' - by trapezoidal rule
      // The ir=0 point is 0.5dr off-axis: Calculate the contribution
      // of this half-cell, taking into account that L at axis is zero.
      // (one 1/2 factor is for average between zero and L_r at (0,iz) and 
      // the 2nd one for half-cell size).
      real_t L_r_int=0.25*dr*L_r(0,iz);

      // Eq. (4.1) for r=0.5dr
      p(0,iz)=p0(iz)+2.0*L_r_int;
      //if (p(0,iz)<0.) std::cout<<"ERR p(0,iz)<0: iz "<<iz<<" p(0,iz) "<<p(0,iz)<<" L_r_int "<<L_r_int
      //                       <<" dr "<<"j_phi "<<j_phi(0,iz)<<" Bz "<<B_z(0,iz)<<"\n";
  // std::cout<<"j_phi(0,iz): "<<j_phi(0,iz)<<" B_z(0,iz): "<<B_z(0,iz)<<" j_z(0,iz):"<<j_z(0,iz)<<" B_phi(0,iz):"<<B_phi(0,iz);
  //  std::cout<<" L_r(0,iz): "<<L_r(0,iz)<<" p0: "<<p(0,iz)<<" L_r_int:"<<L_r_int<<"\n";
      for(int ir=1;ir<fg->xSize();ir++)
	{
	  // trapezoidal integration
    L_r_int+=0.5*(L_r(ir-1,iz)+L_r(ir,iz))*dr; //back
	  //L_r_int+=0.25*(L_r(ir-1,iz)+L_r(ir+1,iz))*dr;//center
    //L_r_int+=0.5*(L_r(ir,iz)+L_r(ir+1,iz))*dr;//forward
	 
	  //Eq. 4.1
	  p(ir,iz)=2.0*L_r_int+p0(iz);
   /* if (p(ir,iz)<0.) 
    {std::cout<<"ERR p(ir,iz)<0: ir "<<ir<<" iz "<<iz<<" p(ir,iz) "<<p(ir,iz)<<" L_r_int "<<L_r_int
                           <<" p0(iz) "<<p0(iz)<<"\n";
     //p(ir,iz)=1e-9;                        
    }   */                     
    RBeta(ir,iz)= p(ir,iz)/(B_z(ir,iz)*B_z(ir,iz)+B_r(ir,iz)*B_r(ir,iz)+B_phi(ir,iz)*B_phi(ir,iz));
    
	}
    }


  //--- Tau (inv. local height-scale)

 for(int iz=0;iz<box->zSize();iz++)
    {
      // Integral_0^r dL_r(r',z)/dz dr' - by trapezoidal rule
      // The ir=0 point is 0.5dr off-axis - see the note above
      // The on emore 1/2 is from centered differenciation of d/dz

      //real_t dLrdz_int=0.125*dr*invDz*(L_r(0,iz+1)-L_r(0,iz-1)); //central
      real_t dLrdz_int=0.25*dr*invDz*(L_r(0,iz+1)-L_r(0,iz)); //forward

      // Eq. (4.2) for r=0.5dr
      tau(0,iz)=(L_z(0,iz)-dLrdz_int+0.5*p0(iz)*tau0(iz))/(0.5*p(0,iz));

      for(int ir=1;ir<fg->xSize();ir++)
	{
	  // trapezoidal integration of dL_r/dz
	  dLrdz_int+=0.25*dr*invDz*(L_r(ir-1,iz+1)-L_r(ir-1,iz-1)+L_r(ir,iz+1)-L_r(ir,iz-1));//central?
    //dLrdz_int+=0.5*dr*invDz*(L_r(ir,iz+1)-L_r(ir,iz)+L_r(ir,iz+1)-L_r(ir,iz));//forward?


	  // Eq. (4.2)
	  tau(ir,iz)=(L_z(ir,iz)-dLrdz_int+0.5*p0(iz)*tau0(iz))/(0.5*p(ir,iz));
    
   /* if (tau(ir,iz)<0.) 
      std::cout<<" ERR: tau<0 ir,iz:"<< ir <<" / "<< iz << " Tau: "<<tau(ir,iz)<<" Lz "<<L_z(ir,iz)<<" dLrdz_int:"<<dLrdz_int<<" RHS:"<<0.5*p0(iz)*tau0(iz)
              << "w_0: "<<w_0<<" L_g: "<<L_g<<" p:"<<p(ir,iz)<<" p0:"<<p0(iz)<<"\n"; 
    */         
	}
    }  

 //--- Density,temperature
 std::cout<<" Errors: \n";
 for(int iz=0;iz<box->zSize();iz++)
   {for(int ir=0;ir<fg->xSize();ir++)
     {
       rho(ir,iz)=L_g/(kappa*w_0)*tau(ir,iz)*p(ir,iz);
       w(ir,iz)=w_0/(L_g*tau(ir,iz));
       //if ((iz==5)&&(ir==40)) std::cout<<"w_0: "<<w_0<<" L_g: "<<L_g<<" Tau: "<<tau(ir,iz)<<" p:"<<p(ir,iz)<<
       //                                   "Lz:"<<L_z(ir,iz)<<"\n";
        
      /* if (((p(ir,iz)<0.)&&(iz<180))||(isnan(p(ir,iz)))) std::cout<<" ir "<<ir<<" iz "<<iz<<" p(ir,iz) "<<p(ir,iz) <<"\n";
         else if (((rho(ir,iz)<0.)&&(iz<180))||(isnan(rho(ir,iz)))) std::cout<<" ir "<<ir<<" iz "<<iz<<" rho(ir,iz) "<<rho(ir,iz) <<"\n";
         else if (((w(ir,iz)<0.)&&(iz<180))||(isnan(w(ir,iz)))) std::cout<<" ir "<<ir<<" iz "<<iz<<" w(ir,iz) "<<w(ir,iz) <<"\n";
      */
     }
   }
/*
  ofstream dd;
  
  dd.open("2d.csv");
 
  dd<<"z;r;Br;Bphi;Bz;Jr;Jphi;Jz;Lr;Lz;p;tau;rho;w\n";
  
 for(int iz=bg->yBegin();iz<bg->yEnd();iz++)
 {   for(int ir=bg->xBegin();ir<bg->xEnd();ir++)
    { std::cout<<"iz"<<iz<<" ir"<<ir<<" BE:"<< bg->xEnd()<<" - "<<bg->yEnd()<<"\n";
      std::cout<<"&&";    
      std::cout <<  iz<<";"<<ir<<";"<<B_r(ir,iz)<<";"<<B_phi(ir,iz)<<";"<<B_z(ir,iz)<<";"  //vCellSize
        <<j_r(ir,iz)<<";"<<j_phi(ir,iz)<<";"<<j_z(ir,iz)<<";"<<L_r(ir,iz)<<";"<<L_z(ir,iz)<<";"
        <<p(ir,iz)<<";"<<tau(ir,iz)<<";"<<rho(ir,iz)<<";"<<w(ir,iz)
         <<"\n";
      
    }
    std::cout<<"@@"; 
 }     
   
  
  std::cout<<"closed";
 dd.close();
  std::cout<<"closed";
  */
 
 std::cout<<"\n JD open \n";

  ofstream jd;
  ofstream hd;
  ofstream TD;
  TD.open("../../../data/2D.csv");
  jd.open("//c/Users/Jan/Desktop/1d_vertical.csv");
  hd.open("//c/Users/Jan/Desktop/1d_horizontal.csv");
  std::cout<<"\n JD opened \n"<<fg->yBegin()<<" "<<fg->yEnd()
           <<" "<<fg->xBegin()<<" "<<fg->xEnd()<<" \n" ;
  jd<<"z;r;Br;Bphi;Bz;Jr;Jphi;Jz;Lr;Lz;p;tau;rho;w\n";
  
 for(int iz=fg->yBegin();iz<fg->yEnd();iz++)
    for(int ir=fg->xBegin();ir<fg->xEnd();ir++)
    { 
      //std::cout<<"\n IZ "<<iz<<" IR"<<ir<<"\n";
      if (ir==100) 
      {
       
       jd<<iz*vCellSize<<";"<<ir*hCellSize<<";"<<B_r(ir,iz)<<";"<<B_phi(ir,iz)<<";"<<B_z(ir,iz)<<";"
         <<j_r(ir,iz)<<";"<<j_phi(ir,iz)<<";"<<j_z(ir,iz)<<";"<<L_r(ir,iz)<<";"<<L_z(ir,iz)<<";"
         <<p(ir,iz)<<";"<<tau(ir,iz)<<";"<<rho(ir,iz)<<";"<<w(ir,iz)
         <<"\n";
         
      }
      if (iz==50) //&&(ir<10)) 
      {
       
       hd<<iz<<";"<<ir<<";"<<B_r(ir,iz)<<";"<<B_phi(ir,iz)<<";"<<B_z(ir,iz)<<";"
         <<j_r(ir,iz)<<";"<<j_phi(ir,iz)<<";"<<j_z(ir,iz)<<";"<<L_r(ir,iz)<<";"<<L_z(ir,iz)<<";"
         <<p(ir,iz
         )<<";"<<tau(ir,iz)<<";"<<rho(ir,iz)<<";"<<w(ir,iz)
         <<"\n";
         //std::cout<<"pressure"<<p(ir,iz)<< "\n";
      }

      if (output_csv_2D==true)
      {
       
       TD<<B_r(ir,iz)<<";"<<B_phi(ir,iz)<<";"<<B_z(ir,iz)<<";"
         <<p(ir,iz)<<";"<<rho(ir,iz)<<"\n";
         
      }
    }
    
  
 jd.close();
 hd.close();
std::cout<<"\n JD closed \n";
 



if (output_csv==true)
{
 ofstream mhdeal("/r/mhs_data.csv", ios::out | ios::ate);
  real_t density;
  real_t pressure;
  real_t Bx;
  real_t By;
  real_t Bz;

   for(int iz=0;iz<box->zSize();iz++)
    {
      real_t z=iz*box->dz();

      for(int iy=0;iy<box->ySize();iy++)
	{
	  real_t y=iy*box->dy()-y0;

	  for(int ix=0;ix<box->xSize();ix++)
	    {
	      real_t x=ix*box->dx()-x0;
        
        
	      real_t r=sqrt(x*x+y*y);
	      real_t r_inv=(real_t)1.0/r;
        
        
        
	       density=(rho.atc(r,z));
                
		     pressure=(p.atc(r,z));
         
        
	      // Vector of magnetic field: Tranform to cartesian components
	       Bx=(r_inv*(-y*B_phi.atc(r,z)+x*B_r.atc(r,z)));
         
	       By=(r_inv*(x*B_phi.atc(r,z)+y*B_r.atc(r,z)));
	       Bz=(B_z.atc(r,z));
       mhdeal << ix<<"; "
              << iy<<"; "
              << iz<<"; "
              << density<<"; "
              << pressure<<"; "
              << Bx<<"; "
              << By<<"; "
              << Bz<<"\n ";
        




       
	    }

	}
    }
  
  mhdeal.close();  
}
if (output_binary==true)
{
 ofstream mhdeal("/r/mhs_data.out", ios::out | ios::binary);
  real_t density;
  real_t pressure;
  real_t Bx;
  real_t By;
  real_t Bz;

   for(int iz=0;iz<box->zSize();iz++)
    {
      real_t z=iz*box->dz();

      for(int iy=0;iy<box->ySize();iy++)
	{
	  real_t y=iy*box->dy()-y0;

	  for(int ix=0;ix<box->xSize();ix++)
	    {
	      real_t x=ix*box->dx()-x0;
        
        
	      real_t r=sqrt(x*x+y*y);
	      real_t r_inv=(real_t)1.0/r;
        
        
        
	       density=(rho.atc(r,z));
                
		     pressure=(p.atc(r,z));
         
        
	      // Vector of magnetic field: Tranform to cartesian components
	       Bx=(r_inv*(-y*B_phi.atc(r,z)+x*B_r.atc(r,z)));
         
	       By=(r_inv*(x*B_phi.atc(r,z)+y*B_r.atc(r,z)));
	       Bz=(B_z.atc(r,z));
         //Bx=x;
         
	       //By=y;
	       //Bz=z;
         mhdeal.write(reinterpret_cast<const char*>(&density),sizeof(real_t));
         mhdeal.write(reinterpret_cast<const char*>(&pressure),sizeof(real_t));
         mhdeal.write(reinterpret_cast<const char*>(&Bx),sizeof(real_t));
         mhdeal.write(reinterpret_cast<const char*>(&By),sizeof(real_t));
         mhdeal.write(reinterpret_cast<const char*>(&Bz),sizeof(real_t));
        




       
	    }

	}
    }
  
  mhdeal.close(); 
  std::cout<<"mhdeal out \n";
}

 //------- Fill the 3D cartesian box
if (output3D==true)

 { 
  std::cout<<"Working on 3D output \n";
  R3VectorField BB(box);
   
  RScalarField dens(box);
  RScalarField pressure(box);
  RScalarField temperature(box);
  
  R3VectorField JJ(box);
  RScalarField RBeta_3D(box);
  RScalarField A_3D(box);


  
  //float Fdens [box->xSize()][box->ySize()][box->zSize()];
  //float help;
  
  for(int iz=0;iz<box->zSize();iz++)
    {
      real_t z=iz*box->dz();

      for(int iy=0;iy<box->ySize();iy++)
	{
	  real_t y=iy*box->dy()-y0;

	  for(int ix=0;ix<box->xSize();ix++)
	    {
	      real_t x=ix*box->dx()-x0;
        
        
	      real_t r=sqrt(x*x+y*y);
	      real_t r_inv=(real_t)1.0/r;
        
        /*if (std::abs((B_phi.atc(r,z)-B_phi.at(r,z)))>1.)   //Testing difference between interpolations
          { 
            std::cout<<"err at (r,z)("<<r<<","<<z<<") \n";  
            cin.get();
          }
        */
        
        
        //Fdens[ix][iy][iz]=(float) rho.atc(r,z);
        
	      dens(ix,iy,iz)=ReverseFloat(rho.atc(r,z));
        temperature(ix,iy,iz)=ReverseFloat(w.atc(r,z));
        
		    pressure(ix,iy,iz)=ReverseFloat(p.atc(r,z));
        RBeta_3D(ix,iy,iz)=ReverseFloat(RBeta.atc(r,z));
	      // Vector of magnetic field: Tranform to cartesian components
	      BB(ix,iy,iz)[1]=ReverseFloat(r_inv*(-y*B_phi.atc(r,z)
						  +x*B_r.atc(r,z)));
         
	      BB(ix,iy,iz)[2]=ReverseFloat(r_inv*(x*B_phi.atc(r,z)
						  +y*B_r.atc(r,z)));
	      BB(ix,iy,iz)[3]=ReverseFloat(B_z.atc(r,z));
        // Vector of current density: Tranform to cartesian components
	      JJ(ix,iy,iz)[1]=ReverseFloat(r_inv*(-y*j_phi.atc(r,z)
						  +x*B_r.atc(r,z)));
	      JJ(ix,iy,iz)[2]=ReverseFloat(r_inv*(x*j_phi.atc(r,z)
						  +y*B_r.atc(r,z)));
	      JJ(ix,iy,iz)[3]=ReverseFloat(j_z.atc(r,z));
        A_3D(ix,iy,iz)=ReverseFloat(A.atc(r,z));

        




       
	    }

	}
    }

 
/*
//------- Write the VTK file

std::cout<<"3D output float";
ofstream out("/r/mhs_data.vtk");

 //--- Header
  out<<"# vtk DataFile Version 2.0"<<endl
     <<"MHD data"<<endl
     <<"BINARY"<<endl
     <<"DATASET STRUCTURED_POINTS"<<endl
     <<"DIMENSIONS "<<box->xSize()<<" "<<box->ySize()<<" "<<box->zSize()<<endl
     <<"ORIGIN "<<0<<" "<<0<<" "<<0<<endl
     <<"SPACING "<<box->dx()<<" "<<box->dy()<<" "<<box->dz()<<endl
     <<"POINT_DATA "<<box->gridExtent()<<endl<<flush;
//--- Data

  // 1. density

    size_t dataSize=box->gridExtent()*sizeof(float);
    
  out<<"SCALARS rho float 1"<<endl
     <<"LOOKUP_TABLE default"<<endl<<flush;

  out.write((const char*)Fdens,dataSize);

  out.flush();
out.close();

*/
  //write input for MHDeal



  //------- Write the VTK file

std::cout<<"3D output ";
ofstream out("/r/mhs_data.vtk");

  
 
  

  


  //--- Header
  out<<"# vtk DataFile Version 2.0"<<endl
     <<"MHD data"<<endl
     <<"BINARY"<<endl
     <<"DATASET STRUCTURED_POINTS"<<endl
     <<"DIMENSIONS "<<box->xSize()<<" "<<box->ySize()<<" "<<box->zSize()<<endl
     <<"ORIGIN "<<0<<" "<<0<<" "<<0<<endl
     <<"SPACING "<<box->dx()<<" "<<box->dy()<<" "<<box->dz()<<endl
     <<"POINT_DATA "<<box->gridExtent()<<endl<<flush;


  //--- Data

  // 1. density

    size_t dataSize=box->gridExtent()*sizeof(real_t);
    
  out<<"SCALARS rho double 1"<<endl
     <<"LOOKUP_TABLE default"<<endl<<flush;

  out.write((const char*)dens.data(),dataSize);

  out.flush();

  
  //out.flush();
  // 2. magnetic field

  out<<"VECTORS magField double"<<endl;

  out.write((const char*)BB.data(),3*dataSize);
  out.flush();
 
    // 3. pressure
  out<<"SCALARS pressure double 1"<<endl<<"LOOKUP_TABLE default"<<endl<<flush;
  out.write((const char*)pressure.data(),dataSize);
  out.flush();  

  // 4. Current density
  out<<"VECTORS J double"<<endl;
  out.write((const char*)JJ.data(),3*dataSize);
  //5. temperature (w)
  out<<"SCALARS temperature double 1"<<endl<<"LOOKUP_TABLE default"<<endl<<flush;
  out.write((const char*)temperature.data(),dataSize);
  out.flush();

  

  //6. beta
  out<<"SCALARS Beta double 1"<<endl<<"LOOKUP_TABLE default"<<endl<<flush;
  out.write((const char*)RBeta_3D.data(),dataSize);
  out.flush();

   // 1. A_phi

   
    
  out<<"SCALARS A_phi double 1"<<endl
     <<"LOOKUP_TABLE default"<<endl<<flush;

  out.write((const char*)A_3D.data(),dataSize);

  out.flush();

  
  out.close();
  
  std::cout<<"done";

 } 
 else std::cout<<"no 3D output \n";

/////////////////////////// FLOAT









//std::cout<<"\n size of: \n BB: "<<sizeof(BB(30,30,30)[1])<<" BB2: "<<sizeof(BB2(30,30,30)[1]);//<<" BBD:"<<sizeof(BBD(1,1,1)[1]);
//std::cout<<"\n LongDouble"<<sizeof(long double)<<" double"<<sizeof(double)<<" float"<<sizeof(float)<<"\n BB2(10,10,10)"<<BB2(10,10,10);
//size_t dataSizeD=box->gridExtent()*sizeof(float);

   // current date/time based on current system
   time_t now = time(0);
   
   // convert now to string form
   char* real_time = ctime(&now);

   
  return 0;
}
