/* 
   class TField

   A template class for physical field on cartesian geometry
   Generic + some specific methods for 1-, 2- and 3-D fields are defined 

   class TFieldAdaptive

   Generalisation for fields on AMR geometries

   20/11/2007 
   barta@asu.cas.cz

*/
 

#ifdef _USES_LIB

#pragma implementation "field"

#include "field.h"

#endif //--> _USES_LIB

#include <stdlib.h>
#include <fstream>
using namespace std;


/*****************************************************************************
                                class TField
*****************************************************************************/

template<typename TData>

#ifndef _USES_LIB
  inline 
#endif

TField<TData>::TField(const TField<TData>& orig)
{
  _geometry=orig.geometry();

  _data=new TData[_geometry->gridExtent()];
 
  for(register int i=0;i<_geometry->gridExtent();i++)
    *(_data+i)=*(orig._data+i);

}


/****************************************************************************/


template<typename TData>

#ifndef _USES_LIB
  inline 
#endif

TField<TData>& TField<TData>::operator=(const TField<TData>& orig)
{
  if(_geometry && (_geometry->gridExtent()!=orig._geometry->gridExtent()))
    {
      if(_data) delete[] _data;
      _data=new TData[orig._geometry->gridExtent()];
    }
 
  _geometry=orig.geometry();

  for(register int i=0;i<_geometry->gridExtent();i++)
    *(_data+i)=*(orig._data+i);

  return *this;
}

/****************************************************************************/

// Set all elements in the field to prescribed value
template<typename TData>

#ifndef _USES_LIB
  inline 
#endif

TField<TData>& TField<TData>::operator=(const TData& value)
{
  for(int i=0;i<_geometry->gridExtent();i++)
    *(_data+i)=value;

  return *this;
}


/****************************************************************************/

// Fast swap operator
template<typename TData>

#ifndef _USES_LIB
  inline 
#endif

// Only the field with the same structure (geometry) can be swapped!
void TField<TData>::swapWith(TField<TData>& that)
{

#ifdef _DEBUG /////////////////// DEBUG //////////////////////////////////////

  if(this->_geometry!=that._geometry)
    {
      cout<<"Swap fields: geometric structure may dissagree!"<<endl;
    }

#endif ////////////////////////// !DEBUG /////////////////////////////////////

  TData* d;

  d=this->_data;
  this->_data=that._data;
  that._data=d;

  return;
}


/****************************************************************************/

// linear interpolation operator on 1D fields

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

TData TField<TData>::at(real_t x) const 
{
  static const real_t dxInv=1.0/_geometry->dx();

  x*=dxInv;

  real_t xLow=floor(x);
  real_t xCell=x-xLow;

  int ix=(int) xLow;

  const TField<TData>& current=*this;

  return current(ix)+xCell*(current(ix+1)-current(ix));
}


/****************************************************************************/

// bi-linear interpolation operator on 2D fields

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

TData TField<TData>::at(real_t x,real_t y) const 
{
  static const real_t invDx=1.0/_geometry->dx();
  static const real_t invDy=1.0/_geometry->dy();

  x*=invDx;
  y*=invDy;

  real_t xLow=floor(x);
  real_t yLow=floor(y);

  real_t xCell=x-xLow;
  real_t yCell=y-yLow;

  int ix=(int) xLow;
  int iy=(int) yLow;

  const TField<TData>& current=*this;

  return current(ix,iy)+xCell*(current(ix+1,iy)-current(ix,iy))
    +yCell*(current(ix,iy+1)+xCell*(current(ix+1,iy+1)-current(ix,iy+1)));
}

/****************************************************************************/

// bi-cubic interpolation operator on 2D fields

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

TData TField<TData>::atc(real_t x,real_t y) const 
{
  /*
    Uses factorized form of Hermite cubic spline interpolation 
    A.K.A. bi-cubic convolution interpolation
    
    See http://en.wikipedia.org/wiki/Bicubic_interpolation
  */

  static const real_t invDx=1.0/_geometry->dx();
  static const real_t invDy=1.0/_geometry->dy();

  x*=invDx;
  y*=invDy;

  real_t xLow=floor(x);
  real_t yLow=floor(y);

  real_t xCell=x-xLow;
  real_t yCell=y-yLow;

  int ix=(int) xLow;
  int iy=(int) yLow;

  const TField<TData>& current=*this;
  TData cxVals[4],cyVals[4];

  for(int yy=-1;yy<3;yy++) // loop over adjacent lines
    {
      for(int xx=-1;xx<3;xx++) // fill the coarse-valued array of neighbours
	cxVals[xx+1]=current(ix+xx,iy+yy);

      cyVals[yy+1]=hcInterp(cxVals,xCell); // interpolate one line
    }

  return hcInterp(cyVals,yCell); // interpolate the column
}


/****************************************************************************/

// tri-linear interpolation operator on 3D fields

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

TData TField<TData>::at(real_t x,real_t y, real_t z) const 
{
  static const real_t invDx=1.0/_geometry->dx();
  static const real_t invDy=1.0/_geometry->dy();
  static const real_t invDz=1.0/_geometry->dz();

  x*=invDx;
  y*=invDy;
  y*=invDz;

  real_t xLow=floor(x);
  real_t yLow=floor(y);
  real_t zLow=floor(z);

  real_t xCell=x-xLow;
  real_t yCell=y-yLow;
  real_t zCell=z-zLow;

  int ix=(int) xLow;
  int iy=(int) yLow;
  int iz=(int) zLow;

  const TField<TData>& current=*this;

  return current(ix,iy,iz)+xCell*(current(ix+1,iy,iz)-current(ix,iy,iz))
    +yCell*(current(ix,iy+1,iz)
	    +xCell*(current(ix+1,iy+1,iz)-current(ix,iy+1,iz)))
    +zCell*(current(ix,iy,iz+1)
	    +xCell*(current(ix+1,iy,iz+1)-current(ix,iy,iz+1))
	    +yCell*(current(ix,iy+1,iz+1)
		    +xCell*(current(ix+1,iy+1,iz+1)-current(ix,iy+1,iz+1))));
}


/****************************************************************************/

// generalisation of the above methods

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

TData TField<TData>::at(const lampa::RVector& r) const 
{
  switch (r.dim())
    {
    case 1: return this->at(r[0]);
    case 2: return this->at(r[0],r[1]);
    case 3: return this->at(r[0],r[1],r[2]);
    
    default: ::exit(-1);
    }
}

/****************************************************************************/

// read data from file

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

int TField<TData>::read(ifstream& inp, bool all)
{
  if(all) // read all the data, including boundaries
    inp.read((char*) _data, geometry()->gridExtent()*elementSize());
  else    // read only interior
    {
      switch(geometry()->dim())
	{
	case 1:
	  inp.read((char*)_data+geometry()->gridPointAt(0),
		   geometry()->xSize()*elementSize()); break;

	case 2:
	  for(int y=0;y<geometry()->ySize();y++)
	    inp.read((char*)_data+geometry()->gridPointAt(0,y),
		     geometry()->xSize()*elementSize()); break;

	case 3:
	  for(int z=0;z<geometry()->zSize();z++)
	    for(int y=0;y<geometry()->ySize();y++)
	      inp.read((char*)_data+geometry()->gridPointAt(0,y,z),
		       geometry()->xSize()*elementSize()); break;

#ifdef _DEBUG /////////////////// DEBUG //////////////////////////////////////

	default:
	  cout<<"TField::read(): only dimenions 1,2, or 3 are supported"
	      <<endl;

#endif ////////////////////////// !DEBUG /////////////////////////////////////

	}
    }
     
  return inp.rdstate();
}


/****************************************************************************/

// write data to file

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

int TField<TData>::write(ofstream& out, bool all)
{
  if(all) // read all the data, including boundaries
    out.write((const char*) _data, geometry()->gridExtent()*elementSize());
  else    // read only interior
    {
      switch(geometry()->dim())
	{
	case 1:
	  out.write((const char*)_data+geometry()->gridPointAt(0),
		    geometry()->xSize()*elementSize()); break;

	case 2:
	  for(int y=0;y<geometry()->ySize();y++)
	    out.write((const char*)_data+geometry()->gridPointAt(0,y),
		      geometry()->xSize()*elementSize()); break;

	case 3:
	  for(int z=0;z<geometry()->zSize();z++)
	    for(int y=0;y<geometry()->ySize();y++)
	      out.write((const char*)_data+geometry()->gridPointAt(0,y,z),
			geometry()->xSize()*elementSize()); break;

#ifdef _DEBUG /////////////////// DEBUG //////////////////////////////////////

	default:
	  cout<<"TField::write(): only dimenions 1,2, or 3 are supported"
	      <<endl;

#endif ////////////////////////// !DEBUG /////////////////////////////////////

	}
    }
     
  return out.rdstate();
}



/****************************************************************************/
// MPI version of I/O routines - has a sense only in paralel programs


#ifdef _USES_MPI ////////////////////  MPI ///////////////////////////////


/****************************************************************************/

// read data from MPI file

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

MPI_Status TField<TData>::mpiRead(MPI_File& inp, MPI_Datatype& MPI_dataStruct)
{
  MPI_Status opStatus;

  MPI_File_read_all(inp,_data,1,MPI_dataStruct,&opStatus);


  return opStatus;
}


/****************************************************************************/

// write data to MPI file

template<class TData>

#ifndef _USES_LIB
  inline 
#endif

MPI_Status TField<TData>::mpiWrite(MPI_File& out, MPI_Datatype& MPI_dataStruct)
{
  MPI_Status opStatus;

  MPI_File_write_all(out,_data,1,MPI_dataStruct,&opStatus);

  return opStatus;
}


#endif ////////////////////////////// !MPI ///////////////////////////////




/******* Hermite spline cubic interpolation. Used for mesh refinement ********/

template<typename TData>
TData hcInterp(TData* coarseVal,real_t iFraction)
{
  TData iResult(0.0);
  real_t iFractPow=1.0;
  TData iFactor[4];

  iFactor[0]=coarseVal[1];
  iFactor[1]=0.5*(coarseVal[2]-coarseVal[0]);
  iFactor[2]=coarseVal[0]-2.5*coarseVal[1]+2.0*coarseVal[2]-0.5*coarseVal[3];
  iFactor[3]=1.5*(coarseVal[1]-coarseVal[2])+0.5*(coarseVal[3]-coarseVal[0]);

  for(int i=0;i<4;i++)
    {
      iResult+=iFractPow*iFactor[i];
      iFractPow*=iFraction;
    }

  return iResult;
}


/*****************************************************************************
                           class TFieldAdaptive
****************************************************************************/


template<typename TData>
TFieldAdaptive<TData>::TFieldAdaptive(const TGeometry* geom):
  TField<TData>(geom),
  _root(0),
  _parent(0),
  _parCoords(0),
  _zoomDepth(0)
{
  _zoomed=new TFieldAdaptive* [this->_geometry->gridExtent()];

  for(register int i=0;i<this->_geometry->gridExtent();i++)
    _zoomed[i]=0;

  _nZoomed=0;
}


template<typename TData>
TFieldAdaptive<TData>::~TFieldAdaptive()
{
  if(_zoomed!=0)
    {
      // delete all children
      for(register int i=0;i<this->_geometry->gridExtent();i++)
	if(_zoomed[i]!=0)
	  delete _zoomed[i];

      // delete list of children
      delete _zoomed;
    }

  if(_parCoords!=0)
    delete _parCoords;
}


template<typename TData> TFieldAdaptive<TData>*
TFieldAdaptive<TData>::zoomIn(int x)
{  
  TZoomedGeometry* geomGlob=(TZoomedGeometry*)this->_geometry;

  TFieldAdaptive* subField=_zoomed[geomGlob->gridPointAt(x)]=
    new TFieldAdaptive(geomGlob->atZoomDepth(_zoomDepth+1));

  subField->_parent=this;
  subField->_root=((_zoomDepth==0)?this:_root);
  subField->_zoomDepth=_zoomDepth+1;
  
  // to be completed!

  return subField;
}


template<typename TData> TFieldAdaptive<TData>*
TFieldAdaptive<TData>::zoomIn(int x,int y,fill_t fill)
{  
  TZoomedGeometry* geomGlob=(TZoomedGeometry*)this->_geometry;
  const TGeometry* geomLoc=geomGlob->atZoomDepth(_zoomDepth+1);

  //create space for new subfield

  TFieldAdaptive* subField=new TFieldAdaptive(geomLoc);

  subField->_parent=this;

  subField->_root=((_zoomDepth==0)?this:_root);
  subField->_zoomDepth=_zoomDepth+1;

  subField->_parCoords=new int[geomGlob->dim()];
  subField->_parCoords[0]=x;
  subField->_parCoords[1]=y;
  
  if(fill!=noFill)
    {
      //create grid in the subfield

      real_t dx=geomGlob->dx();
      real_t dy=geomGlob->dy();
      
      real_t ddx=geomLoc->dx();
      real_t ddy=geomLoc->dy();
      
      real_t xCellCenter=(real_t)x*dx;
      real_t yCellCenter=(real_t)y*dy;
      
      real_t xCellOrig=xCellCenter-0.5*(dx-ddx);
      real_t yCellOrig=yCellCenter-0.5*(dy-ddy);
      
      real_t xNewGridPoint,yNewGridPoint;

      //fill it with the interpolated values

      if(fill==linInterp) // (bi)linear interpolation
	{
	  for(int iy=geomLoc->yBegin();iy<geomLoc->yEnd();iy++)
	    for(int ix=geomLoc->xBegin();ix<geomLoc->xEnd();ix++)
	      {
		xNewGridPoint=xCellOrig+(real_t)ix*ddx;
		yNewGridPoint=yCellOrig+(real_t)iy*ddy;

		(*subField)(ix,iy)=this->at(xNewGridPoint,yNewGridPoint);
	      }
	}
      else // (bi)cubic interpolation
	{
	  for(int iy=geomLoc->yBegin();iy<geomLoc->yEnd();iy++)
	    for(int ix=geomLoc->xBegin();ix<geomLoc->xEnd();ix++)
	      {
		xNewGridPoint=xCellOrig+(real_t)ix*ddx;
		yNewGridPoint=yCellOrig+(real_t)iy*ddy;

		(*subField)(ix,iy)=this->atc(xNewGridPoint,yNewGridPoint);
	      }
	}

    }

 
 //register new subfield in the parent's list(s) of children 

  _zoomed[geomGlob->gridPointAt(x,y)]=subField;
  _zList.push_back(geomGlob->gridPointAt(x,y));
  _nZoomed++;

  return subField;
}


template<typename TData> TFieldAdaptive<TData>*
TFieldAdaptive<TData>::zoomIn(int x,int y,int z)
{  
  TZoomedGeometry* geomGlob=(TZoomedGeometry*)this->_geometry;
  
  TFieldAdaptive* subField=_zoomed[geomGlob->gridPointAt(x,y,z)]=
    new TFieldAdaptive(geomGlob->atZoomDepth(_zoomDepth+1));

  subField->_parent=this;
  subField->_root=((_zoomDepth==0)?this:_root);
  subField->_zoomDepth=_zoomDepth+1;
  
  //to be completed!

  return subField;
}

template<typename TData> void
TFieldAdaptive<TData>::zoomOut(int x,int y)
{  
  int gridPointID=this->geometry()->gridPointAt(x,y);
  TFieldAdaptive* theSubField=_zoomed[gridPointID];

  // de-register the subfield from list of children

  _zoomed[gridPointID]=0;
  _zList.remove(gridPointID);
  _nZoomed--;

  // delete sub-field

  delete theSubField;
  
  return;
}

template<typename TData> TData
TFieldAdaptive<TData>::meanValue() const
{  
  const TGeometry* geomLoc=this->_geometry;

  TData integ(0);

  //integrate over zoomed cell (subfield)

  for(int iy=0;iy<geomLoc->ySize();iy++)
    for(int ix=0;ix<geomLoc->xSize();ix++)
      integ+=(*this)(ix,iy);

  integ*=1.0/((real_t)geomLoc->xSize()*(real_t)geomLoc->xSize());

  return integ;
}


#ifdef _USES_LIB

// Template instances with predefined parameters

template class TField<real_t>;
template class TField<std::complex<real_t> >;

template class TField<lampa::PVector2R>;
template class TField<lampa::PVector3R>;

template class TField<lampa::PVector2C>;
template class TField<lampa::PVector3C>;

template class TField<lampa::SVector<real_t,6> >;
template class TField<lampa::SVector<real_t,8> >;

template class TFieldAdaptive<real_t>;


#endif // _USES_LIB

