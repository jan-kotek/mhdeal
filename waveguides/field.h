/******************************************************************************

                   Main header file for libfield.so

libfield.so implements spatially-discrete representation of physical fields
encapsulating in one object the data (stored in continguous block of memory) 
with their geometrical structure, i.e. the library provides bi-directional 
mapping {physical coordinates} <--> {memory adress} 

Data can be of various types, for example the blocks of double, complex
MVector2D<double>, or MVector3D<complex> (implemented in LAMPa package) provide
real and complex scalar fields, field of real 2-vectors, and field of complex
3-vectros, respectivelly

Jan Kotek
jankotek@asu.cas.cz

Miroslav Barta
barta@asu.cas.cz

******************************************************************************/


#ifndef __LIB_FIELD__ 

#define __LIB_FIELD__


// mesh geometry / memory mapping
#include "field/geometry.h"

// {data,geometry} encapsulation
#include "field/field.h"


/********************* Frequently used predefined types **********************/

// Default floating-point number precision

#ifndef real_t
#define real_t double
#endif

#include "lampa.h"

// real and complex scalar fields

typedef TField<real_t> RScalarField;
typedef TField<std::complex<real_t> > CScalarField;

// real and complex 2D and 3D vector fields

typedef TField<lampa::PVector2R> R2VectorField;
typedef TField<lampa::PVector3R> R3VectorField;

typedef TField<lampa::PVector2C> C2VectorField;
typedef TField<lampa::PVector3C> C3VectorField;

// real and complex 2D and 3D tensor (2nd order) fields

typedef TField<lampa::R2SquareMatrix> R2TensorField;
typedef TField<lampa::R3SquareMatrix> R3TensorField;

typedef TField<lampa::C2SquareMatrix> C2TensorField;
typedef TField<lampa::C3SquareMatrix> C3TensorField;


#endif  // __LIB_FIELD__
