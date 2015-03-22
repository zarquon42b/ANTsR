/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkPolynomialConstrainedSlicewiseTranslationTransform_hxx
#define __itkPolynomialConstrainedSlicewiseTranslationTransform_hxx

#include "itkPolynomialConstrainedSlicewiseTranslationTransform.h"
#include "itkMath.h"

namespace itk
{
// Constructor with default arguments
template <typename TScalar, unsigned int NDimensions>
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::PolynomialConstrainedSlicewiseTranslationTransform() : Superclass(ParametersDimension),
  m_IdentityJacobian(NDimensions, NDimensions)
{
  m_Offset.Fill(0);

  // The Jacobian of this transform is constant.
  // Therefore the m_IdentityJacobian variable can be
  // initialized here and be shared among all the threads.
  this->m_IdentityJacobian.Fill(0.0);
  for( unsigned int i = 0; i < NDimensions; i++ )
    {
    this->m_IdentityJacobian(i, i) = 1.0;
    }
}

// Destructor
template <typename TScalar, unsigned int NDimensions>
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::
~PolynomialConstrainedSlicewiseTranslationTransform()
{
}

// Set the parameters
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>
::SetParameters(const ParametersType & parameters)
{
  // Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  bool modified = false;
  for( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    if( m_Offset[i] != parameters[i] )
      {
      m_Offset[i] = parameters[i];
      modified = true;
      }
    }
  if( modified )
    {
    this->Modified();
    }
}

// Get the parameters
template <typename TScalar, unsigned int NDimensions>
const typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::ParametersType
& PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>
::GetParameters(void) const
  {
  for( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    this->m_Parameters[i] = this->m_Offset[i];
    }
  return this->m_Parameters;
  }

// Print self
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Offset: " << m_Offset << std::endl;
}

// Compose with another affine transformation
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::Compose(const Self *other, bool)
{
  this->Translate(other->m_Offset);
}

// Compose with a translation
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::Translate(const OutputVectorType & offset, bool)
{
  ParametersType newOffset(SpaceDimension);

  for( unsigned int i = 0; i < SpaceDimension; i++ )
    {
    newOffset[i] = m_Offset[i] + offset[i];
    }
  this->SetParameters(newOffset);
}

// Transform a point
template <typename TScalar, unsigned int NDimensions>
typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::OutputPointType
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::TransformPoint(const InputPointType & point) const
{
  OutputPointType outPoint;
  outPoint[0] = point[0];
  outPoint[1] = point[1];
  outPoint[2] = point[2];
  for ( unsigned int i = 0; i <= this->m_PolynomialDegree; i++ )
    {
    ScalarType z = point[2];
    z = ( z - this->m_ZMean ) * this->m_ZRange;
    ScalarType zpow =  vcl_pow( z, static_cast<double>(i) );
    ScalarType paramx = this->m_Parameters[i];
    ScalarType paramy = this->m_Parameters[i + this->m_PolynomialDegree+1 ];
    outPoint[0] += paramx * zpow;
    outPoint[1] += paramy * zpow;
    }
  return outPoint;
}

// Transform a vector
template <typename TScalar, unsigned int NDimensions>
typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::OutputVectorType
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::TransformVector(const InputVectorType & vect) const
{
  return vect;
}

// Transform a vnl_vector_fixed
template <typename TScalar, unsigned int NDimensions>
typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::OutputVnlVectorType
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::TransformVector(const InputVnlVectorType & vect) const
{
  return vect;
}

// Transform a CovariantVector
template <typename TScalar, unsigned int NDimensions>
typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::OutputCovariantVectorType
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::TransformCovariantVector(const InputCovariantVectorType & vect) const
{
  return vect;
}

// return an inverse transformation
template <typename TScalar, unsigned int NDimensions>
bool
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::GetInverse(Self *inverse) const
{
  if( !inverse )
    {
    return false;
    }

  inverse->m_Offset   = -m_Offset;
  return true;
}

// Return an inverse of this transform
template <typename TScalar, unsigned int NDimensions>
typename PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::InverseTransformBasePointer
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>
::GetInverseTransform() const
{
  Pointer inv = New();

  return GetInverse(inv) ? inv.GetPointer() : ITK_NULLPTR;
}

// Compute the Jacobian in one position
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::ComputeJacobianWithRespectToParameters(
  const InputPointType &,
  JacobianType & jacobian) const
{
  // the Jacobian is constant for this transform, and it has already been
  // initialized in the constructor, so we just need to return it here.
  jacobian = this->m_IdentityJacobian;
}

// Compute jacobian with respect to position
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>
::ComputeJacobianWithRespectToPosition(const InputPointType &,
                                       JacobianType & jac) const
{
  jac.SetSize( NDimensions, NDimensions );
  jac.Fill(0.0);
  for( unsigned int dim = 0; dim < NDimensions; dim++ )
    {
    jac[dim][dim] = 1.0;
    }
}

// Set the parameters for an Identity transform of this class
template <typename TScalar, unsigned int NDimensions>
void
PolynomialConstrainedSlicewiseTranslationTransform<TScalar, NDimensions>::SetIdentity()
{
  if ( this->m_Parameters.Size() != this->GetNumberOfParameters() )
    {
    this->m_Parameters.SetSize( this->GetNumberOfParameters() );
    }
  this->m_Parameters.Fill( 0 );
  m_Offset.Fill(0.0);
}

} // namespace

#endif
