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
#ifndef __itkEdgePreservingSmoothingImageFilter_h
#define __itkEdgePreservingSmoothingImageFilter_h

#include "itkBoxImageFilter.h"
#include "itkImage.h"

namespace itk
{
/** \class EdgePreservingSmoothingImageFilter
 * \brief Applies an edge preserving filter to an image
 *
 * This filter requires that the input pixel type provides an operator<()
 * (LessThan Comparable).
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup IntensityImageFilters
 * \ingroup ITKSmoothing
 *
 * \wiki
 * \wikiexample{Smoothing/EdgePreservingSmoothingImageFilter,EdgePreservingSmoothing filter an image}
 * \endwiki
 */
template< typename TInputImage, typename TOutputImage >
class EdgePreservingSmoothingImageFilter:
  public BoxImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef EdgePreservingSmoothingImageFilter                                     Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(EdgePreservingSmoothingImageFilter, BoxImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( SameDimensionCheck,
                   ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
  itkConceptMacro( InputConvertibleToOutputCheck,
                   ( Concept::Convertible< InputPixelType, OutputPixelType > ) );
  itkConceptMacro( InputLessThanComparableCheck,
                   ( Concept::LessThanComparable< InputPixelType > ) );
  // End concept checking
#endif

protected:
  EdgePreservingSmoothingImageFilter();
  virtual ~EdgePreservingSmoothingImageFilter() {}

  /** EdgePreservingSmoothingImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId);

private:
  EdgePreservingSmoothingImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEdgePreservingSmoothingImageFilter.hxx"
#endif

#endif
