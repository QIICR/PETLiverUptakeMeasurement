/*==============================================================================
 
 Program: PETLiverUptakeMeasurement
 
 (c) Copyright University of Iowa All Rights Reserved.
 
 See COPYRIGHT.txt
 or http://www.slicer.org/copyright/copyright.txt for details.
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 
 ==============================================================================*/
 
 /**
\file	itkSegmentBrainFilter.h
\brief	Segments the brain of a fairly specific variety of image.
*/

/*
*itkSegmentBrainFilter.h
*	Segments the brain of a fairly specific variety of image.
*Programmer: Markus Van Tol
*Date: 8/18/11
*
*The input image should be a PET scan of orientation: Brain towards higher Z, right side of body toward lower X
*
*/




#ifndef __itkSegmentBrainFilter_h
#define __itkSegmentBrainFilter_h


#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkGrayscaleFillholeImageFilter.h"

namespace itk
{

/** \class SegmentBrainFilter
 *  \brief Segments the brain of a fairly specific variety of image.
 *  \date 9/13/2011
 *  \author Markus Van Tol
 * A general filter for the specific segmentBrain.m MATLAB program.
 * The input image should be a PET scan of orientation: Brain towards higher Z, right side of body toward lower X
 *
 * The class is templated over the type of the input and output images.
 *
 * Template parameters for class SegmentBrainFilter:
 *
 * - TInputImage = The image type that is going to be segmented.
 * - TOutputImage = The image type that the segmentation will take.
 *
 * \ingroup SegmentationFilters
 */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT SegmentBrainFilter : public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SegmentBrainFilter                         Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  ITK_DISALLOW_COPY_AND_ASSIGN(SegmentBrainFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SegmentBrainFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType;
  typedef typename    InputImageType::PixelType  InputImagePixelType;

  typedef TOutputImage                              OutputImageType;
  typedef typename     OutputImageType::Pointer     OutputImagePointer;
  typedef typename     OutputImageType::RegionType  OutputImageRegionType;
  typedef typename     OutputImageType::PixelType   OutputImagePixelType;

  typedef typename TOutputImage::PointType PointType;
  

  /** ImageDimension enumeration */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  itkSetMacro(UpperThreshold, InputImagePixelType);
  itkGetMacro(UpperThreshold, InputImagePixelType);
  itkSetMacro(LowerThreshold, InputImagePixelType);
  itkGetMacro(LowerThreshold, InputImagePixelType); 
  itkGetMacro(Centroid, PointType);
  itkGetMacro(BoundarySize, PointType);
  itkGetMacro(BoundaryStart, PointType);

  itkSetMacro(MinimumVolume, float);
  itkGetMacro(MinimumVolume, float);
  //itkSetMacro(Verbose, bool);
  //itkGetMacro(Verbose, bool);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ImageDimensionCheck,
      (Concept::SameDimensionOrMinusOne<itkGetStaticConstMacro(InputImageDimension),
                                        itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif



protected:
  SegmentBrainFilter() = default;
  ~SegmentBrainFilter() override = default;
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Apply changes to the output image information. */
  void GenerateOutputInformation() override;

  /** Apply changes to the input image requested region. */
  void GenerateInputRequestedRegion() override;

  /** This method implements the actual creation of the diagonal image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void GenerateData(void) override;

private:
  void FindCentroid(OutputImagePointer outputImage);
  bool MinimumVolume(OutputImagePixelType label, OutputImagePointer image);
  PointType m_Centroid;
  PointType m_BoundarySize;
  PointType m_BoundaryStart;
  OutputImagePixelType m_LowerThreshold{ NumericTraits<OutputImagePixelType>::min() };
  OutputImagePixelType m_UpperThreshold{ NumericTraits<OutputImagePixelType>::max() };
  float m_MinimumVolume;
  //bool m_Verbose;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSegmentBrainFilter.txx"
#endif

#endif
