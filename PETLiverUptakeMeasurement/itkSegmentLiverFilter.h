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
\file itkSegmentLiverFilter.h
\brief Class definition for segmenting the liver of a PET scan image of a certain orientation.
*/
/*
*itkSegmentLiverFilter.h
*Can be used to create a sphere within the liver of an appropriate image.
*Programmer: Markus Van Tol
*Date: 1/6/12
*
*Requires either the location of the brain and distances from it or the location it would get from them.
*Makes a sphere of the specified radius at that location.  The sphere has value 1.  All other pixels have value 0.
*Also identifies a physical bounding box around the centroid to contain the whole segmentation.
*/




#ifndef __itkSegmentLiverFilter_h
#define __itkSegmentLiverFilter_h


#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"

//#include "itkIsolatedConnectedImageFilter.h"


namespace itk
{

/** \class SegmentLiverFilter
 *  \brief Segments the liver of a fairly specific variety of image.
 *  \date 1/6/2012
 *  \author Markus Van Tol
 * A general filter for the specific segmentLiver.m MATLAB program.
 * The input image should be a PET scan of orientation: Brain towards higher Z, right side of body toward lower X
 * Requires a centroid of the brain and a distance values from it, or the interior corner poitn of a bounding box to search for the liver in.  Can use physical or voxel coordinates for these.
 * Makes a sphere (physically or voxel-wise) of the specified radius at that location.  The sphere has value 1.  All other pixels have value 0.
 * Now also erodes the final image by a chosen value, based on the Z-direction spacing.
 * Also identifies a physical bounding box around the centroid to contain the whole segmentation.
 *
 * The class is templated over the type of the input and output images.
 *
 * Template parameters for class SegmentLiverFilter:
 *
 * - TInputImage = The image type that is going to be segmented.
 * - TOutputImage = The image type that the segmentation will take.
 *
 * \ingroup SegmentationFilters
 */


template <class TInputImage, class TOutputImage>
class ITK_EXPORT SegmentLiverFilter : public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SegmentLiverFilter                         Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  ITK_DISALLOW_COPY_AND_ASSIGN(SegmentLiverFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SegmentLiverFilter, ImageToImageFilter);

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

  itkSetMacro(LowerThreshold, InputImagePixelType);
  itkGetMacro(LowerThreshold, InputImagePixelType);
  itkSetMacro(UpperThreshold, InputImagePixelType);
  itkGetMacro(UpperThreshold, InputImagePixelType);  




  /* Set the lower threshold for segmentation */
  //void SetLowerThreshold(InputImagePixelType newThresh)	{	m_LowerThreshold = newThresh;	}

  /* Set the upper threshold for segmentation */
  //void SetUpperThreshold(InputImagePixelType newThresh)	{	m_UpperThreshold = newThresh;	}

  /* Get the lower threshold for segmentation */
  //InputImagePixelType GetLowerThreshold()	{	return m_LowerThreshold;	}

  /* Get the upper threshold for segmentation */
  //InputImagePixelType GetUpperThreshoold()	{	return m_LowerThreshold;	}


  //itkSetMacro(BrainCentroid, PointType);
  //itkSetMacro(SectorLocation, PointType);
  itkGetMacro(BrainCentroid, PointType);
  itkSetMacro(UseBrainCentroid, bool);
    
  /** Set the ROI for the liver region */
    void SetSearchROI(InputImageRegionType searchROI) {m_SearchROI = searchROI; this->m_SetUseSearchROI=true; };
    
  /** Set the physical centroid of the brain */
  void SetBrainCentroid(PointType newCentroid)	{	m_BrainCentroid = newCentroid;	this->SetUseBrainCentroid(true);	}

  /** Set the interior corner of the bounding box directly, as a physical coordinate */
  void SetSectorLocation(PointType newLocation)	{	m_SectorLocation = newLocation;	this->SetUseBrainCentroid(false);	}




  itkSetMacro(RealSpacing, bool);
  itkSetMacro(SpacingDistance, bool);

  /** Set the image to use voxel-based spacing for the sphere and bounding box values */
  void UseVoxelSpacing()	{	this->SetRealSpacing(false);	}

  /** Set the image to use physical spacing for the sphere and bounding box values */
  void UseRealSpacing()		{	this->SetRealSpacing(true);	}

  /** Set the distance from the brain to be physical */
  void UseSpacedDistance()	{	this->SetSpacingDistance(false);	}

  /** Set the distance from the brain and center to be voxel-based */
  void UseVoxelDistance()	{	this->SetSpacingDistance(true);	}

  itkSetMacro(Erosion, float);
  itkGetMacro(Erosion, float);

  itkGetMacro(Radius, double);


  itkSetMacro(MinZFromBrain, double);
  /* Set the minimum Z value below the brain */
  //void SetMinZFromBrain(double newMinZ)		{	m_MinZFromBrain = newMinZ;	}


  itkSetMacro(MinXFromBrain, double);
  /* Set the minimum X value away from the brain*/
  //void SetMinXFromBrain(double newMinX)		{	m_MinXFromBrain = newMinX;	}
  /** Input and output images must be the same dimension, or the output's
      dimension must be one less than that of the input. */
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ImageDimensionCheck,
      (Concept::SameDimensionOrMinusOne<itkGetStaticConstMacro(InputImageDimension),
                                        itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif
  
  /** Get the physical centroid of the liver, after generating data */
  PointType GetCentroid()	{	return m_LiverCentroid;	}

  /** Get the physical size of the liver segmentation seed, after generating data */
  PointType GetBoundarySize();

  /** Get the physical start of the liver segmentation seed, after generating data */
  PointType GetBoundaryStart();

protected:
  SegmentLiverFilter() = default;
  ~SegmentLiverFilter() override = default;
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
  PointType FindCentroid(OutputImagePointer outputImage);
  InputImageRegionType m_SearchROI;
  bool m_SetUseSearchROI{ false };
  PointType m_BrainCentroid;
  PointType m_SectorLocation;
  PointType m_LiverCentroid;	
  bool m_UseBrainCentroid{ false };
  bool m_RealSpacing{ false };
  double m_MinZFromBrain;
  double m_MinXFromBrain;
  OutputImagePixelType m_LowerThreshold{ NumericTraits<OutputImagePixelType>::min() };
  OutputImagePixelType m_UpperThreshold{ NumericTraits<OutputImagePixelType>::max() };
  bool m_SpacingDistance{ false };
  double m_Radius;
  float m_Erosion{ 10 };
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSegmentLiverFilter.txx"
#endif

#endif
