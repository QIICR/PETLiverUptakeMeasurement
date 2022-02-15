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
 
 
/*
*itkSegmentLiverFilter.txx
*Implementation of itkSegmentLiverFilter.h.
*
*Programmer: Markus Van Tol
*Date: 2/29/12
*
*
*/


#ifndef __itkSegmentLiverFilter_txx
#define __itkSegmentLiverFilter_txx

#include "itkSegmentLiverFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <iostream>
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include "itkRegionalMinimaImageFilter.h"
#include <itkImageFileWriter.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
void
SegmentLiverFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  itkDebugMacro("GenerateOutputInformation Start");

  typename TOutputImage::RegionType outputRegion;
  typename TOutputImage::SizeType  outputSize;
  typename TOutputImage::IndexType outputIndex;
  typename TOutputImage::SpacingType outSpacing;
  typename TOutputImage::PointType outOrigin;

  // Get pointers to the input and output
  auto output = this->GetOutput();
  auto input = const_cast< TInputImage * >( this->GetInput() );
  
  if( !input || !output )
    {
	  if (!input)	{itkDebugMacro("!input");}
	  if (!output) {itkDebugMacro("!output");}
    return;
    }
 
  auto inputIndex = input->GetLargestPossibleRegion().GetIndex();
  auto inputSize = input->GetLargestPossibleRegion().GetSize();
  auto inSpacing = input->GetSpacing();
  auto inOrigin = input->GetOrigin();

  // Set the LargestPossibleRegion of the output.
  for(unsigned int i = 0; i<InputImageDimension; i++)
    {
      outputSize[i]  = inputSize[i];
      outputIndex[i] = inputIndex[i];
      outSpacing[i] = inSpacing[i];
      outOrigin[i]  = inOrigin[i];

    }//end for

  outputRegion.SetSize(outputSize);
  outputRegion.SetIndex(outputIndex);
  output->SetOrigin(outOrigin);
  output->SetSpacing(outSpacing);
  output->SetLargestPossibleRegion(outputRegion);

  itkDebugMacro("GenerateOutputInformation End");
}//end GenerateOutputRegion


template <class TInputImage, class TOutputImage>
void
SegmentLiverFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  itkDebugMacro("GenerateInputRequestedRegion Start");
  Superclass::GenerateInputRequestedRegion();

  if ( this->GetInput() )
    {
    typename TInputImage::RegionType RequestedRegion;
    typename TInputImage::SizeType  inputSize;
    typename TInputImage::IndexType inputIndex;
    
    auto inputLargSize = this->GetInput()->GetLargestPossibleRegion().GetSize();
    auto inputLargIndex = this->GetInput()->GetLargestPossibleRegion().GetIndex();

    for(unsigned int i=0; i<TInputImage::ImageDimension; i++)
      {
      	inputSize[i]=inputLargSize[i];
      	inputIndex[i]=inputLargIndex[i];

      }//end for

    RequestedRegion.SetSize(inputSize);
    RequestedRegion.SetIndex(inputIndex);
    InputImagePointer input = const_cast< TInputImage * > ( this->GetInput() );
    input->SetRequestedRegion (RequestedRegion);
    }//end if

  itkDebugMacro("GenerateInputRequestedRegion End");
}//end GenerateInputRequestedRegion


/**
 * GenerateData Performs the diagonal image creation
 */
template <class TInputImage, class TOutputImage>
void
SegmentLiverFilter<TInputImage, TOutputImage>
::GenerateData( void )
{
  itkDebugMacro("GenerateData Start");

  auto inputImage = this->GetInput();
  auto outputImage = this->GetOutput();
  outputImage->SetBufferedRegion( outputImage->GetRequestedRegion() );
  outputImage->Allocate();

  //1. Threshold the input image
  itkDebugMacro("1. Threshold the input image");
  using ThresholdFilter = itk::BinaryThresholdImageFilter<TInputImage, TOutputImage>;
  auto thresholdFilter = ThresholdFilter::New();
  thresholdFilter->SetLowerThreshold(m_LowerThreshold);
  //thresholdFilter->SetUpperThreshold(m_UpperThreshold);
  thresholdFilter->SetInput(inputImage);
  thresholdFilter->SetInsideValue(1);


  itkDebugMacro("2. Fill any 3D holes in the thresholded image");
  //2. Fill any 3D holes in the thresholded image
  using HoleFilter = itk::GrayscaleFillholeImageFilter<TOutputImage, TOutputImage>;
  auto holeFilter = HoleFilter::New();
  holeFilter->SetInput(thresholdFilter->GetOutput());
  

  itkDebugMacro("3. Create a distance map of the entire image");
  //3. Create a distance map of the entire image
  using MaurerImage = itk::Image<float,3>;
  using MaurerFilter = itk::SignedMaurerDistanceMapImageFilter<TOutputImage, MaurerImage>;
  auto distanceMapFilter = MaurerFilter::New();
  distanceMapFilter->SetInput(holeFilter->GetOutput());
  if (m_RealSpacing == true)
  {  distanceMapFilter->UseImageSpacingOn();  }
  else
  {  distanceMapFilter->UseImageSpacingOff();  }
  distanceMapFilter->SquaredDistanceOff();
  distanceMapFilter->Update();	//It appears to making the values negative, so I changed to searching for minima instead of maxima


  itkDebugMacro("4. Manually remove from the distance map values that are part of the background in the thresholded image");
  //4. Manually remove from the distance map values that are part of the background in the thresholded image
  using OutputIteratorType = itk::ImageRegionIterator<TOutputImage>;
  using MaurerIteratorType = itk::ImageRegionIterator<MaurerImage>; 
  auto distanceMap = distanceMapFilter->GetOutput();
  MaurerIteratorType itDist(distanceMap, distanceMap->GetLargestPossibleRegion());
  OutputIteratorType itBin(holeFilter->GetOutput(), holeFilter->GetOutput()->GetLargestPossibleRegion());  
  itDist.GoToBegin();
  itBin.GoToBegin();
  while (!itDist.IsAtEnd())
  {
    if (itBin.Get() == 0)	{	itDist.Set(0);	}
    ++itDist;
    ++itBin;
  }//end while !itDist.IsAtEnd()
  


  itkDebugMacro("5. Find regional minima, the points furthest into the object part of the image");
  //5. Find regional minima, the points furthest into the object part of the image
  using RegionalMinFilterType = itk::RegionalMinimaImageFilter<MaurerImage, MaurerImage>;
  auto regionalMinFilter = RegionalMinFilterType::New();
  regionalMinFilter->SetInput(distanceMap);
  regionalMinFilter->SetForegroundValue(1);
  regionalMinFilter->SetBackgroundValue(0);
  regionalMinFilter->SetFlatIsMinima(false);
  regionalMinFilter->Update();


  itkDebugMacro("6. Find the edge of the sector where the liver may be");
  //6. Find the edge of the sector where the liver may be
  auto regionalMinima = regionalMinFilter->GetOutput();
  using IndexedMaurerIteratorType = itk::ImageRegionIteratorWithIndex<MaurerImage>;
  IndexedMaurerIteratorType itMin(regionalMinima, regionalMinima->GetLargestPossibleRegion());
  itMin.GoToBegin();
  itDist.GoToBegin();

  MaurerImage::IndexType locatorIndex;
  if (m_UseBrainCentroid)
  {
	if (m_SpacingDistance == true)
	{
	    regionalMinima->TransformPhysicalPointToIndex(m_BrainCentroid, locatorIndex);
	    locatorIndex[0] -= m_MinXFromBrain;
	    locatorIndex[2] -= m_MinZFromBrain;
	}//end if m_SpacingDistance == true
	else
	{
 	    PointType locatorPoint = m_BrainCentroid;
	    locatorPoint[0] -= m_MinXFromBrain;
	    locatorPoint[2] -= m_MinZFromBrain;
	    regionalMinima->TransformPhysicalPointToIndex(locatorPoint, locatorIndex);
	}//end else
  }//end if m_UseBrainCentroid
  else
  {
	if (m_SpacingDistance == true)
	{  regionalMinima->TransformPhysicalPointToIndex(m_SectorLocation, locatorIndex);  }
	else
	{
	    locatorIndex[0] = m_SectorLocation[0];
	    locatorIndex[1] = m_SectorLocation[1];
	    locatorIndex[2] = m_SectorLocation[2];
	}//end else

  }//end else


  itkDebugMacro("7. Find the lowest of the regional minima within the sector and set it as the location for the sphere.");
  //7. Find the lowest of the regional minima within the sector and set it as the location for the sphere.
  MaurerImage::IndexType lowestIndex;
  lowestIndex.Fill(0);
  MaurerImage::PointType lowestLocation;
  lowestLocation.Fill(0);

  float lowestValue = 0;
  while (!itMin.IsAtEnd())
  {
    auto minIndex = itMin.GetIndex();
    if ( (m_SetUseSearchROI && m_SearchROI.IsInside(minIndex)) ||
         (itMin.Get() != 0 && itMin.GetIndex()[0] < locatorIndex[0] && itMin.GetIndex()[2] < locatorIndex[2]))
    {
      if (itDist.Get() < lowestValue)
      {
	      if (inputImage->GetPixel(itMin.GetIndex()) < m_UpperThreshold)
  	    {
		      lowestValue = itDist.Get();
		      lowestIndex = itMin.GetIndex();
	      }//end if
      }//end if
    }//end if
    
    ++itMin;
    ++itDist;
  }//end while !itMin.IsAtEnd

  regionalMinFilter->GetOutput()->TransformIndexToPhysicalPoint(lowestIndex, lowestLocation);
  m_LiverCentroid[0] = lowestLocation[0];
  m_LiverCentroid[1] = lowestLocation[1];
  m_LiverCentroid[2] = lowestLocation[2];


  itkDebugMacro("8. Create the sphere by making a new image with just that pixel");
  //8. Create the sphere by making a new image with just that pixel
  auto sphereCenter = TOutputImage::New();
  typename TOutputImage::RegionType sphereCenterRegion;
  sphereCenterRegion = outputImage->GetLargestPossibleRegion();
  sphereCenter->SetRegions(sphereCenterRegion);
  sphereCenter->Allocate();
  if (m_RealSpacing)	{	sphereCenter->SetSpacing(outputImage->GetSpacing());	sphereCenter->SetOrigin(outputImage->GetOrigin());	}

  using IndexedOutputIteratorType = itk::ImageRegionIteratorWithIndex<TOutputImage>;
  IndexedOutputIteratorType itSC(sphereCenter, sphereCenter->GetLargestPossibleRegion());
  itSC.GoToBegin();
  while (!itSC.IsAtEnd())
  {
    if (itSC.GetIndex() == lowestIndex)
    {	itSC.Set(1);	}
    else
    {	itSC.Set(0);	}
    ++itSC;
  }//end while !itSC.IsAtEnd()


  itkDebugMacro("9. Make a distance map of the single-pixel image");
  //9. Make a distance map of the single-pixel image
  auto sphereDistance = MaurerFilter::New();
  sphereDistance->SetInput(sphereCenter);
  if (m_RealSpacing == true)
  {  sphereDistance->UseImageSpacingOn();  }
  else
  {  sphereDistance->UseImageSpacingOff();  }
  sphereDistance->UseImageSpacingOn();
  sphereDistance->SquaredDistanceOff();
  sphereDistance->Update();


  itkDebugMacro("10. Threshold that map with an upper limit of the radius");
  //10. Threshold that map with an upper limit of the radius
  m_Radius = lowestValue * -1;
  using ThresholdSphereFilter = itk::BinaryThresholdImageFilter<MaurerImage, TOutputImage>;
  auto thresholdSphereFilter = ThresholdSphereFilter::New();
  thresholdSphereFilter->SetUpperThreshold(lowestValue*(-1));//Should be- pixel value at the original point		Previously was- m_RadiusSphere
  thresholdSphereFilter->SetInput(sphereDistance->GetOutput());
  thresholdSphereFilter->SetInsideValue(1);
  //thresholdSphereFilter->Update();
  /*typedef itk::ImageFileWriter<TOutputImage> ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(thresholdSphereFilter->GetOutput());
  writer->SetFileName("liver_after_sphereFilter.vtk");
  writer->Update();*/

  itkDebugMacro("11. Erode the edge of the image by the erosion variable");
  //11. Erode the image edge
  using StructuringElementType = itk::BinaryBallStructuringElement<OutputImagePixelType, 3 > ;
  StructuringElementType SE;
  
  //pixels = distance/spacing
  int radius = std::ceil(m_Erosion);
  SE.SetRadius(radius);
  SE.CreateStructuringElement();

  using ErodeImageFilterType = itk::BinaryErodeImageFilter< TOutputImage, TOutputImage, StructuringElementType>;
  auto eroder = ErodeImageFilterType::New();
  eroder->SetKernel( SE );
  eroder->SetErodeValue( 1 );
  eroder->SetInput( thresholdSphereFilter->GetOutput() );
  eroder->Update(); 


  itkDebugMacro("12. Copy that image into the output image");
  //11. Copy that image into the output image
  auto tmp = eroder->GetOutput();
  tmp->SetOrigin(inputImage->GetOrigin());
  tmp->SetSpacing(inputImage->GetSpacing());


  using IteratorType = itk::ImageRegionIterator<TOutputImage> ;
  IteratorType it1(tmp, tmp->GetLargestPossibleRegion());
  IteratorType it2(outputImage, outputImage->GetLargestPossibleRegion());
  it1.GoToBegin();	it2.GoToBegin();
  while (!it1.IsAtEnd())
  {	it2.Set(it1.Get());	++it1;	++it2;	}
  
  itkDebugMacro("GenerateData End");
}//end GenerateData


/*
*GetBoundarySize
*	Returns a physical point that represents the size of the region that makes up the seed.  Cannot be run before data is generated.
*Input:		(member variables)
*Output:	OutputImageType::PointType return - A point representing the physical dimensions of the bounding box around the seed.
*/
template <class TInputImage, class TOutputImage>
typename TOutputImage::PointType SegmentLiverFilter<TInputImage, TOutputImage>::GetBoundarySize()
{
	using OutputIterator = itk::ImageRegionIteratorWithIndex<OutputImageType>;
	OutputImagePointer output = this->GetOutput();
        OutputIterator it(output, output->GetLargestPossibleRegion());
	it.GoToBegin();
	typename OutputImageType::IndexType lowIndex;
	lowIndex[0] = itk::NumericTraits<int>::max();	lowIndex[1] = itk::NumericTraits<int>::max();	lowIndex[2] = itk::NumericTraits<int>::max();
	typename OutputImageType::IndexType highIndex;
	highIndex[0] = 0;	highIndex[1] = 0;	highIndex[2] = 0;
	while (!it.IsAtEnd())
	{
		if (it.Get() == 1)
		{
			if (it.GetIndex()[0] < lowIndex[0])	lowIndex[0] = it.GetIndex()[0];
			if (it.GetIndex()[1] < lowIndex[1])	lowIndex[1] = it.GetIndex()[1];
			if (it.GetIndex()[2] < lowIndex[2])	lowIndex[2] = it.GetIndex()[2];

			if (it.GetIndex()[0] > highIndex[0])	highIndex[0] = it.GetIndex()[0];
			if (it.GetIndex()[1] > highIndex[1])	highIndex[1] = it.GetIndex()[1];
			if (it.GetIndex()[2] > highIndex[2])	highIndex[2] = it.GetIndex()[2];
		}//end if it.Get() == 1
		++it;
	}//end while !it.IsAtEnd()

	PointType boundarySize;
	PointType boundaryStart;
	PointType boundaryEnd;
	output->TransformIndexToPhysicalPoint(lowIndex, boundaryStart);
	output->TransformIndexToPhysicalPoint(highIndex, boundaryEnd);
	boundarySize[0] = boundaryEnd[0] - boundaryStart[0];
	boundarySize[1] = boundaryEnd[1] - boundaryStart[1];
	boundarySize[2] = boundaryEnd[2] - boundaryStart[2];
        if (highIndex[0] == 0)
	{
		boundarySize[0] = 0;
		boundarySize[1] = 0;
		boundarySize[2] = 0;
	}

	return boundarySize;
}//end GetBoundarySize


/*
*GetBoundarySize
*	Returns a physical point that represents the lowest coordinate of the box around the region that makes up the seed.  Cannot be run before data is generated.
*Input:		(member variables)
*Output:	OutputImageType::PointType return - A point representing the the starting location of the bounding box around the seed.
*/
template <class TInputImage, class TOutputImage>
typename TOutputImage::PointType SegmentLiverFilter<TInputImage, TOutputImage>::GetBoundaryStart()
{
	using OutputIterator = itk::ImageRegionIteratorWithIndex<OutputImageType> ;
	OutputImagePointer output = this->GetOutput();
        OutputIterator it(output, output->GetLargestPossibleRegion());
	it.GoToBegin();
	typename OutputImageType::IndexType lowIndex;
	lowIndex[0] = itk::NumericTraits<int>::max();	lowIndex[1] = itk::NumericTraits<int>::max();	lowIndex[2] = itk::NumericTraits<int>::max();
	//Find the single lowest index value in each dimension that has been segmented.
	while (!it.IsAtEnd())
	{
		if (it.Get() == 1)
		{
			if (it.GetIndex()[0] < lowIndex[0])	lowIndex[0] = it.GetIndex()[0];
			if (it.GetIndex()[1] < lowIndex[1])	lowIndex[1] = it.GetIndex()[1];
			if (it.GetIndex()[2] < lowIndex[2])	lowIndex[2] = it.GetIndex()[2];
		}//end if it.Get() == 1
		++it;
	}//end while !it.IsAtEnd()
	PointType boundaryStart;
	output->TransformIndexToPhysicalPoint(lowIndex, boundaryStart);
        if (lowIndex[0] == itk::NumericTraits<int>::max())
	{
		boundaryStart[0] = 0;
		boundaryStart[1] = 0;
		boundaryStart[2] = 0;
	}
	return boundaryStart;
}//end GetBoundaryStart





/* When printed, it will show the midpoint and all 3 vectors. */
template <class TInputImage, class TOutputImage>
void SegmentLiverFilter<TInputImage, TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Brain Centroid: (" << m_BrainCentroid[0] << ',' << m_BrainCentroid[1] << ',' << m_BrainCentroid[2] << ')' << std::endl;
  os << indent << "Min X from brain: " << m_MinXFromBrain << ((m_SpacingDistance) ? " units" : " voxels") << std::endl;
  os << indent << "Min Z from brain: " << m_MinZFromBrain << ((m_SpacingDistance) ? " units" : " voxels") << std::endl;
  os << indent << "Erosion value: " << m_Erosion << " mm" << std::endl;

}//end PrintSelf

} // end namespace itk


#endif
