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
\file 	itkSegmentBrainFilter.txx
\brief	Implementation of itkSegmentBrainFilter.h.
*/
/*
*itkSegmentBrainFilter.txx
*Implementation of itkSegmentBrainFilter.h.
*
*Programmer: Markus Van Tol
*Date: 8/18/11
*
*/

#ifndef __itkSegmentBrainFilter_txx
#define __itkSegmentBrainFilter_txx

#include <cmath>
#include <string>
#include <iostream>
#include "itkSegmentBrainFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkContinuousIndex.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

namespace itk
{
	
template <class TInputImage, class TOutputImage>
void
SegmentBrainFilter<TInputImage, TOutputImage>
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
	if (!input)
		std::cout << "!input" << std::endl;
	if (!output)
		std::cout << "!output" << std::endl;
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
SegmentBrainFilter<TInputImage, TOutputImage>
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
SegmentBrainFilter<TInputImage, TOutputImage>
::GenerateData( void )
{
  using OutputPixelType = typename TOutputImage::PixelType;

  auto inputImage = this->GetInput();
  auto outputImage = this->GetOutput();
  outputImage->SetBufferedRegion( outputImage->GetRequestedRegion() );
  outputImage->Allocate();

  using ThresholdFilter = itk::BinaryThresholdImageFilter<TInputImage, TOutputImage>;
  auto thresholdFilter = ThresholdFilter::New();
  thresholdFilter->SetLowerThreshold(m_LowerThreshold);
  thresholdFilter->SetInput(inputImage);
  thresholdFilter->SetInsideValue(1);
  thresholdFilter->SetOutsideValue(0);
  //image thresholded

  //Need to remove parts below a certain volume
  using TOutputImage2D = itk::Image<OutputPixelType, 2>;

  using HoleFilter = itk::GrayscaleFillholeImageFilter<TOutputImage2D, TOutputImage2D>;
  auto holeFilter = HoleFilter::New();

  using SbSFilterType = itk::SliceBySliceImageFilter<TOutputImage, TOutputImage, HoleFilter, HoleFilter, TOutputImage2D, TOutputImage2D>;
  auto sbsFilter = SbSFilterType::New();

  sbsFilter->SetFilter(holeFilter);
  sbsFilter->SetInput(thresholdFilter->GetOutput());

  using ComponentFilter = itk::ConnectedComponentImageFilter<TOutputImage, TOutputImage>;
  auto componentFilter = ComponentFilter::New();
  componentFilter->SetInput(sbsFilter->GetOutput());
  componentFilter->SetBackgroundValue( 0 );

  componentFilter->Update();
  auto filteredImage = componentFilter->GetOutput();

  using IteratorType = itk::ImageRegionIterator<TOutputImage> ;
  IteratorType it1(filteredImage, filteredImage->GetLargestPossibleRegion());
  IteratorType it2(outputImage, outputImage->GetLargestPossibleRegion());
  it1.GoToBegin();	it2.GoToBegin();
  while (!it1.IsAtEnd())
  {	it2.Set(it1.Get());	++it1;	++it2;	}

  FindCentroid(outputImage);

}//end GenerateData

/*
*FindCentroid
*Locates the centroid for the brain segmentation.
*Also locates the upper and lower points for the boundary box around the brain.
*/
template <class TInputImage, class TOutputImage>
void SegmentBrainFilter<TInputImage, TOutputImage>::
FindCentroid(OutputImagePointer outputImage)
{	
	//go to highest slice and find pixel value for the brain
	//
	//while slices contain that pixel value, keep going down
	//on each slice, get numVoxels and centroid point
	
	using IteratorType = itk::ImageRegionIteratorWithIndex<TOutputImage> ;
	IteratorType it1(outputImage, outputImage->GetLargestPossibleRegion());
	IteratorType it2(outputImage, outputImage->GetLargestPossibleRegion());

	PointType centroid;
	centroid[0] = 0;
	centroid[1] = 0;
	centroid[2] = 0;

	PointType highBoundary;
	highBoundary.Fill(0.0);
	PointType lowBoundary;
	lowBoundary.Fill(0.0);

	typename InputImageType::IndexType highest;
	highest[0] = 0;	highest[1] = 0;	highest[2] = 0;
	typename InputImageType::IndexType lowest;
	lowest[0] = 4294967295;	lowest[1] = 4294967295;	lowest[2] = 4294967295;
	int brainLabel = 0;
	it1.GoToBegin();
	it2.GoToReverseBegin();
	if (it1.GetIndex()[2] > it2.GetIndex()[2])
	{
		itkDebugMacro("normal it");
		//std::cout<<"normal it"<<std::endl;
		//int highestZ = it1.GetIndex()[2];
		while (/*it1.GetIndex()[2] == highestZ &&*/ brainLabel == 0 && it1.IsAtEnd() == false)
		{
			brainLabel = it1.Get();
			if (brainLabel != 0 && ! MinimumVolume(brainLabel, outputImage) )
			{	brainLabel = 0;	}
			++it1;
		}
		if (brainLabel != 0)
		{
			//std::cout<<"found label "<<brainLabel<<" at "<<it1.GetIndex()<<std::endl;
			//it1.GoToBegin();
			int currentZ = it1.GetIndex()[2];
			bool prev_had_brain = true;
			bool current_has_brain = false;
			unsigned long voxelCount = 0;

			while (it1.GetIndex()[2] == currentZ || prev_had_brain == true)
			{
				if (it1.GetIndex()[2] != currentZ)
				{
					prev_had_brain = current_has_brain;
					current_has_brain = false;	
					currentZ = it1.GetIndex()[2];
				}

				if (it1.Get() == brainLabel)
				{
					PointType tempPoint;
					typename InputImageType::IndexType tempIndex = it1.GetIndex();
					outputImage->TransformIndexToPhysicalPoint(tempIndex, tempPoint);
					voxelCount++;
					centroid[0]+=tempPoint[0];
					centroid[1]+=tempPoint[1];
					centroid[2]+=tempPoint[2];
					current_has_brain = true;

					for (int i = 0; i < 3; i++)
					{
						if (tempIndex[i] > highest[i])
						{	highBoundary[i] = tempPoint[i];	highest[i] = tempIndex[i];	}
						if (tempIndex[i] < lowest[i])
						{	lowBoundary[i] = tempPoint[i];	lowest[i] = tempIndex[i];	}
					}//end for
				}//end if
				++it1;
			}//end while
			//std::cout<<"voxelCount="<<voxelCount<<std::endl;
			centroid[0] /= (double) voxelCount;
			centroid[1] /= (double) voxelCount;
			centroid[2] /= (double) voxelCount;
		}
	}//end if
	else
	{
		itkDebugMacro("reverse it");
		//std::cout<<"reverse it"<<std::endl;
		//int highestZ = it2.GetIndex()[2];
		while (/*it2.GetIndex()[2] == highestZ &&*/ brainLabel == 0 && it2.IsAtReverseEnd() == false)
		{
			brainLabel = it2.Get();
			if (brainLabel != 0 && ! MinimumVolume(brainLabel, outputImage) )
			{	brainLabel = 0;	}
			--it2;
		}
		if (brainLabel != 0)
		{
			//std::cout<<"found label "<<brainLabel<<" at "<<it2.GetIndex()<<std::endl;
			//it2.GoToReverseBegin();
			int currentZ = it2.GetIndex()[2];
			bool prev_had_brain = true;
			bool current_has_brain = false;
			unsigned long voxelCount = 0;

			while (it2.GetIndex()[2] == currentZ || prev_had_brain == true)
			{
				if (it2.GetIndex()[2] != currentZ)
				{
					prev_had_brain = current_has_brain;
					current_has_brain = false;	
					currentZ = it2.GetIndex()[2];
				}//end if

				if (it2.Get() == brainLabel)
				{
					PointType tempPoint;
					typename InputImageType::IndexType tempIndex = it2.GetIndex();
					outputImage->TransformIndexToPhysicalPoint(tempIndex, tempPoint);
					voxelCount++;
					centroid[0]+=tempPoint[0];
					centroid[1]+=tempPoint[1];
					centroid[2]+=tempPoint[2];
					current_has_brain = true;

					for (int i = 0; i < 3; i++)
					{
						if (tempIndex[i] > highest[i])
						{	highBoundary[i] = tempPoint[i];	highest[i] = tempIndex[i];	}
						if (tempIndex[i] < lowest[i])
						{	lowBoundary[i] = tempPoint[i];	lowest[i] = tempIndex[i];		}
					}//end for

				}//end if
				--it2;
			}//end while

			//std::stringstream ss;
			//ss << "voxelCount == " << voxelCount;
			//std::cout<<"voxelCount="<<voxelCount<<std::endl;
			//itkDebugMacro(ss.str().c_str());

			centroid[0] /= (double) voxelCount;
			centroid[1] /= (double) voxelCount;
			centroid[2] /= (double) voxelCount;
		
		}//end if
		
	}//end else

	it1.GoToBegin();
	while (!it1.IsAtEnd())
	{
		if (it1.Get() != brainLabel || brainLabel == 0)
		{	it1.Set(0);	}
		else
		{	it1.Set(1);	}
		++it1;
	}//end while

	if (brainLabel != 0)
	{
		m_Centroid = centroid;
		m_BoundarySize[0] = highBoundary[0] - lowBoundary[0];
		m_BoundarySize[1] = highBoundary[1] - lowBoundary[1];
		m_BoundarySize[2] = highBoundary[2] - lowBoundary[2];
		m_BoundaryStart = lowBoundary;
	}
	else
	{
		std::cout<<"error,";
		centroid[0] = 0;
		centroid[1] = 0;
		centroid[2] = 0;
		m_Centroid = centroid;
		m_BoundarySize[0] = 0;	m_BoundarySize[1] = 0;	m_BoundarySize[2] = 0;
		m_BoundaryStart[0] = 0;	m_BoundaryStart[1] = 0;	m_BoundaryStart[2] = 0;
	}

}//end FindCentroid

template <class TInputImage, class TOutputImage>
bool SegmentBrainFilter<TInputImage, TOutputImage>::
MinimumVolume(OutputImagePixelType label, OutputImagePointer image)
{
	float spacePerVoxel = (image->GetSpacing()[0]) * (image->GetSpacing()[1]) * (image->GetSpacing()[2]);
	float voxelCount = 0;
	using IteratorType = itk::ImageRegionIteratorWithIndex<TOutputImage>;
	IteratorType it(image, image->GetLargestPossibleRegion());
	it.GoToBegin();
	while (!it.IsAtEnd())
	{
		if (it.Get() == label)
		{	voxelCount++;	}
		++it;
	}
	//std::cout<<"Label "<<label<<" volume = "<<voxelCount*spacePerVoxel<<std::endl;
	return (voxelCount * spacePerVoxel >= m_MinimumVolume);
	}

/* When printed, it will show the midpoint and all 3 vectors. */
template <class TInputImage, class TOutputImage>
void SegmentBrainFilter<TInputImage, TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf(os,indent);

	os << indent << "Centeroid: (" << m_Centroid[0] << ',' << m_Centroid[1] << ',' << m_Centroid[2] << ')' << std::endl;
	os << indent << "Bounding Box start: (" << m_BoundaryStart[0] << ',' << m_BoundaryStart[1] << ',' << m_BoundaryStart[2] << ')' << std::endl;
	os << indent << "Bounding Box size: (" << m_BoundaryStart[0] << ',' << m_BoundaryStart[1] << ',' << m_BoundaryStart[2] << ')' << std::endl;
}//end PrintSelf

} // end namespace itk

#endif
