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
 
 #include "itkImageFileWriter.h"

#include "itkPluginUtilities.h"

#include "PETLiverUptakeMeasurementCLP.h"

#include "itkSegmentLiverFilter.h"
#include "itkSegmentBrainFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include <iostream>

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
typedef float InputPixelType;
typedef unsigned char OutputPixelType;
typedef itk::Image<InputPixelType,  3> InputImageType;
typedef itk::Image<OutputPixelType, 3> OutputImageType;
typedef InputImageType::RegionType ROIType;

int getLiverRegionWithROI( int argc, char * argv[], InputImageType::Pointer petScan, OutputImageType::Pointer& liverRegion, ROIType roi )
{
    PARSE_ARGS;
    
    typedef itk::SegmentLiverFilter< InputImageType, OutputImageType > SegmentLiverFilterType;
    SegmentLiverFilterType::Pointer segmentLiver = SegmentLiverFilterType::New();
    segmentLiver->SetInput(petScan);
    segmentLiver->SetLowerThreshold(lowerThreshold);
    segmentLiver->SetUpperThreshold(upperThreshold);
    segmentLiver->SetErosion(erosion);
    segmentLiver->UseRealSpacing();
    segmentLiver->UseSpacedDistance();
    
    if (roi.GetSize()[0]>0) // user provided valid ROI for liver reference region search
      segmentLiver->SetSearchROI(roi);
    else // specify ROI for liver reference region search based on brain location
    {
        // identify brain
        const double threshold = 3.0;
        const double minimumVolume = 1000.0;
        typedef itk::Image<long, 3> BrainOutputImageType;
        typedef itk::SegmentBrainFilter< InputImageType, BrainOutputImageType > SegmentBrainFilterType;
        SegmentBrainFilterType::Pointer segmentBrain = SegmentBrainFilterType::New();
        segmentBrain->SetInput(petScan);
        segmentBrain->SetLowerThreshold(threshold);
        segmentBrain->SetMinimumVolume(minimumVolume);
        segmentBrain->Update();
        SegmentBrainFilterType::PointType brainCentroid = segmentBrain->GetCentroid();
        
        // specify ROI based on brain center
        const double minX = 20.0;
        const double minZ = 200.0;
        segmentLiver->SetBrainCentroid(brainCentroid);
        segmentLiver->SetMinXFromBrain(minX);
        segmentLiver->SetMinZFromBrain(minZ);
    }
    
    segmentLiver->Update();
    liverRegion = segmentLiver->GetOutput();

    return EXIT_SUCCESS;
}
    
ROIType convertToROI(std::vector<float> region, InputImageType::Pointer img)
{
    ROIType roi = img->GetLargestPossibleRegion();
    roi.SetSize(0,0); roi.SetSize(1,0); roi.SetSize(2,0); // default size if no user specified ROI was provided
    if (region.size()>=6) // search ROI provided
    {
        ROIType userROI;
        InputImageType::PointType pointA;
        pointA[0] = -region[0]-region[3];
        pointA[1] = -region[1]-region[4];
        pointA[2] = region[2]-region[5];
        InputImageType::PointType pointB;
        pointB[0] = -region[0]+region[3];
        pointB[1] = -region[1]+region[4];
        pointB[2] = region[2]+region[5];
        InputImageType::IndexType idxA;
        img->TransformPhysicalPointToIndex(pointA,idxA);
        InputImageType::IndexType idxB;
        img->TransformPhysicalPointToIndex(pointB,idxB);
        InputImageType::SizeType ROISize;
        ROISize[0] = abs(idxA[0]-idxB[0])+1;
        ROISize[1] = abs(idxA[1]-idxB[1])+1;
        ROISize[2] = abs(idxA[2]-idxB[2])+1;
        InputImageType::IndexType ROIStart;
        ROIStart[0] = std::min(idxA[0], idxB[0]);
        ROIStart[1] = std::min(idxA[1], idxB[1]);
        ROIStart[2] = std::min(idxA[2], idxB[2]);
        userROI.SetIndex(ROIStart);
        userROI.SetSize(ROISize);
        roi = img->GetLargestPossibleRegion();
        roi.Crop(userROI);
    }
    return roi;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  try
    {
        typedef float InputPixelType;
        typedef unsigned char OutputPixelType;
        typedef itk::Image<InputPixelType,  3> InputImageType;
        typedef itk::Image<OutputPixelType, 3> OutputImageType;
        typedef InputImageType::RegionType ROIType;
        
        // read pet scan
        typedef itk::ImageFileReader<InputImageType>  ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( inputVolume.c_str() );
        reader->Update();
        InputImageType::Pointer petScan = reader->GetOutput();
        
        // identify liver reference region
        ROIType liverSearchROI = convertToROI(region, petScan);
        OutputImageType::Pointer liverRegion = OutputImageType::New();
        getLiverRegionWithROI(argc, argv, petScan, liverRegion, liverSearchROI);
        
        // write segmentation result, if requested
        if (outputVolume.size()>0)
        {
          typedef itk::ImageFileWriter<OutputImageType> WriterType;
          WriterType::Pointer writer = WriterType::New();
          writer->SetFileName( outputVolume.c_str() );
          writer->SetInput( liverRegion );
          writer->SetUseCompression(1);
          writer->Update();
        }
        
        // obtain measurements
        typedef itk::LabelStatisticsImageFilter<InputImageType, OutputImageType> StatsFilterType;
        StatsFilterType::Pointer stats = StatsFilterType::New();
        stats->SetInput(petScan);
        stats->SetLabelInput(liverRegion);
        stats->UseHistogramsOn(); // required to calculate median value
        stats->SetHistogramParameters(100000,0,1000); // we need fine grained binning to obtain a reasonably accurate median value
        stats->Update();
        StatsFilterType::LabelPixelType label = 1;
        size_t numVoxels = stats->GetCount(label);
        double meanValue = stats->GetMean(label);
        double stdValue = sqrt(stats->GetVariance(label));
        double medianValue = stats->GetMedian(label);
        
        // write measurements
        std::ofstream writeFile;
        writeFile.open( returnParameterFile.c_str() );
        if (numVoxels==0)
          writeFile << "Mean_s = --" << std::endl;
        else 
          writeFile << "Mean_s = " << meanValue << std::endl;
        if (numVoxels==0)
          writeFile << "Std_s = --" << std::endl;
        else 
          writeFile << "Std_s = " << stdValue << std::endl;
        if (numVoxels==0)
          writeFile << "Median_s = --" << std::endl;
        else 
          writeFile << "Median_s = " << medianValue << std::endl;
        writeFile.close();
        
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
