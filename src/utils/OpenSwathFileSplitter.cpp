// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>



// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <assert.h>

#define OPENSWATH_WORKFLOW_DEBUG

using namespace OpenMS;

static bool SortPairDoubleByFirst(const std::pair<double,double> & left, const std::pair<double,double> & right)
{
  return left.first < right.first;
}

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


class swathfilenew : public SwathFile
{

  public:

    /// Loads a Swath run from a single mzML file
    std::vector<OpenSwath::SwathMap> NewloadMzML(String file, String tmp,
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions = "normal")
    {
      std::cout << "Loading mzML file " << file << " using readoptions " << readoptions << std::endl;
      String tmp_fname = "openswath_tmpfile";

      startProgress(0, 1, "Loading metadata file " + file);
      boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata = populateMetaData_(file);
      exp_meta = experiment_metadata;

      // First pass through the file -> get the meta data
      std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
      std::vector<int> swath_counter;
      int nr_ms1_spectra;
      std::vector<OpenSwath::SwathMap> known_window_boundaries;
      countScansInSwath_(experiment_metadata->getSpectra(), swath_counter, nr_ms1_spectra, known_window_boundaries);
      std::cout << "Determined there to be " << swath_counter.size() <<
        " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
      endProgress();

      FullSwathFileConsumer* dataConsumer;
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      startProgress(0, 1, "Loading data file " + file);
      if (readoptions == "normal")
      {
        dataConsumer = new RegularSwathFileConsumer(known_window_boundaries);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "cache")
      {
        dataConsumer = new CachedSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "split")
      {
        dataConsumer = new MzMLSwathFileConsumer(known_window_boundaries, tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "Unknown or unsupported option " + readoptions);
      }
      LOG_DEBUG << "Finished parsing Swath file " << std::endl; 
      std::vector<OpenSwath::SwathMap> swath_maps;
      dataConsumer->retrieveSwathMaps(swath_maps);
      delete dataConsumer;

      endProgress();
      return swath_maps;
    }

};


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathFileSplitter
  : public TOPPBase
{
public:

  TOPPOpenSwathFileSplitter()
    : TOPPBase("OpenSwathWorkflow", "Split SWATH files", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML')", false);
    setValidFormats_("tr_irt", ListUtils::create<String>("traML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, true);
    registerFlag_("sort_swath_maps", "Sort of input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerFlag_("use_ms1_traces", "Extract the precursor ion trace(s) and use for scoring", true);
    registerFlag_("enable_uis_scoring", "Enable additional scoring of identification assays", true);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

    registerStringOption_("out_tsv", "<file>", "", "TSV output file (mProphet compatible)", false);

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms (chrom.mzML) output", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false, true);
    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window that by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    registerDoubleOption_("mz_extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerFlag_("ppm", "m/z extraction_window is in ppm");

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false, true);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory"));

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));

    // TODO terminal slash !
    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 0, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 500-1000)", false, true);
    setMinInt_("batchSize", 0);

    // registerSubsection_("Scoring", "Scoring parameters section");

    // registerSubsection_("outlierDetection", "Parameters for the outlierDetection for iRT petides. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
  }

  Param getSubsectionDefaults_(const String& name) const
  {}

  void loadSwathFiles(StringList& file_list, bool split_file, String tmp, String readoptions,
    boost::shared_ptr<ExperimentalSettings > & exp_meta,
    std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    swathfilenew swath_file;
    swath_file.setLogType(log_type_);

    if (split_file || file_list.size() > 1)
    {
      // TODO cannot use data reduction here any more ...
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
    }
    else
    {
      FileTypes::Type in_file_type = FileTypes::nameToType(file_list[0]);
      if (in_file_type == FileTypes::MZML || file_list[0].suffix(4).toLower() == "mzml"
        || file_list[0].suffix(7).toLower() == "mzml.gz"  )
      {
        // swath_maps = swath_file.NewloadMzML(file_list[0], tmp, exp_meta, readoptions);
        swath_maps = swath_file.NewloadMzML(file_list[0], tmp, exp_meta, "split");
      }
      else if (in_file_type == FileTypes::MZXML || file_list[0].suffix(5).toLower() == "mzxml"
        || file_list[0].suffix(8).toLower() == "mzxml.gz"  )
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Input file needs to have ending mzML or mzXML");
      }
    }
  }

  ExitCodes main_(int, const char **)
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    //String tr_file = getStringOption_("tr");

    // Param irt_detection_param = getParam_().copy("outlierDetection:", true);

    //tr_file input file type

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");

    /*
    String irt_tr_file = getStringOption_("tr_irt");
    String trafo_in = getStringOption_("rt_norm");
    */

    String out_chrom = getStringOption_("out_chrom");
    bool ppm = getFlag_("ppm");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getFlag_("use_ms1_traces");
    bool enable_uis_scoring = getFlag_("enable_uis_scoring");
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    double mz_extraction_window = getDoubleOption_("mz_extraction_window");
    double rt_extraction_window = getDoubleOption_("rt_extraction_window");
    double extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");
    String swath_windows_file = getStringOption_("swath_windows_file");
    int batchSize = (int)getIntOption_("batchSize");
    Size debug_level = (Size)getIntOption_("debug");

    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");
    String tmp = getStringOption_("tempDirectory");

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }

    /*
    if (trafo_in.empty() && irt_tr_file.empty())
    {
      std::cout << "Since neither rt_norm nor tr_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }
    if ( (out.empty() && out_tsv.empty()) || (!out.empty() && !out_tsv.empty()) )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Either out_features or out_tsv needs to be set (but not both)");
    }
    */

    // Check swath window input
    if (!swath_windows_file.empty())
    {
      LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    /*
    OpenSwathWorkflow::ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = mz_extraction_window;
    cp.ppm                   = ppm;
    cp.rt_extraction_window  = rt_extraction_window,
    cp.extraction_function   = extraction_function;
    cp.extra_rt_extract      = extra_rt_extract;

    OpenSwathWorkflow::ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range
    */

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////

    // (i) Load files
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    loadSwathFiles(file_list, split_file, tmp, readoptions, exp_meta, swath_maps);

    // (ii) Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
    {
      SwathWindowLoader::annotateSwathMapsFromFile(swath_windows_file, swath_maps, sort_swath_maps);
    }

    for (Size i = 0; i < swath_maps.size(); i++)
    {
      LOG_DEBUG << "Found swath map " << i << " with lower " << swath_maps[i].lower
        << " and upper " << swath_maps[i].upper << " and " << swath_maps[i].sptr->getNrSpectra()
        << " spectra." << std::endl;
    }

    // (iii) Sanity check: there should be no overlap between the windows:
    std::vector<std::pair<double, double> > sw_windows;
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      if (!swath_maps[i].ms1)
      {
        sw_windows.push_back(std::make_pair(swath_maps[i].lower, swath_maps[i].upper));
      }
    }
    std::sort(sw_windows.begin(), sw_windows.end(), SortPairDoubleByFirst);

    for (Size i = 1; i < sw_windows.size(); i++)
    {
      double lower_map_end = sw_windows[i-1].second - min_upper_edge_dist;
      double upper_map_start = sw_windows[i].first;
      LOG_DEBUG << "Extraction will go up to " << lower_map_end << " and continue at " << upper_map_start << std::endl;

      if (upper_map_start - lower_map_end > 0.01)
      {
        LOG_INFO << "Extraction will have a gap between " << lower_map_end << " and " << upper_map_start << std::endl;
        if (!force)
        {
          LOG_INFO << "Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

      if (lower_map_end - upper_map_start > 0.01)
      {
        LOG_INFO << "Extraction will overlap between " << lower_map_end << " and " << upper_map_start << std::endl;
        LOG_INFO << "This will lead to multiple extraction of the transitions in this region which should be avoided." << std::endl;
        LOG_INFO << "Please fix this by providing an appropriate extraction file with -swath_windows_file" << std::endl;
        if (!force)
        {
          LOG_INFO << "Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

    }

    // DONE !!! 

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathFileSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
