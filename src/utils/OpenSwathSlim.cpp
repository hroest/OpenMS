// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathSlim
  : public TOPPOpenSwathBase 
{
public:

  TOPPOpenSwathSlim()
    : TOPPOpenSwathBase("OpenSwathSlim", "Complete workflow to run OpenSWATH (simplified version, use OpenSwathWorkflow for full version)", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,sqMass"));

    registerInputFile_("tr", "<file>", "", "Spectral library file ('TraML','tsv','pqp')");
    setValidFormats_("tr", ListUtils::create<String>("traML,tsv,pqp"));

    registerInputFile_("tr_irt", "<file>", "", "Retention time normalization / iRT library file ('TraML')", false);
    setValidFormats_("tr_irt", ListUtils::create<String>("traML,tsv,pqp"));

    registerInputFile_("tr_irt_nonlinear", "<file>", "", "Param estimation library - use these peptides to learn parameters from the data ('TraML')", false);
    setValidFormats_("tr_irt_nonlinear", ListUtils::create<String>("traML,tsv,pqp"));

    registerOutputFile_("out_peaks", "<file>", "", "OpenSwath peakgroup output file (mProphet compatible TSV or OSW file)", false);
    setValidFormats_("out_peaks", ListUtils::create<String>("tsv,osw,featureXML"));

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms output in mzML (chrom.mzML) or sqMass (SQLite format)", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML,sqMass"));

    registerStringOption_("rt_extraction_window", "<double>/auto", "auto", "Retention time extraction window", false, true);
    registerStringOption_("ion_mobility_window", "<double>/auto", "auto", "Extraction window in ion mobility dimension (in milliseconds). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side.", false);
    registerStringOption_("mz_extraction_window", "<double>/auto", "auto", "MZ time extraction window in ppm", false, true);

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,regression_delta_ppm,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));

    registerStringOption_("readOptions", "<name>", "workingInMemory", "Allows OpenSWATH to cache data to disk instead of keepign all data in memory (use cacheWorkingInMemory) in case you run out of memory. If you choose cacheWorkingInMemory, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("workingInMemory,cacheWorkingInMemory"));

    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("uis_threshold_sn",0);
      feature_finder_param.setValue("uis_threshold_peak_area",0);

      feature_finder_param.setValue("Scores:use_ms1_correlation", "true");
      feature_finder_param.setValue("Scores:use_ms1_fullscan", "true");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
      feature_finder_param.setValue("Scores:use_mi_score", "true");
      feature_finder_param.setValue("Scores:use_total_mi_score", "true");
      return feature_finder_param;
    }
    else if (name == "RTNormalization")
    {
      Param p;

      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing.");
      p.setValidStrings("alignmentMethod", ListUtils::create<String>("linear,interpolated,lowess,b_spline"));
      p.setValue("lowess:span", 2.0/3, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);

      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", ListUtils::create<String>("iter_residual,iter_jackknife,ransac,none"));

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", ListUtils::create<String>("true,false"));

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", ListUtils::create<String>("true,false"));

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (name == "Debugging")
    {
      Param p;
      p.setValue("irt_mzml", "", "Chromatogram mzML containing the iRT peptides");
      p.setValue("irt_trafo", "", "Transformation file for RT transform");
      return p;
    }
    else if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  ExitCodes main_(int, const char **) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    Param irt_detection_param = getSubsectionDefaults_("RTNormalization");
    Param feature_finder_param = getSubsectionDefaults_("Scoring");
    Param tsv_reader_param = getSubsectionDefaults_("Library");
    Param debug_params = getParam_().copy("Debugging:", true);

    FileTypes::Type tr_type = FileHandler::getType(tr_file);
    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type for '-tr' !");
      return PARSE_ERROR;
    }

    String output = getStringOption_("out_peaks");
    String out_tsv, out, out_osw;

    FileTypes::Type out_type = FileHandler::getType(output);
    if (out_type == FileTypes::FEATUREXML) out = output;
    else if (out_type == FileTypes::TSV) out = output;
    else if (out_type == FileTypes::OSW) out_osw = output;

    String irt_tr_file = getStringOption_("tr_irt");
    String nonlinear_irt_tr_file = getStringOption_("tr_irt_nonlinear");
    String out_chrom = getStringOption_("out_chrom");

    double mz_extraction_window = getStringOption_("mz_extraction_window").toDouble();
    double rt_extraction_window = getStringOption_("rt_extraction_window").toDouble();
    double im_extraction_window = getStringOption_("ion_mobility_window").toDouble();
    Size debug_level = (Size)getIntOption_("debug");
    double min_upper_edge_dist = 0.0;

    double min_rsq = 0.95;
    double min_coverage = 0.6;

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");
    String tmp = getStringOption_("tempDirectory");

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = true;
    if (readoptions == "cacheWorkingInMemory") readoptions = "cache";
    else if (readoptions == "workingInMemory") readoptions = "normal";

    if (!out_osw.empty() && tr_type != FileTypes::PQP)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OSW output files can only be generated in combination with PQP input files (-tr).");
    }

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = mz_extraction_window;
    cp.ppm                   = true;
    cp.rt_extraction_window  = rt_extraction_window;
    cp.im_extraction_window  = im_extraction_window;
    cp.extraction_function   = "tophat";
    cp.extra_rt_extract      = 0.0;

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range
    cp_irt.mz_extraction_window = 100;
    cp_irt.ppm                  = true;

    ChromExtractParams cp_ms1 = cp;

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(tr_type, tr_file, tsv_reader_param);
    LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, " <<
      transition_exp.getCompounds().size() << " compounds with " << transition_exp.getTransitions().size() << " transitions." << std::endl;

    if (tr_type == FileTypes::PQP)
    {
      remove(out_osw.c_str());
      if (!out_osw.empty())
      {
        std::ifstream  src(tr_file.c_str(), std::ios::binary);
        std::ofstream  dst(out_osw.c_str(), std::ios::binary);

        dst << src.rdbuf();
      }
    }

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    bool split_file = false;
    loadSwathFiles(file_list, exp_meta, swath_maps, split_file, tmp, readoptions, "",  false,
            min_upper_edge_dist, false, false, false);

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    String irt_trafo_out; // = debug_params.getValue("irt_trafo");
    String irt_mzml_out; // = debug_params.getValue("irt_mzml");
    TransformationDescription trafo_rtnorm;
    if (nonlinear_irt_tr_file.empty())
    {
      trafo_rtnorm = loadTrafoFile("", irt_tr_file, swath_maps,
                                   min_rsq, min_coverage, feature_finder_param,
                                   cp_irt, irt_detection_param, mz_correction_function, debug_level,
                                   false, load_into_memory, irt_trafo_out, irt_mzml_out);
    }
    else
    {
      ///////////////////////////////////
      // First perform a simple linear transform, then do a second, nonlinear one
      ///////////////////////////////////

      Param linear_irt = irt_detection_param;
      linear_irt.setValue("alignmentMethod", "linear");
      trafo_rtnorm = loadTrafoFile("", irt_tr_file, swath_maps,
                                   min_rsq, min_coverage, feature_finder_param,
                                   cp_irt, linear_irt, "none", debug_level,
                                   false, load_into_memory, irt_trafo_out, irt_mzml_out);

      cp_irt.rt_extraction_window = 900; // extract some substantial part of the RT range (should be covered by linear correction)
      cp_irt.rt_extraction_window = 600; // extract some substantial part of the RT range (should be covered by linear correction)

      ///////////////////////////////////
      // Get the secondary transformation (nonlinear)
      ///////////////////////////////////
      OpenSwath::LightTargetedExperiment transition_exp_nl;
      transition_exp_nl = loadTransitionList(FileHandler::getType(nonlinear_irt_tr_file), nonlinear_irt_tr_file, tsv_reader_param);

      std::vector< OpenMS::MSChromatogram > chromatograms;
      OpenSwathRetentionTimeNormalization wf;
      wf.setLogType(log_type_);
      wf.simpleExtractChromatograms(swath_maps, transition_exp_nl, chromatograms,
                                    trafo_rtnorm, cp_irt, false, load_into_memory);

      TransformationDescription im_trafo;
      trafo_rtnorm = wf.RTNormalization(transition_exp_nl, chromatograms, im_trafo, min_rsq,
                                        min_coverage, feature_finder_param, irt_detection_param,
                                        swath_maps, mz_correction_function, cp_irt.im_extraction_window,
                                        cp_irt.mz_extraction_window, cp_irt.ppm);

      TransformationDescription im_trafo_inv = im_trafo;
      im_trafo_inv.invert(); // theoretical -> experimental

      // We now modify the library as this is the easiest thing to do
      for (auto & p : transition_exp.getCompounds())
      {
        p.drift_time = im_trafo_inv.apply(p.drift_time);
      }
    }

    ///////////////////////////////////
    // Set up chromatogram output
    // Either use chrom.mzML or sqlite DB
    ///////////////////////////////////
    Interfaces::IMSDataConsumer * chromatogramConsumer;
    prepareChromOutput(&chromatogramConsumer, exp_meta, transition_exp, out_chrom);

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    FeatureMap out_featureFile;
    bool use_ms1_traces = true;

    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0], use_ms1_traces, false, false); // only active if filename not empty
    OpenSwathOSWWriter oswwriter(out_osw, file_list[0], use_ms1_traces, false, false); // only active if filename not empty

    {
      OpenSwathWorkflow wf(use_ms1_traces, true, -1);
      wf.setLogType(log_type_);
      wf.performExtraction(swath_maps, trafo_rtnorm, cp, cp_ms1, feature_finder_param, transition_exp,
          out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, 0, load_into_memory);
    }

    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    delete chromatogramConsumer;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathSlim tool;
  return tool.main(argc, argv);
}

/// @endcond
