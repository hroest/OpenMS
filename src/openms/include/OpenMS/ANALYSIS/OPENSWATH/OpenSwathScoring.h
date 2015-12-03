// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHSCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHSCORING_H

// data access
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace OpenMS
{

  /** @brief A structure to store which scores should be used by the Algorithm
   *
   * This can be used to turn on/off individual scores.
  */
  struct OPENMS_DLLAPI OpenSwath_Scores_Usage
  {
    // Which scores to use
    bool use_coelution_score_;
    bool use_shape_score_;
    bool use_rt_score_;
    bool use_library_score_;
    bool use_elution_model_score_;
    bool use_intensity_score_;
    bool use_total_xic_score_;
    bool use_nr_peaks_score_;
    bool use_sn_score_;
    bool use_dia_scores_;
    
    OpenSwath_Scores_Usage() :
      use_coelution_score_(true),
      use_shape_score_(true),
      use_rt_score_(true),
      use_library_score_(true),
      use_elution_model_score_(true),
      use_intensity_score_(true),
      use_total_xic_score_(true),
      use_nr_peaks_score_(true),
      use_sn_score_(true),
      use_dia_scores_(true)
    {}

    bool use_ms1_correlation;
    bool use_ms1_fullscan;
    bool use_uis_scores;
  };

  /** @brief A structure to hold the different scores computed by OpenSWATH
   *
   * This struct is used to store the individual OpenSWATH (sub-)scores. It
   * also allows to compute some preliminary quality score for a feature by
   * using a predefined combination of the individual scores determined using
   * LDA.
   *
  */
  struct OPENMS_DLLAPI OpenSwath_Scores
  {
    double elution_model_fit_score;
    double library_corr;
    double library_norm_manhattan;
    double library_rootmeansquare;
    double library_sangle;
    double norm_rt_score;
    double isotope_correlation;
    std::string ind_isotope_correlation;
    double isotope_overlap;
    std::string ind_isotope_overlap;
    double massdev_score;
    std::string ind_massdev_score;
    double xcorr_coelution_score;
    std::string ind_xcorr_coelution_score;
    double xcorr_shape_score;
    std::string ind_xcorr_shape_score;
    double yseries_score;
    double bseries_score;
    double log_sn_score;
    std::string ind_log_sn_score;
    int ind_num_transitions;
    std::string ind_transition_names;
    std::string ind_log_intensity;

    double weighted_coelution_score;
    double weighted_xcorr_shape;
    double weighted_massdev_score;
   
    double xcorr_ms1_coelution_score;
    double xcorr_ms1_shape_score;
    double ms1_ppm_score;
    double ms1_isotope_correlation;
    double ms1_isotope_overlap;

    double library_manhattan;
    double library_dotprod;
    double intensity;
    double total_xic;
    double nr_peaks;
    double sn_ratio;

    double rt_difference;
    double normalized_experimental_rt;
    double raw_rt_score;

    double dotprod_score_dia;
    double manhatt_score_dia;

    OpenSwath_Scores() :
      elution_model_fit_score(0),
      library_corr(0),
      library_norm_manhattan(0),
      library_rootmeansquare(0),
      library_sangle(0),
      norm_rt_score(0),
      isotope_correlation(0),
      ind_isotope_correlation(""),
      isotope_overlap(0),
      ind_isotope_overlap(""),
      massdev_score(0),
      ind_massdev_score(""),
      xcorr_coelution_score(0),
      ind_xcorr_coelution_score(""),
      xcorr_shape_score(0),
      ind_xcorr_shape_score(""),
      yseries_score(0),
      bseries_score(0),
      log_sn_score(0),
      ind_log_sn_score(""),
      ind_num_transitions(0),
      ind_transition_names(""),
      ind_log_intensity(""),
      weighted_coelution_score(0),
      weighted_xcorr_shape(0),
      weighted_massdev_score(0),
      xcorr_ms1_coelution_score(0),
      xcorr_ms1_shape_score(0),
      ms1_ppm_score(0),
      ms1_isotope_correlation(0),
      ms1_isotope_overlap(0),
      library_manhattan(0),
      library_dotprod(0),
      intensity(0),
      total_xic(0),
      nr_peaks(0),
      sn_ratio(0),
      dotprod_score_dia(0),
      manhatt_score_dia(0)
    {
    }

    double calculate_lda_prescore(OpenSwath_Scores scores, int mode = 0)
    {
      // NOTE this score is "better" if it is more negative!

      if (mode == 0)
      {

        // LDA average model on 100 2 x Crossvalidated runs (0.91 TPR/0.20 FDR)
        // legacy mode
        return scores.library_corr                     * -0.34664267 +
               scores.library_norm_manhattan           *  2.98700722 +
               scores.norm_rt_score                    *  7.05496384 +
               scores.xcorr_coelution_score            *  0.09445371 +
               scores.xcorr_shape_score                * -5.71823862 +
               scores.log_sn_score                     * -0.72989582 +
               scores.elution_model_fit_score          *  1.88443209;
      }
      else
      {

        // New LDA model computed using decoys on a set of 4723 decoy
        // and 4767 target transition groups. 
        // The following model contains the major important scores, it
        // recovered 1056 targets at 1% FDR (1559 at 10%) which is down from
        // about 1623 targets at 1% (2080 at 10% FDR) for a full model with all
        // scores (including SWATH-MS scores).
        //
        return scores.library_corr                     *   -0.119628004821 +
               scores.library_norm_manhattan           *   11.672611162623 +
               scores.norm_rt_score                    *   12.210070281257 +
               scores.xcorr_coelution_score            *   0.2469041031634 +
               scores.xcorr_shape_score                *   -3.973416923624 +
               scores.log_sn_score                     *   -0.768155423197;
      }
    }

    double calculate_swath_lda_prescore(OpenSwath_Scores scores, int mode = 0)
    {
      // NOTE this score is "better" if it is more negative!

      if (mode == 0)
      {

        // Swath - LDA average model on 100 2 x Crossvalidated runs (0.76 TPR/0.20 FDR) [without elution model]
        // legacy mode
        return scores.library_corr              * -0.19011762 +
               scores.library_norm_manhattan    *  2.47298914 +
               scores.norm_rt_score             *  5.63906731 +
               scores.isotope_correlation       * -0.62640133 +
               scores.isotope_overlap           *  0.36006925 +
               scores.massdev_score             *  0.08814003 +
               scores.xcorr_coelution_score     *  0.13978311 +
               scores.xcorr_shape_score         * -1.16475032 +
               scores.yseries_score             * -0.19267813 +
               scores.log_sn_score              * -0.61712054; 

      } 
      else 
      {

        // Swath - new LDA model computed using decoys on a set of 4723 decoy
        // and 4767 target transition groups. 
        // The following model contains the major important scores, it
        // recovered 1455 targets at 1% FDR (1955 at 10%) which is down from
        // about 1623 targets at 1% (2080 at 10% FDR) for a full model with all
        // scores.

        return scores.library_corr              *      0.274072972  +
               scores.library_norm_manhattan    *    15.4852263888  +
               scores.norm_rt_score             *    15.1274235495  +
               scores.isotope_correlation       *  -  2.6452497564  +
               scores.isotope_overlap           *     1.9356503434  +
               scores.massdev_score             *     0.0927833181  +
               scores.xcorr_coelution_score     *     0.2122083012  +
               scores.xcorr_shape_score         *  -  0.9086736117  +
               scores.yseries_score             *  -  0.0362475331  +
               scores.log_sn_score              *  -  0.3343873037 ;
      }
    }

  };

  /** @brief A class that calls the scoring routines
   *
   * Use this class to invoke the individual OpenSWATH scoring routines.
   * 
  */
  class OPENMS_DLLAPI OpenSwathScoring 
  {
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightTransition TransitionType;

    double rt_normalization_factor_;
    int add_up_spectra_;
    double spacing_for_spectra_resampling_;
    OpenSwath_Scores_Usage su_;

  public:

    /// Constructor
    OpenSwathScoring();

    /// Destructor
    ~OpenSwathScoring();

    /** @brief Initialize the scoring object
     *
     * Sets the parameters for the scoring.
     *
     * @param rt_normalization_factor Specifies the range of the normalized retention time space
     * @param add_up_spectra How many spectra to add up (default 1)
     * @param spacing_for_spectra_resampling Spacing factor for spectra addition
     * @param su Which scores to actually compute
     *
    */
    void initialize(double rt_normalization_factor,
      int add_up_spectra, double spacing_for_spectra_resampling,
      OpenSwath_Scores_Usage & su);

    /** @brief Score a single peakgroup in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. If more information is available, also consider using the
     * library based scoring and the full-spectrum based scoring.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids The list of native ids (giving a canonical ordering of the transitions)
     * @param normalized_library_intensity The weights to be used for each transition (e.g. normalized library intensities)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<std::string>& native_ids,
          const std::vector<double>& normalized_library_intensity,
          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
          OpenSwath_Scores & scores);

    /** @brief Score identification transitions against detection transitions of a single peakgroup 
     * in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. The scores are computed by scoring identification against detection
     * transitions.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids_identification The list of identification native ids (giving a canonical ordering of the transitions)
     * @param native_ids_detection The list of detection native ids (giving a canonical ordering of the transitions)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicIdScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<std::string>& native_ids_identification,
          const std::vector<std::string>& native_ids_detection,
          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
          OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature against a spectral library
     *
     * The spectral library is provided in a set of transition objects and a
     * peptide object. Both contain information about the expected elution time
     * on the chromatography and the relative intensity of the transitions.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param pep The peptide corresponding to the library transitions
     * @param normalized_feature_rt The retention time of the feature in normalized space
     * @param scores The object to store the result
     *
    */
    void calculateLibraryScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<TransitionType> & transitions,
          const PeptideType& pep,
          const double normalized_feature_rt,
          OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_map The SWATH-MS (DIA map) from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param ms1_map The corresponding MS1 (precursor ion map) from which the precursor spectra can be retrieved (optional, may be NULL)
     * @param diascoring DIA Scoring object to use for scoring
     * @param pep The peptide corresponding to the library transitions
     * @param scores The object to store the result
     *
    */
    void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, 
        const std::vector<TransitionType> & transitions,
        OpenSwath::SpectrumAccessPtr swath_map,
        OpenSwath::SpectrumAccessPtr ms1_map,
        OpenMS::DIAScoring & diascoring,
        const PeptideType& pep,
        OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_map The SWATH-MS (DIA map) from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param diascoring DIA Scoring object to use for scoring
     * @param scores The object to store the result
     *
    */
    void calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
        const TransitionType & transition,
        OpenSwath::SpectrumAccessPtr swath_map,
        OpenMS::DIAScoring & diascoring,
        OpenSwath_Scores & scores);

    /** @brief Computing the normalized library intensities from the transition objects
     *
     * The intensities are normalized such that the sum to one.
     *
     * @param transitions The library transition to score the feature against
     * @param normalized_library_intensity The resulting normalized library intensities
     *
    */
    void getNormalized_library_intensities_(const std::vector<TransitionType> & transitions,
        std::vector<double>& normalized_library_intensity);

    /** @brief Returns an averaged spectrum
     *
     * This function will sum up (add) the intensities of multiple spectra
     * around the given retention time and return an "averaged" spectrum which
     * may contain less noise.
     *
     * @param swath_map The map containing the spectra
     * @param RT The target retention time
     * @param nr_spectra_to_add How many spectra to add up
     *
    */
    OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, 
        double RT, int nr_spectra_to_add);

  };
}

#endif
