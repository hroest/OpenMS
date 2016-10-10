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

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow


// basic file operations
#include <iostream>
#include <fstream>

namespace OpenMS
{

  /// Constructor
  OpenSwathScoring::OpenSwathScoring() :
    rt_normalization_factor_(1.0),
    add_up_spectra_(1),
    spacing_for_spectra_resampling_(0.005)
  {
  }

  /// Destructor
  OpenSwathScoring::~OpenSwathScoring()
  {
  }

  void OpenSwathScoring::initialize(double rt_normalization_factor,
    int add_up_spectra, double spacing_for_spectra_resampling,
    OpenSwath_Scores_Usage & su)
  {
    this->rt_normalization_factor_ = rt_normalization_factor;
    this->add_up_spectra_ = add_up_spectra;
    this->spacing_for_spectra_resampling_ = spacing_for_spectra_resampling;
    this->su_ = su;
  }


  void xcorr_compute(std::vector<std::vector<double> >& sonar_profiles, 
                     double& xcorr_coelution_score, double& xcorr_shape_score)
  {
    /// Cross Correlation array
    typedef std::map<int, double> XCorrArrayType;
    /// Cross Correlation matrix
    typedef std::vector<std::vector<XCorrArrayType> > XCorrMatrixType;

    std::vector<double> intensityi, intensityj;
    XCorrMatrixType xcorr_matrix;
    xcorr_matrix.resize(sonar_profiles.size());
    for (std::size_t i = 0; i < sonar_profiles.size(); i++)
    {
      xcorr_matrix[i].resize(sonar_profiles.size());
      for (std::size_t j = i; j < sonar_profiles.size(); j++)
      {
        // compute normalized cross correlation
        xcorr_matrix[i][j] = OpenSwath::Scoring::normalizedCrossCorrelation(
                sonar_profiles[i], sonar_profiles[j], boost::numeric_cast<int>(sonar_profiles[i].size()), 1);
      }
    }

    // coelution (lag score)
    std::vector<int> deltas;
    for (std::size_t i = 0; i < xcorr_matrix.size(); i++)
    {
      for (std::size_t  j = i; j < xcorr_matrix.size(); j++)
      {
        // first is the lag value, should be an int
        deltas.push_back(std::abs(OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr_matrix[i][j])->first));
      }
    }

    {
      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(deltas.begin(), deltas.end(), msc);
      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();
      xcorr_coelution_score = deltas_mean + deltas_stdv;
    }



    // shape score (intensity)
    std::vector<double> intensities;
    for (std::size_t i = 0; i < xcorr_matrix.size(); i++)
    {
      for (std::size_t j = i; j < xcorr_matrix.size(); j++)
      {
        // second is the Y value (intensity)
        intensities.push_back(OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr_matrix[i][j])->second);
      }
    }
    {
      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(intensities.begin(), intensities.end(), msc);
      xcorr_shape_score = msc.mean();
    }
  }

  void sonar_scores(OpenSwath::IMRMFeature* imrmfeature,
                                            const std::vector<OpenSwath::LightTransition> & transitions,
                                            std::vector<OpenSwath::SwathMap> swath_maps,
                                            OpenSwath::SpectrumAccessPtr ms1_map,
                                            OpenMS::DIAScoring & diascoring, 
                                            const OpenSwath::LightCompound& compound, OpenSwath_Scores & scores)
  {
    double precursor_mz = transitions[0].getPrecursorMZ();

    std::ofstream debug_file;
    debug_file.open("debug_sonar_profiles.tsv",  std::fstream::in | std::fstream::out | std::fstream::app);

    String native_id = 0;
    if (transitions.size() > 0)
    {
      native_id = transitions[0].getNativeID();
    }
    debug_file << native_id << "\t" << imrmfeature->getRT() << "\tcentr";
    for (Size it = 0; it < swath_maps.size(); it++)
    {
      debug_file << "\t" << (swath_maps[it].lower + swath_maps[it].upper) / 2.0;
    }
    debug_file << "\n";


    // idea 1: check the elution profile of each SONAR scan ...
    for (int kk = 0; kk < imrmfeature->getNativeIDs().size(); kk++)
    {
      std::vector<double> rt;
      imrmfeature->getFeature(imrmfeature->getNativeIDs()[kk])->getRT(rt);
      // std::cout << " for feathre  " << imrmfeature->getNativeIDs()[kk] << " st: " << rt[0] << " to " << rt.back() << std::endl;
    }

    std::cout << " doing RT " << imrmfeature->getRT() << " using maps: " ;
    for (int i  = 0; i < swath_maps.size() ; i++)
    {
      std::cout << (swath_maps[i].lower + swath_maps[i].upper) / 2 << " " ;
    }
    std::cout << std::endl;



    // idea 2: check the SONAR profile (e.g. in the dimension of) of the best scan (RT apex)
    double RT = imrmfeature->getRT();

    //double dia_extract_window_ = 0.1;
    double dia_extract_window_ = 1.0;
    bool dia_centroided_ = false;

    // Aggregate sonar profiles
    std::vector<std::vector<double> > sonar_profiles;

    std::vector<double> sn_score;
    std::vector<double> diff_score;
    std::vector<double> trend_score;
    std::vector<double> rsq_score;
    std::vector<double> mz_median_score;
    std::vector<double> mz_stdev_score;
    for (Size k = 0; k < transitions.size(); k++)
    {
      String native_id = transitions[k].getNativeID();
      // double rel_intensity = intensities[native_id];
      // If no charge is given, we assume it to be 1
      int putative_fragment_charge = 1;
      if (transitions[k].fragment_charge > 0)
      {
        putative_fragment_charge = transitions[k].fragment_charge;
      }

      std::cout << " transition " << native_id << " at " << RT << " will analyze with " << swath_maps.size() << " maps" << std::endl;

      // Gather profiles 
      std::vector<double> sonar_profile;
      std::vector<double> sonar_mz_profile;
      std::vector<bool> signal_exp;
      for (int swath_idx = 0; swath_idx < swath_maps.size(); swath_idx++)
      {
        OpenSwath::SpectrumAccessPtr swath_map = swath_maps[swath_idx].sptr;
        // std::cout << "  swath_idx " << swath_idx << (swath_maps[swath_idx].lower + swath_maps[swath_idx].upper) / 2.0 << std::endl;
        
        bool expect_signal = false;
        if (swath_maps[swath_idx].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[swath_idx].lower && precursor_mz < swath_maps[swath_idx].upper) 
        {
          // std::cout << "   expect signal... " << std::endl;
          expect_signal = true;
        }
        else
        {
          // std::cout << " expect no signal... " << std::endl;
        }

        // find closest 
        std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
        // std::cout << " try to find RT " << RT << " found nr indices " << indices.size() << std::endl;
        if (indices.empty() )  {continue;}
        int closest_idx = boost::numeric_cast<int>(indices[0]);
        if (indices[0] != 0 &&
            std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
            std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
        {
          closest_idx--;
        }
        OpenSwath::SpectrumPtr spectrum_ = swath_map->getSpectrumById(closest_idx);
        double left = transitions[k].getProductMZ() - dia_extract_window_ / 2.0;
        double right = transitions[k].getProductMZ() + dia_extract_window_ / 2.0; 
        double mz, intensity;
        integrateWindow(spectrum_, left, right, mz, intensity, dia_centroided_);

        sonar_profile.push_back(intensity);
        sonar_mz_profile.push_back(mz);
        signal_exp.push_back(expect_signal);

        // std::cout << "   " << swath_idx << " integrated data " << mz << " int : " << intensity << std::endl;
      }
      sonar_profiles.push_back(sonar_profile);

      debug_file << native_id << "\t" << imrmfeature->getRT() << "\tint";
      for (Size it = 0; it < sonar_profile.size(); it++)
      {
        debug_file << "\t" << sonar_profile[it];
      }
      debug_file << "\n";
      debug_file << native_id << "\t" << imrmfeature->getRT() << "\tmz";
      for (Size it = 0; it < sonar_mz_profile.size(); it++)
      {
        debug_file << "\t" << sonar_mz_profile[it];
      }
      debug_file << "\n";

      // Analyze profiles 
      std::vector<double> sonar_profile_pos;
      std::vector<double> sonar_mz_profile_pos;
      std::vector<double> sonar_profile_neg;
      std::vector<double> sonar_mz_profile_neg;
      for (Size it = 0; it < sonar_profile.size(); it++)
      {
        if (signal_exp[it])
        {
          sonar_profile_pos.push_back(sonar_profile[it]);
          sonar_mz_profile_pos.push_back(sonar_mz_profile[it]);
        }
        else
        {
          sonar_profile_neg.push_back(sonar_profile[it]);
          sonar_mz_profile_neg.push_back(sonar_mz_profile[it]);
        }
      }

      // try to find diff between first and last 
      double sonar_trend = 1.0;
      if (sonar_profile_pos.size() > 1)
      {
        double int_end = sonar_profile_pos[sonar_profile_pos.size()-1] + sonar_profile_pos[sonar_profile_pos.size()-2];
        double int_start = sonar_profile_pos[0] + sonar_profile_pos[1];
        if (int_end > 0.0)
        {
          sonar_trend = int_start / int_end;
        }
        else
        {
          sonar_trend = 0.0;
        }
      }

      // try to find R^2 of a linear regression (optimally, there is no trend)
      std::vector<double> xvals; 
      for (int pr_idx = 0; pr_idx < sonar_profile_pos.size(); pr_idx++) {xvals.push_back(pr_idx);}
      Math::LinearRegression lr;
      lr.computeRegression(0.95, xvals.begin(), xvals.end(), sonar_profile_pos.begin());
      double rsq = lr.getRSquared();

      // try to find largest diff
      double sonar_largediff = 0.0;
      for (int pr_idx = 0; pr_idx < sonar_profile_pos.size()-1; pr_idx++)
      {
        double diff = std::fabs(sonar_profile_pos[pr_idx] - sonar_profile_pos[pr_idx+1]);
        if (diff > sonar_largediff) {sonar_largediff = diff;}
      }

      double sonar_sn = 1.0;
      double pos_med = 1.0;
      double neg_med = 1.0;
        std::cout << " compute from profile sizes " << 
      sonar_profile_pos.size() << " and " <<  sonar_profile_neg.size() << std::endl;

      // from here on, its not sorted any more !!
      if (!sonar_profile_pos.empty() && !sonar_profile_neg.empty())
      {
        pos_med = Math::median(sonar_profile_pos.begin(), sonar_profile_pos.end()); 
        neg_med = Math::median(sonar_profile_neg.begin(), sonar_profile_neg.end()); 

        // compute the relative difference between the medians (or if the
        // medians are zero, compute the difference to the max element)
        if (neg_med > 0.0)
        {
          sonar_sn = pos_med / neg_med;
        }
        else if (*std::max_element(sonar_profile_neg.begin(), sonar_profile_neg.end()) > 0.0)
        {
          sonar_sn = pos_med / *std::max_element(sonar_profile_neg.begin(), sonar_profile_neg.end());
        }

      }

      double median_mz = 0.0;
      double mz_stdev = -1.0;
      if (!sonar_mz_profile_pos.empty())
      {
        median_mz = Math::medianFast(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end()); 

        double sum = std::accumulate(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end(), 0.0);
        double mean = sum / sonar_mz_profile_pos.size();

        double sq_sum = std::inner_product(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end(), sonar_mz_profile_pos.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / sonar_mz_profile_pos.size() - mean * mean);

        mz_stdev = stdev;
      }

      std::cout << " computed SN: " << sonar_sn  <<  "(from " << pos_med << " and neg " << neg_med <<  ")"
        << " large diff: "  << sonar_largediff << " trend " << sonar_trend << std::endl;
      sn_score.push_back(sonar_sn);
      diff_score.push_back(sonar_largediff / pos_med);
      trend_score.push_back(sonar_trend);
      rsq_score.push_back(rsq);

      mz_median_score.push_back(median_mz);
      mz_stdev_score.push_back(mz_stdev);
    }
        
    double xcorr_coelution_score, xcorr_shape_score;
    xcorr_compute(sonar_profiles, xcorr_coelution_score, xcorr_shape_score);

    double sn_av = std::accumulate(sn_score.begin(), sn_score.end(), 0.0) / sn_score.size();
    double diff_av = std::accumulate(diff_score.begin(), diff_score.end(), 0.0) / diff_score.size();
    double trend_av = std::accumulate(trend_score.begin(), trend_score.end(), 0.0) / trend_score.size();
    double rsq_av = std::accumulate(rsq_score.begin(), rsq_score.end(), 0.0) / rsq_score.size();

    double mz_median = std::accumulate(mz_median_score.begin(), mz_median_score.end(), 0.0) / mz_median_score.size();
    double mz_stdev = std::accumulate(mz_stdev_score.begin(), mz_stdev_score.end(), 0.0) / mz_stdev_score.size();

    scores.sonar_sn = sn_av;
    scores.sonar_diff = diff_av;
    scores.sonar_trend = trend_av;
    scores.sonar_rsq = rsq_av;
    scores.sonar_lag = xcorr_coelution_score;
    scores.sonar_shape = xcorr_shape_score;

    debug_file.close();
  }

  void OpenSwathScoring::calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature,
                                            const std::vector<TransitionType> & transitions,
                                            std::vector<OpenSwath::SwathMap> swath_maps,
                                            OpenSwath::SpectrumAccessPtr ms1_map,
                                            OpenMS::DIAScoring & diascoring, 
                                            const CompoundType& compound, OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(transitions.size() > 0, "There needs to be at least one transition.");
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    std::vector<OpenSwath::SwathMap> used_swath_maps;
    if (swath_maps.size() > 1 || transitions.empty())
    {
      // std::cout << " dia scores1 , sonar1 " << std::endl;

      double precursor_mz = transitions[0].getPrecursorMZ();
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        if (swath_maps[i].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[i].lower && precursor_mz < swath_maps[i].upper) 
        {
          used_swath_maps.push_back(swath_maps[i]);
          // std::cout << " will use map  sonar " << swath_maps[i].lower << " -  " << swath_maps[i].upper << std::endl;
        }
      }

      // TODO compute some scores ...
      //sonar_scores(imrmfeature, transitions, used_swath_maps, ms1_map, diascoring, compound, scores);
      sonar_scores(imrmfeature, transitions, swath_maps, ms1_map, diascoring, compound, scores);
    }
    else
    {
      used_swath_maps = swath_maps;
    }

    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum_ = getAddedSpectra_(used_swath_maps, imrmfeature->getRT(), add_up_spectra_);
    OpenSwath::SpectrumPtr* spectrum = &spectrum_;

    // Mass deviation score
    diascoring.dia_massdiff_score(transitions, (*spectrum), normalized_library_intensity,
        scores.massdev_score, scores.weighted_massdev_score);

    // DIA dotproduct and manhattan score based on library intensity
    diascoring.score_with_isotopes((*spectrum), transitions, scores.dotprod_score_dia, scores.manhatt_score_dia);

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    // Currently this is computed for an averagine model of a peptide so its
    // not optimal for metabolites - but better than nothing, given that for
    // most fragments we dont really know their composition
    diascoring.dia_isotope_scores(transitions, (*spectrum), imrmfeature, scores.isotope_correlation, scores.isotope_overlap);

    // Peptide-specific scores
    if (compound.isPeptide())
    {
      // Presence of b/y series score
      OpenMS::AASequence aas;
      int by_charge_state = 1; // for which charge states should we check b/y series
      OpenSwathDataAccessHelper::convertPeptideToAASequence(compound, aas);
      diascoring.dia_by_ion_score((*spectrum), aas, by_charge_state, scores.bseries_score, scores.yseries_score);
    }

    // FEATURE we should not punish so much when one transition is missing!
    scores.massdev_score = scores.massdev_score / transitions.size();

    // Compute precursor-level scores:
    // - compute mass difference in ppm
    // - compute isotopic pattern score
    if (ms1_map && ms1_map->getNrSpectra() > 0) 
    {
      double precursor_mz = transitions[0].precursor_mz;
      OpenSwath::SpectrumPtr ms1_spectrum = getAddedSpectra_(ms1_map, imrmfeature->getRT(), add_up_spectra_);
      diascoring.dia_ms1_massdiff_score(precursor_mz, ms1_spectrum, scores.ms1_ppm_score);

      int precursor_charge = 1;
      if (compound.getChargeState() != 0) {precursor_charge = compound.getChargeState();}

      if (compound.isPeptide())
      {
        diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum,
                                          precursor_charge, scores.ms1_isotope_correlation,
                                          scores.ms1_isotope_overlap);
      }
      else
      {
        diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum,
                                          precursor_charge, scores.ms1_isotope_correlation,
                                          scores.ms1_isotope_overlap, compound.sum_formula);
      }
    }
  }

  void OpenSwathScoring::calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature, 
                                              const TransitionType & transition,
                                              std::vector<OpenSwath::SwathMap> swath_maps,
                                              OpenMS::DIAScoring & diascoring,
                                              OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    std::vector<OpenSwath::SwathMap> used_swath_maps;
    if (swath_maps.size() > 1)
    {
      // std::cout << " dia scores1 , sonar " << std::endl;


      double precursor_mz = transition.getPrecursorMZ();
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        if (swath_maps[i].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[i].lower && precursor_mz < swath_maps[i].upper) 
        {
          used_swath_maps.push_back(swath_maps[i]);
          // std::cout << " will use map  sonar " << swath_maps[i].lower << " -  " << swath_maps[i].upper << std::endl;
        }
      }
    }
    else
    {
      used_swath_maps = swath_maps;
    }

    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum_ = getAddedSpectra_(used_swath_maps, imrmfeature->getRT(), add_up_spectra_);
    OpenSwath::SpectrumPtr* spectrum = &spectrum_;

    // If no charge is given, we assume it to be 1
    int putative_product_charge = 1;
    if (transition.getProductChargeState() > 0)
    {
      putative_product_charge = transition.getProductChargeState();
    }

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    diascoring.dia_ms1_isotope_scores(transition.getProductMZ(), (*spectrum), putative_product_charge, scores.isotope_correlation, scores.isotope_overlap);
    // Mass deviation score
    diascoring.dia_ms1_massdiff_score(transition.getProductMZ(), (*spectrum), scores.massdev_score);
  }

  void OpenSwathScoring::calculateChromatographicScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids,
        const std::vector<double>& normalized_library_intensity,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & scores)
  {
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrMatrix(imrmfeature, native_ids);

    // XCorr score (coelution)
    if (su_.use_coelution_score_)
    {
      scores.xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
      scores.weighted_coelution_score = mrmscore_.calcXcorrCoelutionScore_weighted(normalized_library_intensity);
    }

    // XCorr score (shape)
    // mean over the intensities at the max of the crosscorrelation
    // FEATURE : weigh by the intensity as done by mQuest
    // FEATURE : normalize with the intensity at the peak group apex?
    if (su_.use_shape_score_)
    {
      scores.xcorr_shape_score = mrmscore_.calcXcorrShape_score();
      scores.weighted_xcorr_shape = mrmscore_.calcXcorrShape_score_weighted(normalized_library_intensity);
    }

    // check that the MS1 feature is present and that the MS1 correlation should be calculated
    if (imrmfeature->getPrecursorIDs().size() > 0 && su_.use_ms1_correlation)
    {
      mrmscore_.initializeMS1XCorr(imrmfeature, native_ids, "Precursor_i0"); // perform cross-correlation on monoisotopic precursor
      scores.xcorr_ms1_coelution_score = mrmscore_.calcMS1XcorrCoelutionScore();
      scores.xcorr_ms1_shape_score = mrmscore_.calcMS1XcorrShape_score();
    }

    if (su_.use_nr_peaks_score_) 
    { 
      scores.nr_peaks = boost::numeric_cast<int>(imrmfeature->size());
    }

    // Signal to noise scoring
    if (su_.use_sn_score_)
    {
      scores.sn_ratio = mrmscore_.calcSNScore(imrmfeature, signal_noise_estimators);
      // everything below S/N 1 can be set to zero (and the log safely applied)
      if (scores.sn_ratio < 1) { scores.log_sn_score = 0; }
      else { scores.log_sn_score = std::log(scores.sn_ratio); }
    }
  }

  void OpenSwathScoring::calculateChromatographicIdScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids_identification,
        const std::vector<std::string>& native_ids_detection,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & idscores)
  { 
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrIdMatrix(imrmfeature, native_ids_identification, native_ids_detection);

    if (su_.use_coelution_score_)
    {
      idscores.ind_xcorr_coelution_score = mrmscore_.calcIndXcorrIdCoelutionScore();
    }

    if (su_.use_shape_score_)
    {
      idscores.ind_xcorr_shape_score = mrmscore_.calcIndXcorrIdShape_score();
    }

    // Signal to noise scoring
    if (su_.use_sn_score_)
    {
      idscores.ind_log_sn_score = mrmscore_.calcIndSNScore(imrmfeature, signal_noise_estimators);
    }
  }

  void OpenSwathScoring::calculateLibraryScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<TransitionType> & transitions,
        const CompoundType& pep,
        const double normalized_feature_rt,
        OpenSwath_Scores & scores)
  {
    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    std::vector<std::string> native_ids;
    OpenSwath::MRMScoring mrmscore_;
    for (Size i = 0; i < transitions.size(); i++) {native_ids.push_back(transitions[i].getNativeID());}

    if (su_.use_library_score_)
    {
      mrmscore_.calcLibraryScore(imrmfeature, transitions, 
          scores.library_corr, scores.library_norm_manhattan, scores.library_manhattan, 
          scores.library_dotprod, scores.library_sangle, scores.library_rootmeansquare);
    }

    // Retention time score
    if (su_.use_rt_score_)
    {
      // rt score is delta iRT
      double normalized_experimental_rt = normalized_feature_rt;
      double rt_score = mrmscore_.calcRTScore(pep, normalized_experimental_rt);

      scores.normalized_experimental_rt = normalized_experimental_rt;
      scores.raw_rt_score = rt_score;
      scores.norm_rt_score = rt_score / rt_normalization_factor_;
    }
  }

  void OpenSwathScoring::getNormalized_library_intensities_(const std::vector<TransitionType> & transitions,
      std::vector<double>& normalized_library_intensity)
  {
    normalized_library_intensity.clear();
    for (Size i = 0; i < transitions.size(); i++) 
    {
      normalized_library_intensity.push_back(transitions[i].getLibraryIntensity());
    }
    for (Size i = 0; i < normalized_library_intensity.size(); i++) 
    { 
      // the library intensity should never be below zero
      if (normalized_library_intensity[i] < 0.0) { normalized_library_intensity[i] = 0.0; } 
    } 
    OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0], boost::numeric_cast<int>(normalized_library_intensity.size()));
  }

  OpenSwath::SpectrumPtr OpenSwathScoring::getAddedSpectra_(std::vector<OpenSwath::SwathMap> swath_maps,
                                                            double RT, int nr_spectra_to_add)
  {
    if (swath_maps.size() == 1) 
    {
      return getAddedSpectra_(swath_maps[0].sptr, RT, nr_spectra_to_add);
    }
    else
    {
      std::vector<OpenSwath::SpectrumPtr> all_spectra;
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        OpenSwath::SpectrumPtr spec = getAddedSpectra_(swath_maps[i].sptr, RT, nr_spectra_to_add);
        all_spectra.push_back(spec);
      }
      OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      return spectrum_;
    }
  }

  OpenSwath::SpectrumPtr OpenSwathScoring::getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, 
                                                            double RT, int nr_spectra_to_add)
  {
    std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
    if (indices.empty() ) 
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
    }
    int closest_idx = boost::numeric_cast<int>(indices[0]);
    if (indices[0] != 0 &&
        std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
        std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
    {
      closest_idx--;
    }

    if (nr_spectra_to_add == 1)
    {
      OpenSwath::SpectrumPtr spectrum_ = swath_map->getSpectrumById(closest_idx);
      return spectrum_;
    }
    else
    {
      std::vector<OpenSwath::SpectrumPtr> all_spectra;
      // always add the spectrum 0, then add those right and left
      all_spectra.push_back(swath_map->getSpectrumById(closest_idx));
      for (int i = 1; i <= nr_spectra_to_add / 2; i++) // cast to int is intended!
      {
        if (closest_idx - i >= 0) 
        {
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx - i));
        }
        if (closest_idx + i < (int)swath_map->getNrSpectra()) 
        {
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx + i));
        }
      }
      OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      return spectrum_;
    }
  }

}

