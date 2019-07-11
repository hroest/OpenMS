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

#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>

// scoring
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

// #define DEBUG_IMSCORING

namespace OpenMS
{

  class MSDriftSpectrum 
  {
    std::vector<double> mz_;
    std::vector<double> intens_;
    std::vector<double> im_;

    bool sorted_by_mz = false;
    bool sorted_by_im = false;

		struct DoubleVectorComparator
		{
				const std::vector<double> & value_vector;

				explicit DoubleVectorComparator(const std::vector<double> & val_vec):
						value_vector(val_vec) {}

				bool operator()(int i1, int i2)
				{
						return value_vector[i1] < value_vector[i2];
				}
		};

    void init()
    {
      sorted_by_mz = std::adjacent_find(mz_.begin(), mz_.end(), std::greater<double>()) == mz_.end();
      sorted_by_im = std::adjacent_find(im_.begin(), im_.end(), std::greater<double>()) == im_.end();
    }

  public:

    MSDriftSpectrum () {};

    MSDriftSpectrum (const std::vector<double>& mz, const std::vector<double>& intens, const std::vector<double>& im) :
      mz_(mz), intens_(intens), im_(im)
    {
      init();
    }

    MSDriftSpectrum (std::vector<double>::const_iterator& mz_start,
                     std::vector<double>::const_iterator& mz_end,
                     std::vector<double>::const_iterator& intens_start,
                     std::vector<double>::const_iterator& im_start) :
      mz_(mz_start, mz_end), 
      intens_(intens_start, std::next(intens_start, std::distance(mz_start, mz_end))), 
      im_(im_start, std::next(im_start, std::distance(mz_start, mz_end)))
    {
      init();
    }

    MSDriftSpectrum (std::vector<double>::const_iterator& mz_start,
                     std::vector<double>::const_iterator& mz_end,
                     std::vector<double>::const_iterator& intens_start,
                     std::vector<double>::const_iterator& intens_end,
                     std::vector<double>::const_iterator& im_start,
                     std::vector<double>::const_iterator& im_end) :
      mz_(mz_start, mz_end), 
      intens_(intens_start, intens_end),
      im_(im_start, im_end)
    {
      init();
    }

    void projectMZAxis(double eps = 1e-6)
    {
      sortByMZ();

      std::vector<double> mz, intens;
      auto mz_it = mz_.begin();
      auto intens_it = intens_.begin();
      // auto im_it = im_.begin();
      mz.push_back(*mz_it);
      intens.push_back(*intens_it);
      ++mz_it; ++intens_it;

      auto mz_prev = mz_it;
      while (mz_it != mz_.end())
      {
        if ( fabs(*mz_it - *mz_prev) < eps)
        {
          intens.back() += *intens_it;
        }
        else
        {
          intens.push_back(*intens_it);
          mz.push_back(*mz_it);
        }
        ++mz_it; ++intens_it;
      }
    }

    void sortByMZ()
    {
      if (sorted_by_mz) return;

      // sort by m/z
    	std::sort(im_.begin(), im_.end(), DoubleVectorComparator(mz_));
    	std::sort(intens_.begin(), intens_.end(), DoubleVectorComparator(mz_));
    	std::sort(mz_.begin(), mz_.end() );

      sorted_by_im = false;
      sorted_by_mz = true;
    }

    void sortByIM()
    {
      if (sorted_by_im) return;

      // sort by ion mobility
    	std::sort(mz_.begin(), mz_.end(), DoubleVectorComparator(im_));
    	std::sort(intens_.begin(), intens_.end(), DoubleVectorComparator(im_));
    	std::sort(im_.begin(), im_.end() );

      sorted_by_im = true;
      sorted_by_mz = false;
    }

    void projectIMMZAxis(double eps = 1e-6)
    {
      sortByIM();

      std::vector<double> im, intens;
      // auto mz_it = mz_.begin();
      auto intens_it = intens_.begin();
      auto im_it = im_.begin();
      im.push_back(*im_it);
      intens.push_back(*intens_it);
      ++im_it; ++intens_it;

      auto im_prev = im_it;
      while (im_it != im_.end())
      {
        if ( fabs(*im_it - *im_prev) < eps)
        {
          intens.back() += *intens_it;
        }
        else
        {
          intens.push_back(*intens_it);
          im.push_back(*im_it);
        }
        ++im_it; ++intens_it;
      }
    }

    double getIMSpacing_sub(std::vector<double>& tmp)
    {
      double min_diff = fabs(tmp[tmp.size()-1] - tmp[0]);

      for (Size k = 1; k < tmp.size(); k++) 
      {
        double diff = fabs(tmp[k] - tmp[k-1]);
        if (diff > 1e-5 && diff < min_diff) min_diff = diff;
      }
      return min_diff;
    }

    static std::vector<double> getIMValues(const std::vector<double>& im, double eps = 1e-5)
    {
      OPENMS_PRECONDITION(im.size() >= 2, "Needs at least 2 entries")

      std::vector<double> tmp(im);
      std::sort(tmp.begin(), tmp.end());

      // min_diff = getIMSpacing_sub(tmp);

      std::vector<double> im_values;
      im_values.push_back( tmp[0] );
      for (Size k = 1; k < tmp.size(); k++) 
      {
        double diff = fabs(tmp[k] - tmp[k-1]);
        if (diff > eps)
        {
          im_values.push_back( tmp[k] );
        }
      }
      return im_values;
    }

  };
}

namespace OpenMS
{

  typedef std::pair<double, double> MobilityPeak;
  typedef std::vector< MobilityPeak > IMProfile;

  std::vector<double> computeGrid(const std::vector< IMProfile >& im_profiles, double eps)
  {
    // Extract all ion mobility values across all transitions and produce a
    // grid of all permitted ion mobility values
    std::vector<double> im_grid;
    {
      std::vector< double > mobilityValues;
      for (const auto & im_profile : im_profiles) 
      {
        for (const auto & k : im_profile) mobilityValues.push_back(k.first);
      }

      // sort all extracted values
      std::sort(mobilityValues.begin(), mobilityValues.end());

      // Reduce mobility values to grid (consider equal if closer than eps)
      // 
      // In some cases there are not enough datapoints available (one of the
      // transitions has no datapoints)
      if (!mobilityValues.empty())
      {
        im_grid.push_back( mobilityValues[0] );
        for (Size k = 1; k < mobilityValues.size(); k++) 
        {
          double diff = fabs(mobilityValues[k] - mobilityValues[k-1]);
          if (diff > eps)
          {
            im_grid.push_back( mobilityValues[k] );
          }
        }
      }
    }
    return im_grid;
  }

  /*
   @brief Extracts ion mobility values projected onto a grid

   For a given ion mobility profile and a grid, compute an ion mobilogram
   across the grid for each ion mobility data point. Returns two data arrays
   for the ion mobilogram: intensity (y) and ion mobility (x). Zero values are
   inserted if no data point was found for a given grid value.

   @param profile The ion mobility data
   @param im_grid The grid to be used
   @param al_int_values The intensity vector (y)
   @param al_im_values The ion mobility vector (x)
   @param max_peak_idx The grid position of the maximum

  */
  void alignToGrid(const IMProfile& profile,
               const std::vector<double>& im_grid,
               std::vector< double >& al_int_values,
               std::vector< double >& al_im_values,
               Size & max_peak_idx)
  {
    auto pr_it = profile.begin();
    max_peak_idx = 0;
    double max_int = 0;
    for (Size k = 0; k < im_grid.size(); k++)
    {
      // In each iteration, the IM value of pr_it should be equal to or
      // larger than the master container. If it is equal, we add the current
      // data point, if it is larger we add zero and advance the counter k.
      if (pr_it != profile.end() && fabs(pr_it->first - im_grid[k] ) < 1e-4 ) 
      {
        al_int_values.push_back(pr_it->second);
        al_im_values.push_back(pr_it->first);
        ++pr_it;
      }
      else
      {
        al_int_values.push_back(0.0);
        al_im_values.push_back( im_grid[k] );
      }
      // OPENMS_LOG_DEBUG << "grid position " << im_grid[k] << " profile position " << pr_it->first << std::endl;

      // check that we did not advance past 
      if (pr_it != profile.end() && (im_grid[k] - pr_it->first) > 1e-3)
      {
        std::cout << " This should never happen, pr_it has advanced past the master container: " << im_grid[k]  << "  / " <<  pr_it->first  << std::endl;
        throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      // collect maxima
      if (pr_it != profile.end() && pr_it->second > max_int)
      {
        max_int = pr_it->second;
        max_peak_idx = k;
      }
    }
  }

  /**
    @brief Integrate intensity in an ion mobility spectrum from start to end

    This function will integrate the intensity in a spectrum between mz_start
    and mz_end, returning the total intensity and an intensity-weighted drift
    time value.

    This function also returns the full ion mobility profile in "res".

    @note If there is no signal, mz will be set to -1 and intensity to 0
  */
  void integrateDriftSpectrum(OpenSwath::SpectrumPtr spectrum, 
                              double mz_start,
                              double mz_end,
                              double & im,
                              double & intensity,
                              IMProfile& res, 
                              double drift_start,
                              double drift_end)
  {
    OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");

    // rounding multiplier for the ion mobility value
    // TODO: how to improve this -- will work up to 42949.67296
    double IM_IDX_MULT = 10e5;

    std::map< int, double> im_chrom;
    {
      // get the weighted average for noncentroided data.
      // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
      typedef std::vector<double>::const_iterator itType;

      itType mz_arr_end = spectrum->getMZArray()->data.end();
      itType int_it = spectrum->getIntensityArray()->data.begin();
      itType im_it = spectrum->getDriftTimeArray()->data.begin();

      // this assumes that the spectra are sorted!
      itType mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(),
        spectrum->getMZArray()->data.end(), mz_start);
      itType mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

      // also advance intensity and ion mobility iterator now
      std::iterator_traits< itType >::difference_type iterator_pos = std::distance((itType)spectrum->getMZArray()->data.begin(), mz_it);
      std::advance(int_it, iterator_pos);
      std::advance(im_it, iterator_pos);

      // Iterate from mz start to end, only storing ion mobility values that are in the range
      for (; mz_it != mz_it_end; ++mz_it, ++int_it, ++im_it)
      {
        if ( *im_it >= drift_start && *im_it <= drift_end)
        {
          // std::cout << "IM " << *im_it << " mz " << *mz_it << " int " << *int_it << std::endl;
          im_chrom[ int((*im_it)*IM_IDX_MULT) ] += *int_it;
          intensity += (*int_it);
          im += (*int_it) * (*im_it);
        }
      }

      if (intensity > 0.)
      {
        // std::cout << " before " << im << std::endl;
        im /= intensity;
        // std::cout << " after " << im << std::endl;
      }
      else
      {
        im = -1;
        intensity = 0;
      }

    }

    for (auto k : im_chrom) 
    {
      res.push_back(std::make_pair( k.first / IM_IDX_MULT, k.second ) );
    }

  }


  /// Constructor
  IonMobilityScoring::IonMobilityScoring()
  {
  }

  /// Destructor
  IonMobilityScoring::~IonMobilityScoring()
  {
  }

  void IonMobilityScoring::driftScoringMS1Contrast(OpenSwath::SpectrumPtr spectrum, OpenSwath::SpectrumPtr ms1spectrum, 
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_lower,
                                        const double drift_upper,
                                        const double drift_target,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const bool /* use_spline */,
                                        const double drift_extra)
  {
    OPENMS_PRECONDITION(spectrum != nullptr, "Spectrum cannot be null");
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");
    OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot score drift time if no drift time is available.");

    double eps = 1e-5; // eps for two grid cells to be considered equal

    double drift_width = fabs(drift_upper - drift_lower);
    double drift_lower_used = drift_lower - drift_width * drift_extra;
    double drift_upper_used = drift_upper + drift_width * drift_extra;

    double delta_drift = 0;
    std::vector< IMProfile > im_profiles;
    double computed_im = 0;
    double computed_im_weighted = 0;
    double sum_intensity = 0;
    int tr_used = 0;

    // Step 1: MS2 extraction
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      double im(0), intensity(0);
      IMProfile res;
      const TransitionType transition = transitions[k];
      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      double left(transition.getProductMZ()), right(transition.getProductMZ());
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);

      integrateDriftSpectrum(spectrum, left, right, im, intensity, res, drift_lower_used, drift_upper_used);
      im_profiles.push_back( std::move(res) );
    }

    // Step 2: MS1 extraction
    double im(0), intensity(0);
    IMProfile ms1_profile;
    double left(transitions[0].getPrecursorMZ()), right(transitions[0].getPrecursorMZ());
    DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
    integrateDriftSpectrum(ms1spectrum, left, right, im, intensity, ms1_profile, drift_lower_used, drift_upper_used); // TODO: aggregate over isotopes
    im_profiles.push_back(ms1_profile);

    std::vector<double> im_grid = computeGrid(im_profiles, eps); // ensure grid is based on all profiles!
    im_profiles.pop_back();

    // Step 3: Align the IMProfile vectors to the grid
    std::vector< std::vector< double > > raw_im_profiles_aligned;
    std::vector<double> delta_im;
    for (const auto & profile : im_profiles) 
    {
      std::vector< double > al_int_values; // intensity values
      std::vector< double > al_im_values; // ion mobility values
      Size max_peak_idx = 0;
      alignToGrid(profile, im_grid, al_int_values, al_im_values, max_peak_idx);
      raw_im_profiles_aligned.push_back(al_int_values);
    }

    std::vector< double > ms1_int_values; // intensity values
    std::vector< double > ms1_im_values; // ion mobility values
    Size max_peak_idx = 0;

    alignToGrid(ms1_profile, im_grid, ms1_int_values, ms1_im_values, max_peak_idx);

    // Step 4: MS1 contrast scores
    {
      OpenSwath::MRMScoring mrmscore_;
      mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, raw_im_profiles_aligned);
      OPENMS_LOG_DEBUG << "all-all: Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastCoelutionScore() << " / shape  precursor " << 
        mrmscore_.calcXcorrPrecursorContrastShapeScore() << std::endl;
      scores.im_ms1_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
      scores.im_ms1_contrast_shape = mrmscore_.calcXcorrPrecursorContrastShapeScore();
    }

    // Step 5: contrast precursor vs summed fragment ions
    std::vector<double> fragment_values;
    fragment_values.resize(ms1_int_values.size(), 0);
    for (Size k = 0; k < fragment_values.size(); k++)
    {
      for (Size i = 0; i < raw_im_profiles_aligned.size(); i++)
      {
        fragment_values[k] += raw_im_profiles_aligned[i][k];
      }
    }
    {
      OpenSwath::MRMScoring mrmscore_;
      mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, {fragment_values});
      OPENMS_LOG_DEBUG << "Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastCoelutionScore() << " / shape  precursor " << 
        mrmscore_.calcXcorrPrecursorContrastShapeScore() << std::endl;
      scores.im_ms1_sum_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
      scores.im_ms1_sum_contrast_shape = mrmscore_.calcXcorrPrecursorContrastShapeScore();
    }

  }

  void IonMobilityScoring::driftScoringMS1(OpenSwath::SpectrumPtr spectrum, 
                                           const std::vector<TransitionType> & transitions,
                                           OpenSwath_Scores & scores,
                                           const double drift_lower,
                                           const double drift_upper,
                                           const double drift_target,
                                           const double dia_extract_window_,
                                           const bool dia_extraction_ppm_,
                                           const bool /* use_spline */,
                                           const double drift_extra)
  {
    OPENMS_PRECONDITION(spectrum != nullptr, "Spectrum cannot be null");
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");
    OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot score drift time if no drift time is available.");

    double drift_width = fabs(drift_upper - drift_lower);
    double drift_lower_used = drift_lower - drift_width * drift_extra;
    double drift_upper_used = drift_upper + drift_width * drift_extra;

    double im(0), intensity(0);
    double left(transitions[0].getPrecursorMZ()), right(transitions[0].getPrecursorMZ());
    DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
    DIAHelpers::integrateDriftSpectrum(spectrum, left, right, im, intensity, drift_lower_used, drift_upper_used);

    // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
    scores.im_ms1_delta_score = fabs(drift_target - im);
  }

  void IonMobilityScoring::driftScoring(OpenSwath::SpectrumPtr spectrum, 
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_lower,
                                        const double drift_upper,
                                        const double drift_target,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const bool use_spline,
                                        const double drift_extra)
  {
    OPENMS_PRECONDITION(spectrum != nullptr, "Spectrum cannot be null");
    OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot score drift time if no drift time is available.");

    double eps = 1e-5; // eps for two grid cells to be considered equal

    double drift_width = fabs(drift_upper - drift_lower);
    double drift_lower_used = drift_lower - drift_width * drift_extra;
    double drift_upper_used = drift_upper + drift_width * drift_extra;

    auto im_range = MSDriftSpectrum::getIMValues(spectrum->getDriftTimeArray()->data);

    // IMProfile: a data structure that holds points <im_value, intensity>
    double delta_drift = 0;
    std::vector< IMProfile > im_profiles;
    double computed_im = 0;
    double computed_im_weighted = 0;
    double sum_intensity = 0;
    int tr_used = 0;

    // Step 1: MS2 extraction
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType transition = transitions[k];
      IMProfile res;
      double im(0), intensity(0);

      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      double left(transition.getProductMZ()), right(transition.getProductMZ());
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      integrateDriftSpectrum(spectrum, left, right, im, intensity, res, drift_lower_used, drift_upper_used);
      im_profiles.push_back(res);

      // TODO what do to about those that have no signal ?
      if (intensity <= 0.0) {continue;} // note: im is -1 then

      tr_used++;

      delta_drift += fabs(drift_target - im);
      OPENMS_LOG_DEBUG << "  -- have delta drift time " << fabs(drift_target -im ) << " with im " << im << std::endl;
      computed_im += im;
      computed_im_weighted += im * intensity;
      sum_intensity += intensity;
      // delta_drift_weighted += delta_drift * normalized_library_intensity[k];
      // weights += normalized_library_intensity[k];
    }
    OPENMS_LOG_DEBUG << " Scoring delta drift time " << delta_drift / tr_used << std::endl;
    scores.im_delta_score = delta_drift / tr_used;

    if (tr_used != 0)
    {
      computed_im /= tr_used;
      computed_im_weighted /= sum_intensity;
    }
    else
    {
      computed_im = -1;
      computed_im_weighted = -1;
    }

    OPENMS_LOG_DEBUG << " Scoring weighted delta drift time " << computed_im_weighted << " -> get difference " << std::fabs(computed_im_weighted - drift_target)<< std::endl;
    scores.im_drift = computed_im;
    scores.im_drift_weighted = computed_im_weighted;

    // Step 2: Align the IMProfile vectors to the grid
    std::vector<double> im_grid = computeGrid(im_profiles, eps);
    std::vector< std::vector< double > > raw_im_profiles_aligned;
    std::vector<double> delta_im;
    for (const auto & profile : im_profiles) 
    {
      std::vector< double > al_int_values; // intensity values
      std::vector< double > al_im_values; // ion mobility values
      Size max_peak_idx = 0;
      alignToGrid(profile, im_grid, al_int_values, al_im_values, max_peak_idx);

      // im_profiles_aligned.push_back(aligned_profile);
      raw_im_profiles_aligned.push_back(al_int_values);

      if (!use_spline) {continue;}

      // We need at least 2 datapoints for the cubic spline
      // (sometimes there are just not enough datapoints available in the spectra)
      if (al_im_values.size() < 2)
      {
        double max_diff = std::max( std::fabs(drift_target - drift_lower), std::fabs(drift_target - drift_upper) );
        delta_im.push_back(max_diff);
        continue;
      }

      // Use cubic spline interpolation to find exact minima / maxima
      CubicSpline2d peak_spline (al_im_values, al_int_values);
      double spline_im(0);
      if (max_peak_idx > 0 && max_peak_idx < al_im_values.size() )
      {
        double spline_int(0);
        OpenMS::Math::spline_bisection(peak_spline, al_im_values[ max_peak_idx - 1], al_im_values[ max_peak_idx + 1], spline_im, spline_int);
      }
      delta_im.push_back(fabs(drift_target - spline_im));
    }

    // Step 3: Compute cross-correlation scores based on ion mobilograms
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrMatrix(raw_im_profiles_aligned);

    double xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
    double xcorr_shape_score = mrmscore_.calcXcorrShapeScore(); // can be nan!

    if (use_spline)
    {
      OpenSwath::mean_and_stddev delta_m;
      delta_m = std::for_each(delta_im.begin(), delta_im.end(), delta_m);
      scores.im_delta_score = delta_m.mean(); // can be nan!
    }

    scores.im_xcorr_coelution_score = xcorr_coelution_score;
    scores.im_xcorr_shape_score = xcorr_shape_score;
  }

}

