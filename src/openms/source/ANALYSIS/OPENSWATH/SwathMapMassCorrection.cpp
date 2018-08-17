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
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>

// Classes
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>

// Functions
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoringHelper.h>

#define SWATHMAPMASSCORRECTION_DEBUG

namespace OpenMS
{

  void integrateDriftSpectrum_x(OpenSwath::SpectrumPtr spectrum, 
                              double mz_start,
                              double mz_end,
                              double & im,
                              double & intensity,
                              std::vector<std::pair<double, double> >& res, 
                              double drift_start,
                              double drift_end)
  {
    OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");
    im = 0; intensity = 0;

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
          // std::cout << "IM " << *im_it << " mz " << *mz_it << " int " << *int_it << " cum " << im << std::endl;
          im_chrom[ int((*im_it)*IM_IDX_MULT) ] += *int_it;
          intensity += (*int_it);
          im += (*int_it) * (*im_it);
        }
      }

      if (intensity > 0.)
      {
        im /= intensity;
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

  void SwathMapMassCorrection::correctMZ(
    const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    TransformationDescription& im_trafo,
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    const std::string& corr_type,
    const double mz_extr_window,
    const bool ppm)
  {
    LOG_DEBUG << "SwathMapMassCorrection::correctMZ with type " << corr_type << " and window " << mz_extr_window << " in ppm " << ppm << std::endl;

    bool is_ppm = bool(corr_type == "quadratic_regression_delta_ppm" ||
                       corr_type == "weighted_quadratic_regression_delta_ppm" ||
                       corr_type == "regression_delta_ppm");

    if (corr_type == "none")
    {
      return;
    }

#ifdef SWATHMAPMASSCORRECTION_DEBUG
    std::cout.precision(16);
    std::ofstream os("debug_ppmdiff.txt");
    os.precision(writtenDigits(double()));

    std::cout.precision(16);
    std::ofstream os_im("debug_imdiff.txt");
    os_im.precision(writtenDigits(double()));
#endif

    TransformationDescription::DataPoints data_all;
    std::vector<double> weights;
    std::vector<double> exp_mz;
    std::vector<double> theo_mz;
    std::vector<double> delta_ppm;

    TransformationDescription::DataPoints data_im;
    std::vector<double> exp_im;
    std::vector<double> theo_im;
    for (auto trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); ++trgroup_it)
    {

      // we need at least one feature to find the best one
      auto transition_group = trgroup_it->second;
      if (transition_group->getFeatures().size() == 0)
      {
        continue;
      }

      // Find the feature with the highest score
      double bestRT = -1;
      double highest_score = -1000;
      for (auto mrmfeature = transition_group->getFeatures().begin(); mrmfeature != transition_group->getFeatures().end(); ++mrmfeature)
      {
        if (mrmfeature->getOverallQuality() > highest_score)
        {
          bestRT = mrmfeature->getRT();
          highest_score = mrmfeature->getOverallQuality();
        }
      }

      // Get the corresponding SWATH map(s), for SONAR there will be more than one map
      std::vector<OpenSwath::SwathMap> used_maps;
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        if (swath_maps[i].lower < transition_group->getTransitions()[0].precursor_mz &&
            swath_maps[i].upper >= transition_group->getTransitions()[0].precursor_mz)
        {
          used_maps.push_back(swath_maps[i]);
        }
      }

      if (used_maps.empty())
      {
        continue;
      }

      // Get the spectrum for this RT and extract raw data points for all the
      // calibrating transitions (fragment m/z values) from the spectrum
      OpenSwath::SpectrumPtr sp = OpenSwathScoring().fetchSpectrumSwath(used_maps, bestRT, 1, 0, 0);
      for (std::vector< OpenMS::MRMFeatureFinderScoring::TransitionType >::const_iterator
          tr = transition_group->getTransitions().begin();
          tr != transition_group->getTransitions().end(); ++tr)
      {
        double mz, intensity;
        double left = tr->product_mz - mz_extr_window / 2.0;
        double right = tr->product_mz + mz_extr_window / 2.0;
        bool centroided = false;

        auto pepref = tr->getPeptideRef();
        auto kk = const_cast<OpenSwath::LightTargetedExperiment&>(targeted_exp);

        auto pep = kk.getPeptideByRef(pepref);

        if (ppm)
        {
          left = tr->product_mz - mz_extr_window / 2.0  * tr->product_mz * 1e-6;
          right = tr->product_mz + mz_extr_window / 2.0 * tr->product_mz * 1e-6;
        }

        // integrate spectrum at the position of the theoretical mass
        DIAHelpers::integrateWindow(sp, left, right, mz, intensity, centroided);

        // IMProfile: a data structure that holds points <im_value, intensity>
        typedef std::vector< std::pair<double, double> > IMProfile;
        double delta_drift = 0;
        std::vector< IMProfile > im_profiles;

        double drift_lower(0), drift_upper(0), drift_target(0);
        if (!transition_group->getChromatograms().empty())
        {
          auto & prec = transition_group->getChromatograms()[0].getPrecursor();
          drift_lower = prec.getDriftTime() - prec.getDriftTimeWindowLowerOffset();
          drift_upper = prec.getDriftTime() + prec.getDriftTimeWindowUpperOffset();
          drift_target = prec.getDriftTime(); 
        }
        drift_target = pep.getDriftTime();

        // double left(transition->getProductMZ()), right(transition->getProductMZ());
        // adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
        IMProfile res;
        double im(0);
        intensity = 0;
        double drift_lower_used(drift_target - 0.1), drift_upper_used(drift_target + 0.1);
        integrateDriftSpectrum_x(sp, left, right, im, intensity, res, drift_lower_used, drift_upper_used);
        // std::cout << " extract " << pepref << " at " << drift_target << " from " << drift_lower << " to " << drift_upper << std::endl;
        // std::cout << " found " << im << " vs " << drift_target << std::endl;

        // skip empty windows
        if (mz == -1)
        {
          continue;
        }

        // store result masses
        data_im.push_back(std::make_pair(im, drift_target));
        exp_im.push_back(im);
        theo_im.push_back(drift_target);

        // store result masses
        data_all.push_back(std::make_pair(mz, tr->product_mz));
        // regression weight is the log2 intensity
        weights.push_back( log(intensity) / log(2.0) );
        exp_mz.push_back( mz );
        // y = target = theoretical
        theo_mz.push_back( tr->product_mz );
        double diff_ppm = (mz - tr->product_mz) * 1000000 / mz;
        // y = target = delta-ppm
        delta_ppm.push_back(diff_ppm);

#ifdef SWATHMAPMASSCORRECTION_DEBUG
        os << mz << "\t" << tr->product_mz << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;
        os_im << mz << "\t" << im << "\t" << drift_target << "\t" << bestRT << std::endl;
        // std::cout << mz << "\t" << tr->product_mz << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;
#endif
        LOG_DEBUG << mz << "\t" << tr->product_mz << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;

      }
    }

    std::vector<double> regression_params;
    if (corr_type == "none" || data_all.size() < 3)
    {
      return;
    }
    else if (corr_type == "unweighted_regression")
    {
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_regression")
    {
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegressionWeighted(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "quadratic_regression")
    {
      // Quadratic fit
      Math::QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "weighted_quadratic_regression")
    {
      // Quadratic fit (weighted)
      Math::QuadraticRegression qr;
      qr.computeRegressionWeighted(exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      Math::QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "regression_delta_ppm")
    {
      // Regression fit using ppm differences
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      Math::QuadraticRegression qr;
      qr.computeRegressionWeighted(exp_mz.begin(), exp_mz.end(), delta_ppm.begin(), weights.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Unknown correction type " + corr_type);
    }

    printf("# mz regression parameters: Y = %g + %g X + %g X^2\n",
           regression_params[0],
           regression_params[1],
           regression_params[2]);

    LOG_DEBUG << "# mz regression parameters: Y = " << regression_params[0] << " + " <<
      regression_params[1] << " X + " << regression_params[2] << " X^2" << std::endl;

#ifdef SWATHMAPMASSCORRECTION_DEBUG
    os.close();
    double s_ppm_before = 0;
    double s_ppm_after = 0;
    for (TransformationDescription::DataPoints::iterator d = data_all.begin(); d != data_all.end(); ++d)
    {
      double ppm_before = (d->first - d->second) * 1000000 / d->first;
      double predict = d->first*d->first*regression_params[2] + d->first*regression_params[1]+regression_params[0];
      double ppm_after = ( predict - d->second) * 1000000 / d->first;
      if (is_ppm)
      {
        double new_mz = d->first - predict*d->first/1000000;
        ppm_after = ( new_mz - d->second) * 1000000 / d->first;
      }
      s_ppm_before += std::fabs(ppm_before);
      s_ppm_after += std::fabs(ppm_after);
    }
    std::cout <<" sum residual sq ppm before " << s_ppm_before << " / after " << s_ppm_after << std::endl;
#endif

    // Replace the swath files with a transforming wrapper.
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      swath_maps[i].sptr = boost::shared_ptr<OpenSwath::ISpectrumAccess>(
        new SpectrumAccessQuadMZTransforming(swath_maps[i].sptr,
          regression_params[0], regression_params[1], regression_params[2], is_ppm));
    }

    {
      std::vector<double> im_regression_params;
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_im.begin(), exp_im.end(), theo_im.begin()); // to convert exp_im -> theoretical im
      im_regression_params.push_back(lr.getIntercept());
      im_regression_params.push_back(lr.getSlope());
      im_regression_params.push_back(0.0);

      std::cout << "# im regression parameters: Y = " << im_regression_params[0] << " + " <<
        im_regression_params[1] << " X + " << im_regression_params[2] << " X^2" << std::endl;

      // store IM transformation, using the selected model
      im_trafo.setDataPoints(data_im);
      Param model_params;
      model_params.setValue("symmetric_regression", "false");
      // model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
      // model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
      String model_type = "linear";
      im_trafo.fitModel(model_type, model_params);

      // std::cout << " trafo transforms 1 to " << im_trafo.apply(1.0) << std::endl; // converts experimental IM to theorteical IM
    }

    LOG_DEBUG << "SwathMapMassCorrection::correctMZ done." << std::endl;
  }

/* display debug output in R

# os << mz << "\t" << tr->product_mz << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;
 df = read.csv("debug_ppmdiff.txt", sep="\t" , header=F)


 pdf("/tmp/mz.pdf")
 df = read.csv("/tmp/debug_ppmdiff.txt", sep="\t" , header=F)
 colnames(df) = c("mz", "theomz", "dppm", "int", "rt")
 df$absd = df$theomz-df$mz
 plot(df$mz, df$absd, main="m/z vs absolute deviation")
 plot(df$mz, df$dppm, main="m/z vs ppm deviation")
 hist(df$absd, main="absolute deviation (before)", breaks=45)
 abline(v=median(df$absd), col="red")
 hist(df$dppm, main="deviation ppm (before)", breaks=45)
 abline(v=median(df$dppm), col="red")

 linm = lm(theomz ~ mz, df)
 df$x2 = df$mz*df$mz
 quadm = lm(theomz ~ mz + x2, df)

 plot(df$mz, df$dppm, cex=0.5)
 quadm_ppm = lm(dppm ~ mz + x2, df)

 plot(df$mz, df$dppm, main="PPM difference and model fit")
 points(df$mz, predict(quadm_ppm), cex=0.3, col="red")

 plot(resid(linm), main="m/z residual (after linear fit)")
 plot(resid(quadm), main="m/z residual (after quadratic fit)")
 plot(resid(quadm_ppm), main="m/z residual (after quadratic ppm fit)")

 hist(resid(linm), main="m/z residual (after linear fit)", breaks=45)
 hist(resid(quadm), main="m/z residual (after quadratic fit)", breaks=45)
 hist(resid(quadm_ppm), main="m/z residual (after quadratic fit)", breaks=45)

 summary(resid(quadm_ppm))
 quantile(resid(quadm_ppm), probs=c(0.05, 0.95))

 



 sd(df$dppm)
 mean(df$dppm)
 sd(df$ppmdiff_pred)
 mean(df$ppmdiff_pred)
 sd(df$dppm - predict(quadm_ppm))
 mean(df$dppm - predict(quadm_ppm))



*/

/* display debug output in R

 pdf("/tmp/im_calibration.pdf")
 df = read.csv("/tmp/debug_imdiff.txt", sep="\t" , header=F)

 colnames(df) = c("mz", "im", "theo_im", "rt")
 df$absd = df$theo_im - df$im
 linm = lm(theo_im ~ im, df)

 plot(df$im, df$theo_im, main="Ion Mobility (before)")
 abline(0, 1)
 abline(lm(theo_im ~ im, df), col="red")

 plot(df$im, df$absd, main="Ion Mobility residual (before)")
 hist(df$absd, main="Ion Mobility residual (before)", breaks=45)

 plot(resid(linm), main="Ion Mobility residual (after fit)")
 hist(resid(linm), main="Ion Mobility residual (after fit)", breaks=45)
 hist(resid(linm), main="Ion Mobility residual (after fit)", breaks=45, xlim=c(-0.03, 0.03))

 scale = max(df$theo_im) - min(df$theo_im)
 hist(resid(linm)/scale * 100, main="Ion Mobility residual (after fit)", breaks=45)
 hist(resid(linm)/scale * 100, main="Ion Mobility residual (after fit)", breaks=45, xlim=c(-3,3))
 summary(resid(linm)/scale * 100)

 summary(resid(linm))
 quantile(resid(linm), probs=c(0.05, 0.95))

 dev.off()

*/

/*


# mz regression parameters: Y = -0.0156725 + 1.00005 X + -3.7127e-08 X^2
sum residual sq ppm before 394.0529252540293 / after 245.4905283249099
rsq: 0.9993098240025888 points: 7

# mz regression parameters: Y = 42.2354 + -0.107071 X + 6.87643e-05 X^2
sum residual sq ppm before 394.0529252540293 / after 205.778656429406
rsq: 0.9993098240025888 points: 7


-- done [took 4.53 s (CPU), 4.55 s (Wall)] --
OpenSwathWorkflow took 23.92 s (wall), 23.73 s (CPU), 0.00 s (system), 23.73 s (user).

-- done [took 1.41 s (CPU), 1.41 s (Wall)] --
OpenSwathWorkflow took 21.07 s (wall), 20.90 s (CPU), 0.00 s (system), 20.90 s (user).

*/

}

