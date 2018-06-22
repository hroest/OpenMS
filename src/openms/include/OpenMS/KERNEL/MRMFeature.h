// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>

namespace OpenMS
{
  /**
    @brief A multi-chromatogram MRM feature 

    An MRM feature contains corresponding features in multiple chromatograms,
    it is thus a representation of a peak group. The individual features in
    each chromatogram are represented by OpenMS Features.

  */
  class OPENMS_DLLAPI MRMFeature :
    public Feature
  {
public:

    ///Type definitions
    //@{
    /// Feature list type
    typedef std::vector<Feature> FeatureListType;
    //@}

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    MRMFeature();

    /// Copy constructor
    MRMFeature(const MRMFeature &rhs);

    /// Assignment operator
    MRMFeature & operator=(const MRMFeature & rhs);

    /// Destructor
    ~MRMFeature() override;
    //@}

    ///@name Accessors
    //@{

    const double & getExpectedRT() const
    {
      return expected_rt_;
    }
    void setExpectedRT(double rt)
    {
      expected_rt_ = rt;
    }

    const double & getPeakApex() const
    {
      return peak_apex_;
    }
    void setPeakApex(double pa)
    {
      peak_apex_ = pa;
    }

    /// get all peakgroup scores
    const OpenSwath_Scores & getScores() const;

    /// get all peakgroup scores
    OpenSwath_Scores & getScores();

    /// get a specified feature
    Feature & getFeature(const String& key);

    /// get a specified feature (const)
    const Feature & getFeature(const String& key) const;

    /// set all peakgroup scores
    void setScores(const OpenSwath_Scores & scores);

    /// Adds an feature from a single chromatogram into the feature.
    void addFeature(const Feature & feature, const String& key);

    /// get a list of features
    const std::vector<Feature> & getFeatures() const;

    /// get a list of IDs of available features
    void getFeatureIDs(std::vector<String> & result) const;

    /// Adds a precursor feature from a single chromatogram into the feature.
    void addPrecursorFeature(const Feature & feature, const String& key);

    /// get a list of IDs of available precursor features
    void getPrecursorFeatureIDs(std::vector<String> & result) const;

    /// get a specified precursor feature
    Feature & getPrecursorFeature(String key);

    /// get a specified precursor feature (const)
    const Feature & getPrecursorFeature(String key) const;

    void scoresAsMetaValue(bool ms1only, const OpenSwath_Scores_Usage&);

    void IDScoresAsMetaValue(bool decoy, const OpenSwath_Ind_Scores& idscores);
    //@}

protected:

    FeatureListType features_;

    double expected_rt_;
    double peak_apex_;

    FeatureListType precursor_features_;

    /// OpenSwathScores
    OpenSwath_Scores scores_;

    /// map native ids to the features
    std::map<String, int> feature_map_;

    /// map native ids to the precursor features
    std::map<String, int> precursor_feature_map_;

  };
}

