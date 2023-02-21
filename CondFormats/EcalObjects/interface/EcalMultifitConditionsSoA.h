#ifndef CondFormats_EcalObjects_EcalMultifitConditionsSoA_h
#define CondFormats_EcalObjects_EcalMultifitConditionsSoA_h

#include <Eigen/Dense>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"


using PulseShapeArray = std::array<float, EcalPulseShape::TEMPLATESAMPLES>;
using SampleCorrelationArray = std::array<float, ecalPh1::sampleSize>;

using CovarianceMatrix = Eigen::Matrix<float, EcalPulseShape::TEMPLATESAMPLES, EcalPulseShape::TEMPLATESAMPLES>;



GENERATE_SOA_LAYOUT(EcalMultifitConditionsSoALayout,
                    SOA_COLUMN(uint32_t, rawid),
                    SOA_COLUMN(float, pedestals_mean_x12),
                    SOA_COLUMN(float, pedestals_mean_x6),
                    SOA_COLUMN(float, pedestals_mean_x1),
                    SOA_COLUMN(float, pedestals_rms_x12),
                    SOA_COLUMN(float, pedestals_rms_x6),
                    SOA_COLUMN(float, pedestals_rms_x1),
                    SOA_COLUMN(PulseShapeArray, pulse_shapes),
                    SOA_EIGEN_COLUMN(CovarianceMatrix, pulse_covariance),
                    SOA_COLUMN(float, gainRatios_gain12Over6),
                    SOA_COLUMN(float, gainRatios_gain6Over1),
                    // Sample correlation scalar: array of 10 values for each gain in EB and EE
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EB_G12),
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EB_G6),
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EB_G1),
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EE_G12),
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EE_G6),
                    SOA_SCALAR(SampleCorrelationArray, sample_correlation_EE_G1)

  )

using EcalMultifitConditionsSoA = EcalMultifitConditionsSoALayout<>;


#endif
