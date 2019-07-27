// -*- C++ -*-
#ifndef BUNDLE_LARGE_COMMON_H
#define BUNDLE_LARGE_COMMON_H

#include "Geometry/v3d_cameramatrix.h"

//#define USE_TUKEYS_BIWEIGHT

namespace V3D
{

   struct SimpleDistortionFunction
   {
         double k1, k2;

         SimpleDistortionFunction()
            : k1(0), k2(0)
         { }

         Vector2d operator()(Vector2d const& xu) const
         {
            double const r2 = xu[0]*xu[0] + xu[1]*xu[1];
            double const r4 = r2*r2;
            double const kr = 1 + k1*r2 + k2*r4;

            Vector2d xd;
            xd[0] = kr * xu[0];
            xd[1] = kr * xu[1];
            return xd;
         }

         Matrix2x2d derivativeWrtRadialParameters(Vector2d const& xu) const
         {
            double const r2 = xu[0]*xu[0] + xu[1]*xu[1];
            double const r4 = r2*r2;

            Matrix2x2d deriv;

            deriv[0][0] = xu[0] * r2; // d xd/d k1
            deriv[0][1] = xu[0] * r4; // d xd/d k2
            deriv[1][0] = xu[1] * r2; // d yd/d k1
            deriv[1][1] = xu[1] * r4; // d yd/d k2
            return deriv;
         }

         Matrix2x2d derivativeWrtUndistortedPoint(Vector2d const& xu) const
         {
            double const r2 = xu[0]*xu[0] + xu[1]*xu[1];
            double const r4 = r2*r2;
            double const kr = 1 + k1*r2 + k2*r4;
            double const dkr = 2*k1 + 4*k2*r2;

            Matrix2x2d deriv;
            deriv[0][0] = kr + xu[0] * xu[0] * dkr; // d xd/d xu
            deriv[0][1] =      xu[0] * xu[1] * dkr; // d xd/d yu
            deriv[1][0] = deriv[0][1];              // d yd/d xu
            deriv[1][1] = kr + xu[1] * xu[1] * dkr; // d yd/d yu
            return deriv;
         }
   }; // end struct SimpleDistortionFunction

//**********************************************************************

   inline double sqr(double const x) { return x*x; }

#if !defined(USE_TUKEYS_BIWEIGHT)
   inline double psi(double const tau2, double const r2) { return (r2 < tau2) ? r2*(2.0 - r2/tau2)/4.0f : tau2/4; }
   inline double psi_weight(double const tau2, double const r2) { return std::max(0.0, 1.0 - r2/tau2); }
   inline double psi_hat(double const tau2, double const r2, double const w2) { return w2*r2 + tau2/2.0*(w2-1)*(w2-1); }
#else
   inline double psi(double const tau2, double const r2)
   {
      double const r4 = r2*r2, tau4 = tau2*tau2;
      return (r2 < tau2) ? r2*(3.0 - 3*r2/tau2 + r4/tau4)/6.0f : tau2/6.0;
   }
   inline double psi_weight(double const tau2, double const r2)
   {
      return sqr(std::max(0.0, 1.0 - r2/tau2));
   }
   inline double psi_hat(double const tau2, double const r2, double const w2)
   {
      double const w = sqrt(w2);
      return w2*r2 + tau2/3.0*sqr(w-1)*(2*w+1);
   }
#endif


// lifted schur cost define start

#if !defined(USE_TUKEYS_BIWEIGHT)
   inline double kappa(double const tau, double const w2) { return 0.70710678118*tau*(w2 - 1); }
   inline double dkappa(double const tau, double const w2) { return 0.70710678118*tau; }
#else
   inline double kappa(double const tau, double const w2)
   {
      double const f = sqrt(1.0/3);
      double const w = sqrt(w2);
      return f * tau * (w - 1)*sqrt(2*w+1);
   }
   inline double dkappa(double const tau, double const w2)
   {
      double const f = sqrt(3.0)/2;
      double const w = sqrt(w2);
      return f * tau / sqrt(2*w+1);
   }
#endif

// The mapping w |-> h(w)
#if 1
   inline double omega(double const w) { return w; }
   inline double domega(double const w) { return 1.0; }
   inline double omega2_inv(double const w) { return sqrt(w); }
#elif 0
   //double const slope_w = 0.25; //3.0;
   double const slope_w = 4.0;
   inline double omega(double const w) { return slope_w*w; }
   inline double domega(double const w) { return slope_w; }
   inline double omega2_inv(double const w) { return sqrt(w)/slope_w; }
#elif 1
   inline double omega(double const w) { return 0.5*(w/sqrt(1+w*w)+1); }
   inline double domega(double const w) { return 0.5*sqrt(1+w*w)/(1+w*w); }
   inline double omega2_inv(double const w) { double const ww = sqrt(w); return (2*ww-1.0)/sqrt(1.0 - sqr(2*ww-1)); }
// #elif 1
//    inline double omega(double const w) { return w*w; }
//    inline double domega(double const w) { return 2.0*w; }
//    inline double omega2_inv(double const w) { return sqrt(sqrt(w)); }
// #elif 0
//    inline double omega(double const w) { return sqr(std::max(0.0, w)); }
//    inline double domega(double const w) { return (w > 0) ? 2.0*w : 0.0; }
//    inline double omega2_inv(double const w) { return sqrt(sqrt(w)); }
#elif 1
   // h(w) = exp(w)
   double const exp_const = 1.0/4.0;
   inline double omega(double const w) { return exp(w * exp_const); }
   inline double domega(double const w) { return exp_const * exp(w * exp_const); }
   inline double omega2_inv(double const w) { return log(w)/2.0/exp_const; }
#else
   // h(w) = soft-max(0, w), logistic
   double const exp_const = 4.0;
   inline double omega(double const w) { return exp_const * log(1.0 + exp(w/exp_const)); }
   inline double domega(double const w) { double const exp_w = exp(w/exp_const); return exp_w/(1.0 + exp_w); }
   inline double omega2_inv(double const w) { return exp_const * log(exp(sqrt(w)/exp_const) - 1.0); }
#endif

//  lifted schur cost define end


   inline double
   showErrorStatistics(double const avg_focal_length, double const inlierThreshold,
                       vector<CameraMatrix> const& cams,
                       vector<SimpleDistortionFunction> const& distortions,
                       vector<Vector3d> const& Xs,
                       vector<Vector2d> const& measurements,
                       vector<int> const& correspondingView,
                       vector<int> const& correspondingPoint)
   {
      int const K = measurements.size();

      int nInliers = 0;

      double meanReprojectionError = 0.0, inlierReprojectionError = 0.0;
      for (int k = 0; k < K; ++k)
      {
         int const i = correspondingView[k];
         int const j = correspondingPoint[k];
         Vector2d p = cams[i].projectPoint(distortions[i], Xs[j]);

         double const reprojectionError = avg_focal_length * norm_L2(p - measurements[k]);
         meanReprojectionError += reprojectionError;
         if (reprojectionError <= inlierThreshold)
         {
            ++nInliers;
            inlierReprojectionError += reprojectionError;
         }
      }
      cout << "mean reprojection error = " << meanReprojectionError/K << endl;
      cout << "inlier mean reprojection error = " << inlierReprojectionError/nInliers << " " << nInliers << " / " << K << " inliers." << endl;
      return double(nInliers) / K;
   }

   inline double
   showObjective(double const avg_focal_length, double const inlierThreshold,
                 vector<CameraMatrix> const& cams,
                 vector<SimpleDistortionFunction> const& distortions,
                 vector<Vector3d> const& Xs,
                 vector<Vector2d> const& measurements,
                 vector<int> const& correspondingView,
                 vector<int> const& correspondingPoint)
   {
      
      int const K = measurements.size();

      double const tau2 = inlierThreshold*inlierThreshold;
      double const avg_focal_length2 = avg_focal_length*avg_focal_length;

      double obj = 0.0;
      double E_data = 0.0, E_reg = 0.0;
      
      for (int k = 0; k < K; ++k)
      {
         int const i = correspondingView[k];
         int const j = correspondingPoint[k];
         Vector2d p = cams[i].projectPoint(distortions[i], Xs[j]);

         double const r2 = avg_focal_length2 * sqrNorm_L2(p - measurements[k]);
         obj += psi(tau2, r2);

         // for lifted schur objective
         
         //double const r2 = avg_focal_length2*sqrNorm_L2(p - measurements[k]);
         //double const weights = psi_weight(tau2, r2);
         //double const w2 = sqr(omega(weights));
         //E_data += w2*r2;
         //E_reg += sqr(kappa(inlierThreshold, w2));
         //E_data /= 2.0; E_reg /= 2.0;
         //obj += psi_hat(tau2, r2, w2);
         //obj += psi(tau2, r2);
         //obj += r2; 
         //obj += E_data + E_reg;
      }
      //cout << "true objective = " << obj << endl;
      return obj;
      

      // Computing lifted schur objective

      /*
      double E_data = 0.0, E_reg = 0.0;

      for (int k = 0; k < _nMeasurements; ++k)
      {
         double const w2 = sqr(omega(_weights[k]));
         //double const w2 = sqr(omega(_weights[k]));
         E_data += w2*sqrNorm_L2(residuals[k]);
         E_reg += sqr(kappa(_inlierThreshold, w2));
      }
      E_data /= 2.0; E_reg /= 2.0;
      //cout << "evalCost(): E_data = " << E_data << " E_reg = " << E_reg << endl;
      obj = E_data + E_reg;

      return obj; 
      */
   }

//**********************************************************************

   //double const AVG_FOCAL_LENGTH = 1000.0;
   double const AVG_FOCAL_LENGTH = 1.0;

   enum
   {
      FULL_BUNDLE_METRIC = 0,
      FULL_BUNDLE_FOCAL_LENGTH = 1, // f
      FULL_BUNDLE_RADIAL = 2,       // f, k1, k2
   };

   //int const bundle_mode = FULL_BUNDLE_METRIC;
   //int const bundle_mode = FULL_BUNDLE_FOCAL_LENGTH;
   int const bundle_mode = FULL_BUNDLE_RADIAL;

   //double const inlier_threshold = 1.0;
   double const inlier_threshold = 0.5;

} // end namespace V3D

#endif
