/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2008-2016 Andrea Maggiulli (a.maggiulli@gmail.com)

 This file is part of QLNet Project https://github.com/amaggiulli/qlnet

 QLNet is free software: you can redistribute it and/or modify it
 under the terms of the QLNet license.  You should have received a
 copy of the license along with this program; if not, license is
 available at <https://github.com/amaggiulli/QLNet/blob/develop/LICENSE>.

 QLNet is a based on QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 The QuantLib license is available online at http://quantlib.org/license.shtml.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/
using System;


namespace QLNet
{
   public class Brent : Solver1D
   {
      protected override double solveImpl(ISolver1d f, double xAccuracy)
      {
         /* The implementation of the algorithm was inspired by Press, Teukolsky, Vetterling, and Flannery,
            "Numerical Recipes in C", 2nd edition, Cambridge University Press */

         double min1, min2;
         double froot, p, q, r, s, xAcc1, xMid;
         // dummy assignements to avoid compiler warning
         double d = Const.ZERO_DOUBLE, e = Const.ZERO_DOUBLE;

         root_ = xMax_;
         froot = fxMax_;
         while (evaluationNumber_ <= maxEvaluations_)
         {
            if ((froot > Const.ZERO_DOUBLE && fxMax_ > Const.ZERO_DOUBLE) ||
                (froot < Const.ZERO_DOUBLE && fxMax_ < Const.ZERO_DOUBLE))
            {

               // Rename xMin_, root_, xMax_ and adjust bounds
               xMax_ = xMin_;
               fxMax_ = fxMin_;
               e = d = root_ - xMin_;
            }
            if (Math.Abs(fxMax_) < Math.Abs(froot))
            {
               xMin_ = root_;
               root_ = xMax_;
               xMax_ = xMin_;
               fxMin_ = froot;
               froot = fxMax_;
               fxMax_ = fxMin_;
            }
            // Convergence check
            xAcc1 = Const.TWO_DOUBLE * Const.QL_EPSILON * Math.Abs(root_) + Const.FIFTY_PERCENT * xAccuracy;
            xMid = (xMax_ - root_) / Const.TWO_DOUBLE;
            if (Math.Abs(xMid) <= xAcc1 || Utils.close(froot, Const.ZERO_DOUBLE))
               return root_;
            if (Math.Abs(e) >= xAcc1 &&
                Math.Abs(fxMin_) > Math.Abs(froot))
            {

               // Attempt inverse quadratic interpolation
               s = froot / fxMin_;
               if (Utils.close(xMin_, xMax_))
               {
                  p = Const.TWO_DOUBLE * xMid * s;
                  q = Const.ONE_DOUBLE - s;
               }
               else
               {
                  q = fxMin_ / fxMax_;
                  r = froot / fxMax_;
                  p = s * (Const.TWO_DOUBLE * xMid * q * (q - r) - (root_ - xMin_) * (r - Const.ONE_DOUBLE));
                  q = (q - Const.ONE_DOUBLE) * (r - Const.ONE_DOUBLE) * (s - Const.ONE_DOUBLE);
               }
               if (p > Const.ZERO_DOUBLE)
                  q = -q;  // Check whether in bounds
               p = Math.Abs(p);
               min1 = Const.THREE_DOUBLE * xMid * q - Math.Abs(xAcc1 * q);
               min2 = Math.Abs(e * q);
               if (Const.TWO_DOUBLE * p < (min1 < min2 ? min1 : min2))
               {
                  e = d;                // Accept interpolation
                  d = p / q;
               }
               else
               {
                  d = xMid;  // Interpolation failed, use bisection
                  e = d;
               }
            }
            else
            {
               // Bounds decreasing too slowly, use bisection
               d = xMid;
               e = d;
            }
            xMin_ = root_;
            fxMin_ = froot;
            if (Math.Abs(d) > xAcc1)
               root_ += d;
            else
               root_ += Math.Abs(xAcc1) * Math.Sign(xMid);
            froot = f.value(root_);
            evaluationNumber_++;
         }
         Utils.QL_FAIL("maximum number of function evaluations (" + maxEvaluations_ + ") exceeded",
                       QLNetExceptionEnum.MaxNumberFuncEvalExceeded);
         return Const.ZERO_INT;
      }
   }
}
