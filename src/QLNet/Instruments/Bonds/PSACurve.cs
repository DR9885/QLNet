/*
 Copyright (C) 2008, 2009 , 2010, 2011, 2012  Andrea Maggiulli (a.maggiulli@gmail.com)

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
   public class PSACurve : IPrepayModel
   {

      public PSACurve(Date startdate)
         : this(startdate, Const.ONE_INT) {}

      public PSACurve(Date startdate, double multiplier)
      {
         _startDate = startdate;
         _multi = multiplier;
      }

      public double getCPR(Date valDate)
      {
         Thirty360 dayCounter = new Thirty360();
         int d = dayCounter.dayCount(_startDate, valDate) / Const.THIRTY_INT + Const.ONE_INT;

         return (d <= Const.THIRTY_INT ? Const.SIX_PERCENT * (d / Const.THIRTY_DOUBLE) : Const.SIX_PERCENT) * _multi;
      }

      public double getSMM(Date valDate)
      {
         return Const.ONE_INT - Math.Pow((Const.ONE_INT - getCPR(valDate)), (Const.ONE_INT / Const.TWELVE_DOUBLE));
      }

      private Date _startDate;
      private double _multi;
   }
}
