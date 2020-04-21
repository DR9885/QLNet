/*
 Copyright (C) 2008, 2009 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2008-2013 Andrea Maggiulli (a.maggiulli@gmail.com)

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
using System.Collections.Generic;
using System.Linq;

namespace QLNet
{
   //! %Coupon paying a fixed interest rate
   public class FixedRateCoupon : Coupon
   {
      // constructors
      public FixedRateCoupon(Date paymentDate, double nominal, double rate, DayCounter dayCounter,
                             Date accrualStartDate, Date accrualEndDate,
                             Date refPeriodStart = null, Date refPeriodEnd = null, Date exCouponDate = null)
         : base(paymentDate, nominal, accrualStartDate, accrualEndDate, refPeriodStart, refPeriodEnd, exCouponDate)
      {
         rate_ = new InterestRate(rate, dayCounter, Compounding.Simple, Frequency.Annual);
      }

      public FixedRateCoupon(Date paymentDate, double nominal, InterestRate interestRate,
                             Date accrualStartDate, Date accrualEndDate,
                             Date refPeriodStart = null, Date refPeriodEnd = null, Date exCouponDate = null, double? amount = null)
      : base(paymentDate, nominal, accrualStartDate, accrualEndDate, refPeriodStart, refPeriodEnd, exCouponDate)
      {
         amount_ = amount;
         rate_ = interestRate;
      }

      //! CashFlow interface
      public override double amount()
      {
         if (amount_ != null)
            return amount_.Value;

         return nominal() * (rate_.compoundFactor(accrualStartDate_, accrualEndDate_, refPeriodStart_, refPeriodEnd_) - Const.ONE_DOUBLE);
      }

      //! Coupon interface
      public override double rate() { return rate_.rate(); }
      public InterestRate interestRate() { return rate_; }
      public override DayCounter dayCounter() { return rate_.dayCounter(); }
      public override double accruedAmount(Date d)
      {
         if (d <= accrualStartDate_ || d > paymentDate_)
            return Const.ZERO_INT;
         else if (tradingExCoupon(d))
         {
            return -nominal() * (rate_.compoundFactor(d,
                                                      accrualEndDate_,
                                                      refPeriodStart_,
                                                      refPeriodEnd_) - Const.ONE_DOUBLE);
         }
         else
            return nominal() * (rate_.compoundFactor(accrualStartDate_, Date.Min(d, accrualEndDate_),
                                                     refPeriodStart_, refPeriodEnd_) - Const.ONE_DOUBLE);
      }

      private InterestRate rate_;
      private double? amount_;

   }

   //! helper class building a sequence of fixed rate coupons
   public class FixedRateLeg : RateLegBase
   {
      // properties
      private List<InterestRate> couponRates_ = new List<InterestRate>();
      private DayCounter firstPeriodDC_, lastPeriodDC_ ;
      private Calendar calendar_;
      private Period exCouponPeriod_;
      private   Calendar exCouponCalendar_;
      private   BusinessDayConvention exCouponAdjustment_;
      private   bool exCouponEndOfMonth_;

      // constructor
      public FixedRateLeg(Schedule schedule)
      {
         schedule_ = schedule;
         calendar_ = schedule.calendar();
         paymentAdjustment_ = BusinessDayConvention.Following;
      }

      // other initializers
      public FixedRateLeg withCouponRates(double couponRate, DayCounter paymentDayCounter)
      {
         return withCouponRates(couponRate, paymentDayCounter, Compounding.Simple, Frequency.Annual);
      }
      public FixedRateLeg withCouponRates(double couponRate, DayCounter paymentDayCounter, Compounding comp)
      {
         return withCouponRates(couponRate, paymentDayCounter, comp, Frequency.Annual);
      }

      public FixedRateLeg withCouponRates(double couponRate, DayCounter paymentDayCounter,
                                          Compounding comp, Frequency freq)
      {
         couponRates_.Clear();
         couponRates_.Add(new InterestRate(couponRate, paymentDayCounter, comp, freq));
         return this;
      }


      public FixedRateLeg withCouponRates(List<double> couponRates, DayCounter paymentDayCounter)
      {
         return withCouponRates(couponRates, paymentDayCounter, Compounding.Simple, Frequency.Annual);
      }
      public FixedRateLeg withCouponRates(List<double> couponRates, DayCounter paymentDayCounter, Compounding comp)
      {
         return withCouponRates(couponRates, paymentDayCounter, comp, Frequency.Annual);
      }

      public FixedRateLeg withCouponRates(List<double> couponRates, DayCounter paymentDayCounter,
                                          Compounding comp, Frequency freq)
      {
         couponRates_.Clear();
         foreach (double r in couponRates)
            couponRates_.Add(new InterestRate(r, paymentDayCounter, comp, freq));
         return this;
      }

      public FixedRateLeg withCouponRates(InterestRate couponRate)
      {
         couponRates_.Clear();
         couponRates_.Add(couponRate);
         return this;
      }

      public FixedRateLeg withCouponRates(List<InterestRate>couponRates)
      {
         couponRates_ = couponRates;
         return this;
      }

      public FixedRateLeg withFirstPeriodDayCounter(DayCounter dayCounter)
      {
         firstPeriodDC_ = dayCounter;
         return this;
      }

      public FixedRateLeg withLastPeriodDayCounter(DayCounter dayCounter)
      {
         lastPeriodDC_ = dayCounter;
         return this;
      }

      public FixedRateLeg withPaymentCalendar(Calendar cal)
      {
         calendar_ = cal;
         return this;
      }

      public FixedRateLeg withExCouponPeriod(Period period, Calendar cal, BusinessDayConvention convention, bool endOfMonth = false)
      {
         exCouponPeriod_ = period;
         exCouponCalendar_ = cal;
         exCouponAdjustment_ = convention;
         exCouponEndOfMonth_ = endOfMonth;
         return this;
      }

      // creator
      public override List<CashFlow> value()
      {

         if (couponRates_.Count == Const.ZERO_INT)
            throw new ArgumentException("no coupon rates given");
         if (notionals_.Count == Const.ZERO_INT)
            throw new ArgumentException("no nominals given");

         List<CashFlow> leg = new List<CashFlow>();

         Calendar schCalendar = schedule_.calendar();

         // first period might be short or long
         Date start = schedule_[Const.ZERO_INT], end = schedule_[Const.ONE_INT];
         Date paymentDate = calendar_.adjust(end, paymentAdjustment_);
         Date exCouponDate = null;
         InterestRate rate = couponRates_[Const.ZERO_INT];
         double nominal = notionals_[Const.ZERO_INT];

         if (exCouponPeriod_ != null)
         {
            exCouponDate = exCouponCalendar_.advance(paymentDate,
                                                     -exCouponPeriod_,
                                                     exCouponAdjustment_,
                                                     exCouponEndOfMonth_);
         }
         if (schedule_.isRegular(Const.ONE_INT))
         {
            if (!(firstPeriodDC_ == null || firstPeriodDC_ == rate.dayCounter()))
               throw new ArgumentException("regular first coupon does not allow a first-period day count");
            leg.Add(new FixedRateCoupon(paymentDate, nominal, rate, start, end, start, end, exCouponDate));
         }
         else
         {
            Date refer = end - schedule_.tenor();
            refer = schCalendar.adjust(refer, schedule_.businessDayConvention());
            InterestRate r = new InterestRate(rate.rate(),
                                              (firstPeriodDC_ == null || firstPeriodDC_.empty()) ? rate.dayCounter() : firstPeriodDC_,
                                              rate.compounding(), rate.frequency());
            leg.Add(new FixedRateCoupon(paymentDate, nominal, r, start, end, refer, end, exCouponDate));
         }

         // regular periods
         for (int i = Const.TWO_INT; i < schedule_.Count - Const.ONE_INT; ++i)
         {
            start = end; end = schedule_[i];
            paymentDate = calendar_.adjust(end, paymentAdjustment_);
            if (exCouponPeriod_ != null)
            {
               exCouponDate = exCouponCalendar_.advance(paymentDate,
                                                        -exCouponPeriod_,
                                                        exCouponAdjustment_,
                                                        exCouponEndOfMonth_);
            }
            if ((i - Const.ONE_INT) < couponRates_.Count)
               rate = couponRates_[i - Const.ONE_INT];
            else
               rate = couponRates_.Last();
            if ((i - Const.ONE_INT) < notionals_.Count)
               nominal = notionals_[i - Const.ONE_INT];
            else
               nominal = notionals_.Last();

            leg.Add(new FixedRateCoupon(paymentDate, nominal, rate, start, end, start, end, exCouponDate));
         }

         if (schedule_.Count > Const.TWO_INT)
         {
            // last period might be short or long
            int N = schedule_.Count;
            start = end; end = schedule_[N - Const.ONE_INT];
            paymentDate = calendar_.adjust(end, paymentAdjustment_);
            if (exCouponPeriod_ != null)
            {
               exCouponDate = exCouponCalendar_.advance(paymentDate,
                                                        -exCouponPeriod_,
                                                        exCouponAdjustment_,
                                                        exCouponEndOfMonth_);
            }

            if ((N - Const.TWO_INT) < couponRates_.Count)
               rate = couponRates_[N - Const.TWO_INT];
            else
               rate = couponRates_.Last();
            if ((N - Const.TWO_INT) < notionals_.Count)
               nominal = notionals_[N - Const.TWO_INT];
            else
               nominal = notionals_.Last();

            InterestRate r = new InterestRate(rate.rate(),
                                              lastPeriodDC_ == null ? rate.dayCounter() : lastPeriodDC_, rate.compounding(), rate.frequency());
            if (schedule_.isRegular(N - Const.ONE_INT))
               leg.Add(new FixedRateCoupon(paymentDate, nominal, r, start, end, start, end, exCouponDate));
            else
            {
               Date refer = start + schedule_.tenor();
               refer = schCalendar.adjust(refer, schedule_.businessDayConvention());
               leg.Add(new FixedRateCoupon(paymentDate, nominal, r, start, end, start, refer, exCouponDate));
            }
         }
         return leg;
      }
   }
}
