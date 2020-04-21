/*
 Copyright (C) 2008, 2009 , 2010 Andrea Maggiulli (a.maggiulli@gmail.com)

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
using System.Collections.Generic;

namespace QLNet
{
   public abstract class PrincipalLegBase
   {
      protected Schedule schedule_;
      protected List<double> notionals_;
      protected BusinessDayConvention paymentAdjustment_;
      protected DayCounter dayCounter_;
      protected int sign_;
      public static implicit operator List<CashFlow>(PrincipalLegBase o) { return o.value(); }
      public abstract List<CashFlow> value();


      // initializers
      public PrincipalLegBase withNotionals(double notional)
      {
         notionals_ = new List<double>() { notional };
         return this;
      }
      public PrincipalLegBase withNotionals(List<double> notionals)
      {
         notionals_ = notionals;
         return this;
      }
      public PrincipalLegBase withPaymentAdjustment(BusinessDayConvention convention)
      {
         paymentAdjustment_ = convention;
         return this;
      }
      public PrincipalLegBase withPaymentDayCounter(DayCounter dayCounter)
      {
         dayCounter_ = dayCounter;
         return this;
      }

      public PrincipalLegBase withSign(int sign)
      {
         sign_ = sign;
         return this;
      }

   }

   //! helper class building a Bullet Principal leg
   public class PricipalLeg : PrincipalLegBase
   {
      // constructor
      public PricipalLeg(Schedule schedule, DayCounter paymentDayCounter)
      {
         schedule_ = schedule;
         paymentAdjustment_ = BusinessDayConvention.Following;
         dayCounter_ = paymentDayCounter;
      }

      // creator
      public override List<CashFlow> value()
      {
         List<CashFlow> leg = new List<CashFlow>();

         // the following is not always correct
         Calendar calendar = schedule_.calendar();

         // first period
         Date start = schedule_[Const.ZERO_INT], end = schedule_[schedule_.Count - Const.ONE_INT];
         Date paymentDate = calendar.adjust(start, paymentAdjustment_);
         double nominal = notionals_[Const.ZERO_INT];
         double quota = nominal / (schedule_.Count - Const.ONE_INT);

         leg.Add(new Principal(nominal * sign_, nominal, paymentDate, start, end, dayCounter_,  start, end));

         if (schedule_.Count == Const.TWO_INT)
         {
            paymentDate = calendar.adjust(end, paymentAdjustment_);
            leg.Add(new Principal(nominal * sign_ * Const.NEGATIVE_ONE_DOUBLE, Const.ZERO_INT, paymentDate, start, end, dayCounter_, start, end));
         }
         else
         {
            end = schedule_[Const.ZERO_INT];
            // regular periods
            for (int i = Const.ONE_INT; i <= schedule_.Count - Const.ONE_INT; ++i)
            {
               start = end; end = schedule_[i];
               paymentDate = calendar.adjust(start, paymentAdjustment_);
               nominal -= quota;

               leg.Add(new Principal(quota * sign_ * Const.NEGATIVE_ONE_DOUBLE, nominal, paymentDate, start, end, dayCounter_, start, end));
            }
         }

         return leg;
      }
   }

   //! helper class building a Bullet Principal leg
   public class BulletPricipalLeg : PrincipalLegBase
   {
      // constructor
      public BulletPricipalLeg(Schedule schedule)
      {
         schedule_ = schedule;
         paymentAdjustment_ = BusinessDayConvention.Following;
      }

      // creator
      public override List<CashFlow> value()
      {
         List<CashFlow> leg = new List<CashFlow>();

         // the following is not always correct
         Calendar calendar = schedule_.calendar();

         // first period might be short or long
         Date start = schedule_[Const.ZERO_INT], end = schedule_[Const.ONE_INT];
         Date paymentDate = calendar.adjust(start, paymentAdjustment_);
         double nominal = notionals_[Const.ZERO_INT];

         leg.Add(new Principal(nominal, nominal, paymentDate, start, end, dayCounter_, start, end));

         paymentDate = calendar.adjust(end, paymentAdjustment_);
         leg.Add(new Principal(nominal * Const.NEGATIVE_ONE_DOUBLE, Const.ZERO_INT, paymentDate, start, end, dayCounter_, start, end));

         return leg;
      }
   }
}
