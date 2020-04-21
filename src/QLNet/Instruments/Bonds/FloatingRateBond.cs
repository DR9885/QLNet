/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)

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

namespace QLNet
{
   //! floating-rate bond (possibly capped and/or floored)
   //! \test calculations are tested by checking results against cached values.
   public class FloatingRateBond : Bond
   {
      public FloatingRateBond(int settlementDays, double faceAmount, Schedule schedule, IborIndex index,
                              DayCounter paymentDayCounter)
         : this(settlementDays, faceAmount, schedule, index, paymentDayCounter, BusinessDayConvention.Following,
                Const.ZERO_INT, new List<double>() { Const.ONE_INT }, new List<double>() { Const.ZERO_INT }, new List < double? >(), new List < double? >(),
      false, Const.ONE_HUNDRED_INT, null) { }
      public FloatingRateBond(int settlementDays, double faceAmount, Schedule schedule, IborIndex index,
                              DayCounter paymentDayCounter, BusinessDayConvention paymentConvention, int fixingDays,
                              List<double> gearings, List<double> spreads)
         : this(settlementDays, faceAmount, schedule, index, paymentDayCounter, BusinessDayConvention.Following,
                fixingDays, gearings, spreads, new List < double? >(), new List < double? >(), false, Const.ONE_HUNDRED_INT, null) { }
      public FloatingRateBond(int settlementDays, double faceAmount, Schedule schedule, IborIndex index, DayCounter paymentDayCounter,
                              BusinessDayConvention paymentConvention, int fixingDays, List<double> gearings, List<double> spreads,
                              List < double? > caps, List < double? > floors, bool inArrears, double redemption, Date issueDate)
      : base(settlementDays, schedule.calendar(), issueDate)
      {
         maturityDate_ = schedule.endDate();
         cashflows_ = new IborLeg(schedule, index)
         .withPaymentDayCounter(paymentDayCounter)
         .withFixingDays(fixingDays)
         .withGearings(gearings)
         .withSpreads(spreads)
         .withCaps(caps)
         .withFloors(floors)
         .inArrears(inArrears)
         .withNotionals(faceAmount)
         .withPaymentAdjustment(paymentConvention);

         addRedemptionsToCashflows(new List<double>() { redemption });

         Utils.QL_REQUIRE(cashflows().Count != Const.ZERO_INT, () => "bond with no cashflows!");
         Utils.QL_REQUIRE(redemptions_.Count == Const.ONE_INT, () => "multiple redemptions created");

         index.registerWith(update);
      }

      //public FloatingRateBond(int settlementDays, double faceAmount, Date startDate, Date maturityDate, Frequency couponFrequency,
      //                        Calendar calendar, IborIndex index, DayCounter accrualDayCounter,
      //                        BusinessDayConvention accrualConvention = Following,
      //                        BusinessDayConvention paymentConvention = Following,
      //                        int fixingDays = Null<Natural>(),
      //                        List<double> gearings = std::vector<Real>(1, 1.0),
      //                        List<double> spreads = std::vector<Spread>(1, 0.0),
      //                        List<double> caps = std::vector<Rate>(),
      //                        List<double> floors = std::vector<Rate>(),
      //                        bool inArrears = false,
      //                        double redemption = 100.0,
      //                        Date issueDate = Date(),
      //                        Date stubDate = Date(),
      //                        DateGeneration.Rule rule = DateGeneration::Backward,
      //                        bool endOfMonth = false)
      public FloatingRateBond(int settlementDays, double faceAmount, Date startDate, Date maturityDate, Frequency couponFrequency,
                              Calendar calendar, IborIndex index, DayCounter accrualDayCounter,
                              BusinessDayConvention accrualConvention, BusinessDayConvention paymentConvention,
                              int fixingDays, List<double> gearings, List<double> spreads, List < double? > caps,
                              List < double? > floors, bool inArrears, double redemption, Date issueDate,
                              Date stubDate, DateGeneration.Rule rule, bool endOfMonth)
      : base(settlementDays, calendar, issueDate)
      {

         maturityDate_ = maturityDate;

         Date firstDate = null, nextToLastDate = null;
         switch (rule)
         {
            case DateGeneration.Rule.Backward:
               firstDate = null;
               nextToLastDate = stubDate;
               break;
            case DateGeneration.Rule.Forward:
               firstDate = stubDate;
               nextToLastDate = null;
               break;
            case DateGeneration.Rule.Zero:
            case DateGeneration.Rule.ThirdWednesday:
            case DateGeneration.Rule.Twentieth:
            case DateGeneration.Rule.TwentiethIMM:
               Utils.QL_FAIL("stub date (" + stubDate + ") not allowed with " + rule + " DateGeneration::Rule");
               break;
            default:
               Utils.QL_FAIL("unknown DateGeneration::Rule (" + rule + ")");
               break;
         }

         Schedule schedule = new Schedule(startDate, maturityDate_, new Period(couponFrequency), calendar_,
                                          accrualConvention, accrualConvention, rule, endOfMonth, firstDate, nextToLastDate);

         cashflows_ = new IborLeg(schedule, index)
         .withPaymentDayCounter(accrualDayCounter)
         .withFixingDays(fixingDays)
         .withGearings(gearings)
         .withSpreads(spreads)
         .withCaps(caps)
         .withFloors(floors)
         .inArrears(inArrears)
         .withNotionals(faceAmount)
         .withPaymentAdjustment(paymentConvention);

         addRedemptionsToCashflows(new List<double>() { redemption });

         Utils.QL_REQUIRE(cashflows().Count != Const.ZERO_INT, () => "bond with no cashflows!");
         Utils.QL_REQUIRE(redemptions_.Count == Const.ONE_INT, () => "multiple redemptions created");

         index.registerWith(update);
      }

   }
}
