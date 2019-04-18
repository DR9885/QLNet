/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2008 Toyin Akin (toyin_akin@hotmail.com)
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

namespace QLNet
{
   // we need only one instance of the class
   // we can not derive it from IObservable because the class is static
   public class Settings : Singleton<Settings>
   {

      public Settings()
      {}
      private ObservableValue<Date> evaluationDate_;
      private bool includeReferenceDateEvents_;
      private bool enforcesTodaysHistoricFixings_;
      private bool? includeTodaysCashFlows_;

      public Date evaluationDate()
      {
         if (evaluationDate_ == null)
            evaluationDate_ = Date.Today;
         return evaluationDate_;
      }


      public void setEvaluationDate(Date d)
      {
         evaluationDate_ = d;
         notifyObservers();
      }

      public bool enforcesTodaysHistoricFixings
      {
         get
         {
            return enforcesTodaysHistoricFixings_;
         }
         set
         {
            enforcesTodaysHistoricFixings_ = value;
         }
      }

      public bool includeReferenceDateEvents
      {
         get
         {
            return includeReferenceDateEvents_;
         }
         set
         {
            includeReferenceDateEvents_ = value;
         }
      }

      public bool? includeTodaysCashFlows
      {
         get
         {
            return includeTodaysCashFlows_;
         }
         set
         {
            includeTodaysCashFlows_ = value;
         }
      }

      ////////////////////////////////////////////////////
      // Observable interface
      private static readonly WeakEventSource eventSource = new WeakEventSource();
      public static event Callback notifyObserversEvent
      {
         add
         {
            eventSource.Subscribe(value);
         }
         remove
         {
            eventSource.Unsubscribe(value);
         }
      }

      public static void registerWith(Callback handler) { notifyObserversEvent += handler; }
      public static void unregisterWith(Callback handler) { notifyObserversEvent -= handler; }
      private static void notifyObservers()
      {
         eventSource.Raise();
      }

      public static void clearObservers()
      {
         eventSource.Clear();
      }
   }

   // helper class to temporarily and safely change the settings
   public class SavedSettings : IDisposable
   {
      private Date evaluationDate_;
      private bool enforcesTodaysHistoricFixings_;
      private bool includeReferenceDateEvents_;
      private bool? includeTodaysCashFlows_;

      public SavedSettings()
      {
         evaluationDate_ = Settings.Instance.evaluationDate();
         enforcesTodaysHistoricFixings_ = Settings.Instance.enforcesTodaysHistoricFixings;
         includeReferenceDateEvents_ = Settings.Instance.includeReferenceDateEvents;
         includeTodaysCashFlows_ = Settings.Instance.includeTodaysCashFlows;
      }

      public void Dispose()
      {
         if (evaluationDate_ != Settings.Instance.evaluationDate())
            Settings.Instance.setEvaluationDate(evaluationDate_);
         Settings.Instance.enforcesTodaysHistoricFixings = enforcesTodaysHistoricFixings_;
         Settings.Instance.includeReferenceDateEvents = includeReferenceDateEvents_;
         Settings.Instance.includeTodaysCashFlows = includeTodaysCashFlows_;
      }
   }
}
