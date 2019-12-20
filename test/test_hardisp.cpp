#include <iostream>
#include "iers2010/iers2010.hpp"
#include "ggdatetime/dtcalendar.hpp"

/// A simple program to test the implementation of hardisp_impl
/// function, based on the test case provided in the respective FOTRAN
/// code.

int main()
{
  // set the date to performthe calculation
  // 2009 6 25 1 10 45
  auto epoch = ngpt::datetime<ngpt::seconds>{ngpt::year(2009), ngpt::month(6), 
      ngpt::day_of_month(25), ngpt::hours(1), 
      ngpt::minutes(10), ngpt::seconds(45)};


}
