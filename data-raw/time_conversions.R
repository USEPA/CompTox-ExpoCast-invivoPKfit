# Time table data frame for quicker access to time conversions
# This was made in light of recent changes to lubridate::duration()
# that rendered the convert_time() function unstable
# Definitions:
# picoseconds - seconds (classic metric definition)
# minute = 60 seconds
# hour = 60 minutes
# day = 24 hours
# week = 7 days
# month = 1/12 of a year
# year = 365 days

time_conversions <- data.frame(
  TimeFrom = c(
    "picoseconds", "picoseconds", "picoseconds", "picoseconds", "picoseconds", "picoseconds",
    "picoseconds", "picoseconds", "picoseconds", "picoseconds", "picoseconds",
    "nanoseconds", "nanoseconds", "nanoseconds", "nanoseconds", "nanoseconds", "nanoseconds",
    "nanoseconds", "nanoseconds", "nanoseconds", "nanoseconds", "nanoseconds",
    "microseconds", "microseconds", "microseconds", "microseconds", "microseconds", "microseconds",
    "microseconds", "microseconds", "microseconds", "microseconds", "microseconds",
    "milliseconds", "milliseconds", "milliseconds", "milliseconds", "milliseconds", "milliseconds",
    "milliseconds", "milliseconds", "milliseconds", "milliseconds", "milliseconds",
    "seconds", "seconds", "seconds", "seconds", "seconds", "seconds",
    "seconds", "seconds", "seconds", "seconds", "seconds",
    "minutes", "minutes", "minutes", "minutes", "minutes", "minutes",
    "minutes", "minutes", "minutes", "minutes", "minutes",
    "hours", "hours", "hours", "hours", "hours", "hours",
    "hours", "hours", "hours", "hours", "hours",
    "days", "days", "days", "days", "days", "days",
    "days", "days", "days", "days", "days",
    "weeks", "weeks", "weeks", "weeks", "weeks", "weeks",
    "weeks", "weeks", "weeks", "weeks", "weeks",
    "months", "months", "months", "months", "months", "months",
    "months", "months", "months", "months", "months",
    "years", "years", "years", "years", "years", "years",
    "years", "years", "years", "years", "years"
  ),
  TimeTo = c(
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years",
    "picoseconds", "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes",
    "hours", "days", "weeks", "months", "years"
    ),
  conversion = c(
    # from picoseconds
    1, 1/(1E3), 1/(1E6), 1/(1E9), 1/(1E12), 1/(60*1E12),
    1/(60*60*1E12), 1/(24*60*60*1E12), 1/(7*24*60*60*1E12), 12/(365*24*60*60*1E12), 1/(365*24*60*60*1E12),
    # from nanoseconds
    1E3, 1, 1/(1E3), 1/(1E6), 1/(1E9), 1/(60*1E9),
    1/(60*60*1E9), 1/(24*60*60*1E9), 1/(7*24*60*60*1E9), 12/(365*24*60*60*1E9), 1/(365*24*60*60*1E9),
    # from microseconds
    1E6, 1E3, 1, 1/(1E3), 1/(1E6), 1/(60*1E6),
    1/(60*60*1E6), 1/(24*60*60*1E6), 1/(7*24*60*60*1E6), 12/(365*24*60*60*1E6), 1/(365*24*60*60*1E6),
    # from milliseconds
    1E9, 1E6, 1E3, 1, 1/(1E3), 1/(60*1E3),
    1/(60*60*1E3), 1/(24*60*60*1E3), 1/(7*24*60*60*1E3), 12/(365*24*60*60*1E3), 1/(365*24*60*60*1E3),
    # from seconds
    1E12, 1E9, 1E6, 1E3, 1, 1/(60),
    1/(60*60), 1/(24*60*60), 1/(7*24*60*60), 12/(365*24*60*60), 1/(365*24*60*60),
    # from minutes
    60*1E12, 60*1E9, 60*1E6, 60*1E3, 60, 1,
    1/60, 1/(24*60), 1/(7*24*60), 12/(365*24*60), 1/(365*24*60),
    # from hours
    60*60*1E12, 60*60*1E9, 60*60*1E6, 60*60*1E3, 60*60, 60,
    1, 1/24, 1/(7*24), 12/(365*24), 1/(365*24),
    # from days
    24*60*60*1E12, 24*60*60*1E9, 24*60*60*1E6, 24*60*60*1E3, 24*60*60, 24*60,
    24, 1, 1/7, 12/365, 1/365,
    # from weeks
    7*24*60*60*1E12, 7*24*60*60*1E9, 7*24*60*60*1E6, 7*24*60*60*1E3, 7*24*60*60, 7*24*60,
    7*24, 7, 1, (12*7)/365, 7/365,
    # from months
    (365/12)*24*60*60*1E12, (365/12)*24*60*60*1E9, (365/12)*24*60*60*1E6, (365/12)*24*60*60*1E3, (365/12)*24*60*60, (365/12)*24*60,
    (365/12)*24, 365/12, 365/(12*7), 1, 1/12,
    # from years
    365*24*60*60*1E12, 365*24*60*60*1E9, 365*24*60*60*1E6, 365*24*60*60*1E3, 365*24*60*60, 365*24*60,
    365*24, 365, 365/7, 12, 1
    )

)

usethis::use_data(time_conversions, overwrite = TRUE)
