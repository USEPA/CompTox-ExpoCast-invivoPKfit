# list of allowable time units:
# these are the units understood by lubridate::duration()
# Only seconds minutes hours days weeks are now understood by lubridate::duration()
# See convert_time_table
time_units <- c("picoseconds",
                "nanoseconds",
                "microseconds",
                "milliseconds",
                "seconds",
                "minutes",
                "hours",
                "days",
                "weeks",
                "months",
                "years")

usethis::use_data(time_units, overwrite = TRUE)
