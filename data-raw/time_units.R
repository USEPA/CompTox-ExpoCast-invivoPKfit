# list of allowable time units:
# these are the units understood by lubridate::duration()

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
