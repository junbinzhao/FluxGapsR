#' Gap-fill using EM
#'
#' This function automatically gap-fills the missing data points (marked as "NA") in the flux dataset
#' using expectation-maximization (EM) algorithm with up to 3 parallel measured reference flux time series.
#' The function is based on the algorithms in the package `mtsdi`.
#'
#' @param data a data frame that includes the flux (with NA indicating the missing data)
#' @param ref1 a data frame that includes the parallel measured reference flux time series #1,
#' does not require to have the same length as the target data to be filled
#' @param ref2 a data frame that includes the parallel measured reference flux time series #2 (optional),
#' does not require to have the same length as the target data to be filled. Default: NULL
#' @param ref3 a data frame that includes the parallel measured reference flux time series #3 (optional),
#' does not require to have the same length as the target data to be filled. Default: NULL
#' @param Flux a string indicates the column name of the flux variable to be gap-filled
#' @param Flux1 a string indicates the column name of the reference time series in ref1. Default: same as Flux
#' @param Flux2 a string indicates the column name of the reference time series in ref2. Default: same as Flux
#' @param Flux3 a string indicates the column name of the reference time series in ref3. Default: same as Flux
#' @param Date a string indicates the column name for the date in data, ref1, ref2 and ref3,
#' and it HAS to include the time information. Note that all the data frames should have the same name for the date column.
#' @param Date_form a string indicates the format of the date in data, ref1, ref2 and ref3,
#' either "ymd_hms" (default), "mdy_hms" or "dmy_hms". Note that all the data frames should have the same date format.
#' @param win a number indicates the required sampling window length around each gap (total number in two sides), unit: days (default: 5)
#' @param interval a number indicates the temporal resolution of the measurements in the dataset, unit: minutes (default: 10)
#' @param ts logical. TRUE if it is time series. Default: TRUE
#' @param method a string indicates the method for univariate time series filtering,
#' either "spline" (default),"arima", or "gam". See details in the package `mtsdi`.
#' @param sp_df an integer indicates the degrees of freedom to be used for the splines (Default: 10).
#' In case set to NULL, the degrees of freedom will be chosen by cross-validation.
#' See details in the package `mtsdi`.
#' @param fail a string or a number indicates what to do when model fails to converge:
#' 1. use the mean value in the sampling window to fill the gap ("ave", default), or
#' 2. use any value assigned here to fill the gap (e.g., 9999, NA, etc.)
#' @param ... other arguments pass to `mnimput`
#' @return A data frame that includes the original data, gap-filled data ("filled")
#' and a "mark" column that indicates the value in each row of the "filled" is either:
#' 1. original, 2. gap-filled, or 3. failed to converge
#' @examples
#' # read example data
#' df <- read.csv(file = system.file("extdata", "Soil_resp_example.csv", package = "FluxGapsR"),header = T)
#' df_ref <- read.csv(file = system.file("extdata", "Soil_resp_ref_example.csv", package = "FluxGapsR"),header = T)
#' df_filled <- Gapfill_em(data = df,ref1 = df_ref)
#' # visualize the gapfilled results
#' plot(df_filled$filled,col="red")
#' points(df_filled$Flux)
#' @export
Gapfill_em <- function(data,
                       ref1,
                       ref2 = NULL,
                       ref3 = NULL,
                       Flux = "Flux",
                       Flux1 = Flux,
                       Flux2 = Flux,
                       Flux3 = Flux,
                       Date = "Date",
                       Date_form = "ymd_hms",
                       win = 5,
                       interval = 10,
                       ts = TRUE,
                       method = "spline",
                       sp_df = 10,
                       fail = "ave",
                       ...
){
  # # define the pipe from the package "magrittr"
  `%>%` <- magrittr::`%>%`
  ### add sequence mark to the gaps -------
  mt <- is.na(data[,Flux])
  ind <- 1 # index for marking the gaps
  mk <- vector()
  for (i in 1:length(mt)) {
    if (mt[i]==FALSE){
      mk[i] <- 0 # non-gaps are marked as 0
    } else {
      if (mt[i]==TRUE){
        mk[i] <- ind # gaps are marked as the value of ind
        if (mt[i+1]==FALSE) {
          ind <- ind+1 # when reached the end of a gap, change add 1 to ind
        }
      }
    }
  }
  print(paste0(max(mk)," gaps are marked")) # display the total number of gaps

  ### prepare data for gapfilling -----
  # the sampling window length
  pt_h <- 60/interval # how many data points per hour
  winID <- win/2*pt_h*24 # how many data points for the sampling window at EACH side of the gap
  # create vector to save the predicted gapfilled data
  gap <- rep(NA,nrow(data))
  # extract the data needed for gap-filling
  dft <- data[,c(Date,Flux)]
  names(dft) <- c("Date","Flux")
  dft <- dft %>%
    dplyr::mutate(Date=switch(Date_form, # convert the date into right format depending on the input format
                              "ymd_hms" = lubridate::ymd_hms(Date),
                              "mdy_hms" = lubridate::mdy_hms(Date),
                              "dmy_hms" = lubridate::dmy_hms(Date),
                              stop("Invalid date format for data")))
  # extract the data for reference1
  ref1 <- ref1[,c(Date,Flux1)]
  names(ref1) <- c("Date","Flux1")
  ref1 <- ref1 %>%
    dplyr::mutate(Date=switch(Date_form, # convert the date into right format depending on the input format
                              "ymd_hms" = lubridate::ymd_hms(Date),
                              "mdy_hms" = lubridate::mdy_hms(Date),
                              "dmy_hms" = lubridate::dmy_hms(Date),
                              stop("Invalid date format for ref1")))
  dft <- dplyr::left_join(dft,ref1,by = "Date")
  formula <- as.formula("~Flux+Flux1")
  # if the second reference is present
  if (!is.null(ref2)){
    ref2 <- ref2[,c(Date,Flux2)]
    names(ref2) <- c("Date","Flux2")
    ref2 <- ref2 %>%
      dplyr::mutate(Date=switch(Date_form, # convert the date into right format depending on the input format
                                "ymd_hms" = lubridate::ymd_hms(Date),
                                "mdy_hms" = lubridate::mdy_hms(Date),
                                "dmy_hms" = lubridate::dmy_hms(Date),
                                stop("Invalid date format for ref2")))
    dft <- dplyr::left_join(dft,ref2,by = "Date")
    formula <- as.formula("~Flux+Flux1+Flux2")
  }
  # if the third reference is present
  if (!is.null(ref3)){
    ref3 <- ref3[,c(Date,Flux3)]
    names(ref3) <- c("Date","Flux3")
    ref3 <- ref3 %>%
      dplyr::mutate(Date=switch(Date_form, # convert the date into right format depending on the input format
                                "ymd_hms" = lubridate::ymd_hms(Date),
                                "mdy_hms" = lubridate::mdy_hms(Date),
                                "dmy_hms" = lubridate::dmy_hms(Date),
                                stop("Invalid date format for ref3")))
    dft <- dplyr::left_join(dft,ref3,by = "Date")
    formula <- as.formula("~Flux+Flux1+Flux2+Flux3")
  }

  # a vector for marks of each gap
  mark <- rep(0,nrow(dft))
  # a number to record the number of failed regression
  nf <- 0

  ### gap filling by the marked index of each gap ----------
  for (i in 1:max(mk)) {
    indx <- which(mk==i) # index of the gap
    # define the sampling window
    wind_st <- ifelse(min(indx)-winID>=0,min(indx)-winID,1) # use the beginning of time series if not enough sample points are present
    wind_ed <- ifelse(max(indx)+winID>nrow(data),nrow(data),max(indx)+winID) # use the end if not enough
    # extract data to fit the model
    df_em <- data.frame(dft[wind_st:wind_ed,],
                        mk=mk[wind_st:wind_ed])

    # EM imputation
    if (method=="spline"){ # use spline as default
      df_out <- try(mtsdi::mnimput(formula = formula,
                      df_em,
                      ts=ts, method="spline",
                      sp.control=list(df=sp_df),
                      ...
                      ),
              silent = T)
    } else { # or use other methods
      df_out <- try(mtsdi::mnimput(formula = formula,
                                   df_em,
                                   ts=ts, method=method,
                                   ...
      ),
      silent = T)
    }

    # if error
    if (class(df_out)=="try-error"){
      if (fail == "ave"){ # use average in the sampling window
        gap[indx] <- mean(dft$Flux[wind_st:wind_ed],na.rm = T)
        mark[indx] <- 2 # failed to filled gap
        nf <- nf+1 # add up the failed times
        print(paste0("#",i," out of ",max(mk)," gaps: Failed...")) # for checking progress
      } else { # or use the designated value
        gap[indx] <- fail
        mark[indx] <- 2 # failed to filled gap
        nf <- nf+1 # add up the failed times
        print(paste0("#",i," out of ",max(mk)," gaps: Failed...")) # for checking progress
      }
    } else { # if imputation succeed
      gap[indx] <- df_out$filled.dataset$Flux[df_em$mk==i]
      mark[indx] <- 1 # filled gap
      print(paste0("#",i," out of ",max(mk)," gaps: succeed!!")) # for checking progress
    }
  } # end of the loop
  df_new <- data.frame(data,
                       filled = gap,
                       tem = dft[,"Flux"],
                       mark) %>%
    dplyr::mutate(filled = ifelse(mark==0,tem,filled)) %>%
    dplyr::select(-tem) # drop the temperory column

  # print a summary of the gapfilling ------------
  stat <- table(mk)[-1] # number of data points in each gap
  # print using "cat" for break into lines
  cat(paste0("","\n",
             "##### Summary #####","\n",
             "","\n",
             "Total gaps:       ",max(mk),"\n",
             "< 1 day:          ",sum(stat<pt_h*24),"\n",
             ">= 1 & < 7 days:  ",sum(stat>=pt_h*24 & stat<pt_h*24*7),"\n",
             ">= 7 & < 15 days: ",sum(stat>=pt_h*24*7 & stat<pt_h*24*15),"\n",
             ">= 15 days:       ",sum(stat>=pt_h*24*15),"\n",
             "Failed gaps:      ",nf
  ))
  # return the output data frame
  return(df_new)
}
