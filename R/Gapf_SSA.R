#' Gap-fill using SSA
#'
#' This function automatically gap-fills the missing data points (marked as "NA") in the flux dataset
#' using sigular spectrum analysis (SSA). The SSA is based on the algorithms in the package `spectral.methods` and `Rssa`.
#'
#' @param data a data frame that includes the flux (with NA indicating the missing data)
#' @param Flux a string indicates the column name for the flux variable to be gap-filled
#' @param Date a string indicates the column name for the date, which HAS to include the time information
#' @param Date_form a string indicates the format of the date, either "ymd_hms" (default), "mdy_hms" or "dmy_hms"
#' @param borders.wl list of numeric vectors indicate orders of the different periodicity bands to extract.
#' Units are the sampling frequency of the series
#' (needs one vector per step; default: list(a = c(0,10, Inf)) (see details in the package `spectral.methods`)).
#' @param win a number indicates the required sampling window length around each gap (total number in two sides), unit: days (default: 5)
#' @param interval a number indicates the temporal resolution of the measurements in the dataset, unit: minutes (default: 10)
#' @param fail a string or a number indicates what to do when model fails to converge:
#' 1. use the mean value in the sampling window to fill the gap ("ave", default), or
#' 2. use any value assigned here to fill the gap (e.g., 9999, NA, etc.)
#' @param ... other arguments pass to `gapfillSSA`
#' @return A data frame that includes the original data, gap-filled data ("filled")
#' and a "mark" column that indicates the value in each row of the "filled" is either:
#' 1. original, 2. gap-filled, or 3. failed to converge
#' @examples
#' # read example data
#' df <- read.csv(file = system.file("extdata", "Soil_resp_example.csv", package = "FluxGapsR"),header = T)
#' df_filled <- Gapfill_ssa(data = df)
#' # visualize the gapfilled results
#' plot(df_filled$filled,col="red")
#' points(df_filled$Flux)
#' @export
Gapfill_ssa <- function(data,
                        Flux = "Flux",
                        Date = "Date",
                        Date_form = "ymd_hms",
                        borders.wl = list(a = c(0,10, Inf)),
                        win = 5,
                        interval = 10,
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
                              stop("Invalid date format")))

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
    df_ssa <- data.frame(dft[wind_st:wind_ed,],
                         mk=mk[wind_st:wind_ed])

    # SSA fit
    dft_1 <- try(spectral.methods::gapfillSSA(series = df_ssa$Flux,print.stat=F),silent = TRUE) # use SSA to preliminarily fill the gaps
    # try different settings
    if (class(dft_1)=="try-error"){
      dft_1 <- spectral.methods::gapfillSSA(series = df_ssa$Flux,amnt.artgaps = c(0,0),size.biggap = 0,print.stat=F) # in case too many NAs around the gap
    }
    # error still, then claim a fail
    if (class(dft_1)=="try-error"){
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
    } else { # if low frequency successfully extracted
      # extract the high and low frequency
      # window length
      L <- ifelse(sum(!is.na(dft_1$filled.series))>1040,# if available data series is longer than 1000 + 40
                  1000, # use 1000 as window length
                  sum(!is.na(dft_1$filled.series))-40-1) # otherwise, use series length - n.comp (40) - 1
      data.decomposed <- try(spectral.methods::filterTSeriesSSA(series = dft_1$filled.series,
                                                                borders.wl = borders.wl, M = L,
                                                                # n.comp = c(50),
                                                                #harmonics = c(1, 0, 0),
                                                                plot.spectra = F, open.plot = F,print.stat=F),silent = T)

      if (class(data.decomposed)=="try-error"){ # if decomposition failed, change n.comp to 10
        data.decomposed <- try(spectral.methods::filterTSeriesSSA(series = dft_1$filled.series,
                                                                  borders.wl = borders.wl, M = L,
                                                                  n.comp = 10,
                                                                  #harmonics = c(1, 0, 0),
                                                                  plot.spectra = F, open.plot = F,print.stat=F),silent = T)
      }
      # if still fails, then use low frequency only
      if (class(data.decomposed)=="try-error"){ # if decomposition failed, change n.comp to 10
        gap[indx] <- dft_1$filled.series[df_ssa$mk==i]
        mark[indx] <- 1 # filled gap
        print(paste0("#",i," out of ",max(mk)," gaps: use low frequency only!!")) # for checking progress
      } else { # if decomposition succeed
        # combine data
    dft_2 <- data.frame(df_ssa,
                        filled=dft_1$filled.series,
                        high=data.decomposed$dec.series[1,],
                        low=data.decomposed$dec.series[2,]) %>%
      # add time info
      dplyr::mutate(DOY=lubridate::yday(Date),
             time=strftime(Date, format="%H:%M:%S") )
    dft_h <-
      dft_2 %>%  # use all the data in the sampling window for the high frequency
      dplyr::filter(mk!=i) %>% # drop those lines where gaps locate
      dplyr::group_by(time) %>%
      dplyr::summarise(H=mean(high,na.rm = T))
    # add the high frequency column to each day
    dft_2 <- dplyr::left_join(dft_2,dft_h,by="time")
    # calculate the gapfilled value and add to the "gap" vector
    gap[indx] <- rowSums(dft_2[which(dft_2$mk==i),c("filled","H")],na.rm = T) # in case high frequency is NAs
    mark[indx] <- 1 # filled gap
    print(paste0("#",i," out of ",max(mk)," gaps: succeed!!")) # for checking progress
      } # end of decomposition suceed
    } # end of low frequency extracted

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
