#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param data a data frame that includes the flux (with NA indicating the missing data) and soil temperature
#' @param Flux a string indicates the column name for the flux variable to be gap-filled
#' @param Ts a string indicates the column name for the soil temperature
#' @param win a number indicates total sampling window length around each gap, unit: days (default: 5)
#' @param interval a number indicates the temporal resolution of the dataset, unit: minutes (default: 10)
#' @param R10 the start value for the parameter R10 in the regression (default: 1000)
#' @param E0 the start value for the parameter E0 in the regression (default: 400)
#' @param fail a string or a number, what to do when model fails to converge:
#' 1. use the mean value in the sampling window to fill the gap ("ave",default), or
#' 2. use any value assigned here to fill the gap (e.g., 9999, NA, etc.)
#' @return A data frame that includes the original data, gap-filled data ("Filled")
#' and a "Mark" column that indicates the value in each row of "Filled" is either:
#' 1. original, 2. gap-filled, or 3. failed to converge
#' @examples
#' df_filled <- Gapfill_nls(data = df)
#' @export
Gapfill_nls <- function(data,
                        Flux = "Flux",
                        Ts = "Ts",
                        win = 5,
                        interval = 10,
                        R10 = 1000,
                        E0 = 400,
                        fail = "ave"
                        ){
  # # define the pipe from the package "magrittr"
  # `%>%` <- magrittr::`%>%`
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
  dft <- data[,c(Flux,Ts)]
  names(dft) <- c("Flux","Ts")
  # a vector for marks of each gap
  mark <- rep(0,length(dft))
  # a number to record the number of failed regression
  nf <- 0

  ### gap filling by the marked index of each gap ----------
  for (i in 1:max(mk)) {
    indx <- which(mk==i) # index of the gap
    # define the sampling window
    wind_st <- ifelse(min(indx)-winID>=0,min(indx)-winID,1) # use the beginning of time series if not enough sample points are present
    wind_ed <- ifelse(max(indx)+winID>nrow(data),nrow(data),max(indx)+winID) # use the end if not enough

    # fit the nls model
    fit <- try(minpack.lm::nlsLM(Flux~a*exp(b*(1/(283.15-227.13)-1/(Ts-227.13))),
                                 start = list(a=R10,b=E0),
                                 data = dft[wind_st:wind_ed,],
                                 control=nls.lm.control(maxiter = 1000)
                                 ),
               silent = TRUE)
    if (class(fit)=="try-error"){# it failed to fit the Lloyd-Taylor model try the basic model
      fit <- try(minpack.lm::nlsLM(Flux~a*exp(b*Ts),
                                   start = list(a=1,b=0.1),
                                   data = dft[wind_st:wind_ed,],
                                   control=nls.lm.control(maxiter = 1000)
      ),
      silent = TRUE)
    }

    # predict the gaps
    if (class(fit)!="try-error"){ # if the fit converged
      gap[indx] <- predict(fit,newdata=dft[indx,])
      mark[indx] <- 1 # filled gap
      print(paste0("#",i," out of ",max(mk)," gaps: succeed!!")) # for checking progress
    } else {
      if (fail == "ave"){ # use average in the sampling window
        gap[indx] <- mean(dft[wind_st:wind_ed,"Flux"],na.rm = T)
        mark[indx] <- 2 # failed to filled gap
        print(paste0("#",i," out of ",max(mk)," gaps: Failed...")) # for checking progress
      } else { # or use the designated value
        gap[indx] <- fail
        mark[indx] <- 2 # failed to filled gap
        nf <- nf+1 # add up the failed times
        print(paste0("#",i," out of ",max(mk)," gaps: Failed...")) # for checking progress
      }
    }
  } # end of the loop
  df_new <- data.frame(data,gap,mark)
  return(df_new)
  # print a summary of the gapfilling ------------
  stat <- table(mk)[-1] # number of data points in each gap
  # print using "cat" for break into lines
  cat(paste0("Total gaps:       ",max(mk),"\n",
             "< 1 day:          ",sum(stat<pt_h*24),"\n",
             ">= 1 & < 7 days:  ",sum(stat>=pt_h*24 & stat<pt_h*24*7),"\n",
             ">= 7 & < 15 days: ",sum(stat>=pt_h*24*7 & stat<pt_h*24*15),"\n",
             ">= 15 days:       ",sum(stat>=pt_h*24*15),"\n",
             "Failed gaps:      ",nf
             ))
}
