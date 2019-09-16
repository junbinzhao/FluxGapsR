#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param data a data frame that includes the flux and soil temperature
#' @param Flux a string indicates the column name for the flux variable, such as soil respiration rate
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
  winID <- win/2*6*24 # calculate how many data points need to be sample at each side of the gap
  gap <- rep(NA,nrow(data)) # create vector to save the predicted gapfilled data
  # F5h_nls3 <- rep(NA,nrow(plot1_gap))
  # q <- enquo(Gap_type) # quote the column name for dplyr function
  # qs <- quo_name(q) # convert the column name into a string
  for (i in 1:max(data[,Gap_type])) {
    dft <- data %>% # create a dataset from the origial one for gapfilling with data in the artificial gap removed
      # mutate(Flux=ifelse(!!q==i,NA,Flux)) %>% # remove the data in the gap
      mutate(Flux=ifelse(get(Gap_type)==i,NA,Flux)) %>% # remove the data in the gap
      select(Flux,Ts,Ta)
    indx <- which(data[,Gap_type]==i) # index of the gap
    #################### use static time window
    # define the window
    wind_st <- ifelse(min(indx)-winID>=0,min(indx)-winID,1) # use the beginning of time series if not enough sample points are present
    wind_ed <- ifelse(max(indx)+winID>nrow(data),nrow(data),max(indx)+winID) # use the end if not enough
    fit1 <- try(nlsLM(Flux~a*exp(b*Ts), # a basic model
                      start = list(a=1,b=0.1),
                      data = dft[wind_st:wind_ed,],# use the data within 30 days window
                      control=nls.lm.control(maxiter = 1000)),
                silent = T)
    fit2 <- try(nlsLM(Flux~a*exp(b/(227.13-(Ts+273.15))), # Lloyd-Taylor model
                      start = list(a=1000,b=400),
                      data = dft[wind_st:wind_ed,],# use the data within 30 days window
                      control=nls.lm.control(maxiter = 1000)),
                silent = T)
    fit3 <- try(lm(Flux~poly(Ts,2),
                   data = dft[wind_st:wind_ed,]),
                silent = T)
    fit4 <- try(lm(Flux~Ts,
                   data = dft[wind_st:wind_ed,]),
                silent = T)

    ######################
    # for (j in indx) { # fill each point based on the index
    # fit1 <- try(nlsLM(Flux~a*exp(b*Ts), # a basic model
    #             start = list(a=1,b=0.1),
    #             data = dft[(j-2160):(j+2160),]),# use the data within 30 days window
    #             silent = T)
    # fit2 <- try(nlsLM(Flux~a*exp(b/(227.13-(Ts+273.15))), # Lloyd-Taylor model
    #             start = list(a=1000,b=400),
    #             data = dft[(j-2160):(j+2160),]),# use the data within 30 days window
    #             silent = T)
    # fit3 <- lm(Flux~poly(Ts,2),
    #            data = dft[(j-2160):(j+2160),])
    # predict the flux based on the soil temperature for gapfilling
    if (class(fit2)!="try-error"){ # if fit2 is not error, use it
      gap[indx] <- predict(fit2,newdata=dft[indx,])
    } else { # but if fit2 is error, then ...
      if (class(fit1)!="try-error"){ # then use fit1
        gap[indx] <- predict(fit1,newdata=dft[indx,])
      } else { # but if fit1 is also error, then ........
        if (class(fit3)!="try-error"){ # then use fit3
          gap[indx] <- predict(fit3,newdata=dft[indx,])
        } else { # but if fit3 is also error, then use fit4
          gap[indx] <- predict(fit4,newdata=dft[indx,])
        }
      }
    }
    #     if (class(fit1)=="try-error"){
    #       if (class(fit3)=="try-error"){
    #         gap[j] <- predict(fit4,newdata=dft[j,])# if fit1-3 are error, use fit4
    #       } else {
    #         gap[j] <- predict(fit3,newdata=dft[j,]) # if both fit1 and 2 are error, use fit3
    #         } else {
    #           gap[j] <- predict(fit1,newdata=dft[j,]) # if fit2 is error and fit1 is not, use fit1
    #           }
    #     } else {
    #     if (class(fit1)=="try-error"){
    #       gap[j] <- predict(fit2,newdata=dft[j,]) # if fit1 is error but not fit2, use fit2
    #     } else {
    #       if (sum(abs(resid(fit1)))<sum(abs(resid(fit2)))){# if both fit 1 and 2 are NOT error, use the one with smaller residuals
    #       gap[j] <- predict(fit1,newdata=dft[j,])
    #     } else{
    #       gap[j] <- predict(fit2,newdata=dft[j,])
    #     }
    #     }
    #
    #   }
    # }
    # }
    print(c(i,max(data[,Gap_type]),Gap_type)) # for checking progress
  }
  df_new <- data.frame(data,gap,mark)
  return(df_new) # output the vector with all the gapfilled data
}
