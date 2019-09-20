#' Gap-fill using ANN
#'
#' This function automatically gap-fills the missing data points (marked as "NA") in the flux dataset
#' using artificial neural networks (ANN) that take up to three variables as inputs. The ANN algorithms are based on
#' the package `neuralnet`.
#'
#' @param data a data frame that includes the flux (with NA indicating the missing data) and independent variables
#' @param Flux a string indicates the column name for the flux variable to be gap-filled
#' @param var1 a string indicates the column name for the first variable
#' @param var2 a string indicates the column name for the second variable, default: NULL
#' @param var3 a string indicates the column name for the third variable, default: NULL
#' @param win a number indicates the required sampling window length around each gap (total number in two sides), unit: days (default: 5)
#' @param interval a number indicates the temporal resolution of the measurements in the dataset, unit: minutes (default: 10)
#' @param threshold a number specifies the threshold for the partial derivatives of the error function as stopping criteria for the ANN model (default: 1)
#' @param hidden a vector of integers specifies the number of hidden neurons (vertices) in each layer in the ANN model (default: c(2), i.e. one layer with 2 neurons)
#' @param fail a string or a number indicates what to do when model fails to converge:
#' 1. use the mean value in the sampling window to fill the gap ("ave", default), or
#' 2. use any value assigned here to fill the gap (e.g., 9999, NA, etc.)
#' @param ... other arguments pass to `neuralnet`
#' @return A data frame that includes the original data, gap-filled data ("filled")
#' and a "mark" column that indicates the value in each row of the "filled" is either:
#' 1. original, 2. gap-filled, or 3. failed to converge
#' @examples
#' # read example data
#' df <- read.csv(file = system.file("extdata", "Soil_resp_example.csv", package = "FluxGapsR"),header = T)
#' df_filled <- Gapfill_ann(data = df,var1 = "Ts",var2 = "Ta",var3 = "Moist")
#' # visualize the gapfilled results
#' plot(df_filled$filled,col="red")
#' points(df_filled$Flux)
#' @export
Gapfill_ann <- function(data,
                        Flux = "Flux",
                        var1,
                        var2 = NULL,
                        var3 = NULL,
                        win = 5,
                        interval = 10,
                        threshold = 1,
                        hidden = 2,
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

  #  based on variable numbers
  if (is.null(var2)){ # if one variable
    # extract the data needed for gap-filling
    dft <- data[,c(Flux,var1)]
    names(dft) <- c("Flux","var1")
    # scale and normalize the input variables
    dft <- dft %>%
      mutate(var1=scale(var1))
    formula <- as.formula("Flux~var1")
  } else {
    if (is.null(var3)){ # if two variables
      # extract the data needed for gap-filling
      dft <- data[,c(Flux,var1,var2)]
      names(dft) <- c("Flux","var1","var2")
      # scale and normalize the input variables
      dft <- dft %>%
        mutate(var1=scale(var1),
               var2=scale(var2))
      formula <- as.formula("Flux~var1+var2")
    } else { # if three variables
      # extract the data needed for gap-filling
      dft <- data[,c(Flux,var1,var2,var3)]
      names(dft) <- c("Flux","var1","var2","var3")
      # scale and normalize the input variables
      dft <- dft %>%
        mutate(var1=scale(var1),
               var2=scale(var2),
               var3=scale(var3))
      formula <- as.formula("Flux~var1+var2+var3")
    }
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
    df_ann <- dft[wind_st:wind_ed,] %>%
      na.omit(.) # remove data in the gap

    # ANN model
    nn <- try(neuralnet::neuralnet(formula = formula,
                                   data = df_ann,
                                   # data = dft[sample(c(1:nrow(dft)),size = 2000),], ## sample a fraction of data for test
                                   threshold = threshold, # increase the threshold to improve the chance of converge
                                   stepmax = 1e+07, # increase the max step to improve the chance of converge
                                   hidden = hidden, #
                                   linear.output = T,...), # regression, not classification
                                   silent = TRUE)

    # predict the gaps
    if (class(fit)!="try-error"){ # if the fit converged
      gap[indx] <- predict(nn,newdata=dft[indx,])
      mark[indx] <- 1 # filled gap
      print(paste0("#",i," out of ",max(mk)," gaps: succeed!!")) # for checking progress
    } else {
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
    }
  } # end of the loop
  df_new <- data.frame(data,
                       filled = gap,
                       tem = dft[,"Flux"],
                       mark) %>%
    mutate(filled = ifelse(mark==0,tem,filled)) %>%
    select(-tem) # drop the temperory column

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
