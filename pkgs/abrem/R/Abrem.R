# R package 'abrem'
# Abernethy Reliability Methods
# Implementations of lifetime data analysis methods described in
# 'The New Weibull Handbook, Fifth edition' by Dr. Robert B. Abernethy.
# April 2014, Jurgen Symynck
# Copyright 2014, Jurgen Symynck
#
# For more info, visit http://www.openreliability.org/
#
# For the latest version of this file, check the Subversion repository at
# http://r-forge.r-project.org/projects/abernethy/
#
# Disclaimer:
#    The author is not affiliated with Dr. Abernethy or Wes Fulton - CEO of
#    Fulton Findings(TM) and author of the software package SuperSMITH
#-------------------------------------------------------------------------------
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# +-----------------------------------+
# |  execute this software with R:    |
# |  http://www.r-project.org/        |
# +-----------------------------------+

Abrem <- function(x,...){
    #arg <- list(...)
    arg <- splitargs(...)

#    arg2 <- arg[!(names(arg) %in% names(options.abrem()))]
        # extract the arguments that are NOT abrem options.
    opa <- modifyList(options.abrem(), arg$opa)
    ret <- list()
    class(ret) <- "abrem"
    timeorder <- c()
#    lowest <- 1
    if(!missing(x)){
        ret$data <- NULL
        if(is.vector(x) || is.numeric(x)){
            # TODO: the effects of adding is.numeric(x) are not tested thoroughly yet...
            if(opa$verbosity >= 2)message(
                'Abrem: Argument \"x\" is a (numeric) vector of (life-)time observations...')
            if(any(is.na(x))) timeorder <- 1:length(x)
            else timeorder <- order(x)
                # the above is tp prevens ordering attempts when NA values are
                # present in the lifetime observation vector.
                # having NA values implies that the data must be ordered.
            ret$data <- data.frame(time=x[timeorder],event=1)
        }
        if(is.data.frame(x)){
            if(!is.null(x$time) && !is.null(x$event)){
                if(opa$verbosity >= 2)message(
                    'Abrem: Argument \"x\" is a dataframe with $time and $event ',
                        'columns...')
                if(any(is.na(x$time))) timeorder <- 1:length(x$time)
                else timeorder <- order(x$time)
                ret$data  <- as.data.frame(x[timeorder,])
#                ret$data$event <- 1
#                    # temporarily set event vector to 1
            }else{
                stop(': Argument \"x\" is missing $time and/or ",
                    "$event columns...')
            }
        }
    }else{
        ti <- c(arg$rem$time,arg$rem$fail)
        if(xor(!is.null(arg$rem$time), !is.null(arg$rem$fail))){
            if(is.vector(ti)){
                if(opa$verbosity >= 2)message(match.call()[[1]],
                    'Abrem: Argument \"time\" or \"fail\" is vector of complete (life-)time observations...')
                if(any(is.na(ti))) timeorder <- 1:length(arg$rem$time)
                else timeorder <- order(ti)
                ret$data  <- data.frame(time=ti[timeorder],event=1)
            }else{stop('Argument \"time\" or fail\" must be vector.')}
        }
        if(!is.null(arg$rem$susp)){
            if(is.vector(arg$rem$susp)){
                if(opa$verbosity >= 2)message(
                    'Abrem: Argument \"susp\" is vector of right-censored (suspended) (life-)time observations...')
                timeorder <- order(c(ti,arg$rem$susp))
                ret$data  <- data.frame(time=c(ti,arg$rem$susp)[timeorder],
                    event=c(rep(1,length(ti)),rep(0,length(arg$rem$susp)))[timeorder])
            }else{stop('Argument \"susp\" must be a vector.')}
        }
        #else{stop("No (life-)time observations were provided.")}
    }
    
    ### setting the event vector correctly ###
    if(!is.null(arg$rem$event) && !is.null(ret$data)){
        if(is.vector(arg$rem$event)){
            if(opa$verbosity >= 2)message(
                'Abrem: Argument \"event\" is event vector...')
            ret$data$event <- arg$rem$event[timeorder]
        }
    }
    if("median" %in% opa$pp){
            if(opa$verbosity >= 2)message(
            'Abrem: Adding exact median ranks to (life-)time observations ...')
        ret$data <- cbind(ret$data,rank.median=NA)
        ret$data[ret$data$event==1,'rank.median'] <-
            .Call("medianRank1",ret$data$event, PACKAGE= "pivotals")
#        lowest <- min(c(lowest, na.omit(ret$data$rank.median)))
    }
    if("benard" %in% opa$pp){
            if(opa$verbosity >= 2)message(
            'Abrem: Adding Benards ranks to (life-)time observations ...')
        ret$data <- cbind(ret$data,rank.benard=NA)
        ret$data[ret$data$event==1,'rank.benard'] <-
            .Call("medianRank",ret$data$event, PACKAGE= "pivotals")
#        lowest <- min(c(lowest, na.omit(ret$data$rank.benard)))
    }



    ret$n    <- length(ret$data$time)
        # TODO: this assumes that any NA time (in any present
        # in the time column is there for a good reason:
        # accompanied with a censoring indicator (0 or FALSE)
        # TODO: check if the above code  is still valid!
    ret$fail <- sum(ret$data$event)
    ret$cens <- ret$n-ret$fail


    ret$options <- opa
        # always store a full copy of the options.abrem structure here
    ret
    # TODO: check what to do with the automatically added row names that are sometimes out of order
}