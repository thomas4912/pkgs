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

abrem.conf <- function(x,which="all",...){
    # x is a single Abrem or a list of Abrem objects
#    supported_blifeconf  <-  c("mcpivotals","bbb")
    if(missing(x)){
        stop('Argument \"x\" is missing.')
    }else{
        if(identical(class(x),"abrem")) x <- list(x)
        if(!all(sapply(x,function(x)identical(class(x),"abrem")))){
            stop('Argument \"x\" is not of class \"abrem\" or ",
            "a list of \"abrem\" objects.')
        }
        dr <- findMaxDataRange(x,0)
            # for reasonably large confidence bounds
        calculateConfsInAbrem <- function(abrem){
            if(!is.null(abrem$fit)){
                abrem$fit <- lapply(abrem$fit,calculateSingleConf,
                    opadata=abrem$options,datarange=dr,...)
                    # TODO: add error detection
            }
            abrem
        }
        abremlist <- lapply(x,calculateConfsInAbrem)
    }
    if(length(abremlist)==1) abremlist[[1]] else abremlist
}

#            oppdata <- modifyList(oppdata,x$options.abremplot)
#            if(!is.null(x$fit)){
#                for.each <- function(fit){
#                    if(!is.null(fit$options.abremplot)){
#                        oppfit <- modifyList(oppdata,fit$options.abremplot)
#                    }else{oppfit <- oppdata}
#                    if(length(arg$opp)>0){
#                        oppconf <- modifyList(oppfit, arg$opp)
#                    }else{oppconf <- oppfit}
#                    if(!is.null(fit$options$dist)){
#                    }else{
#                        stop(match.call()[[1]],": Distribution type was not provided.")
#                    }
#                    fit
#                }
#                # if(identical(tolower(which),"all"))
#                # lapply(x$fit,conf,type)
#                x$fit[[length(x$fit)]] <- for.each(x$fit[[length(x$fit)]])
#                    # TODO: add support for which="all"
#            }else{
#                if(opa$verbosity >= 1)
#                    stop(match.call()[[1]],": Object does not contain anything under \"$fit\".")
#            }
#        }else{stop("Argument \"x\" is not of class \"abrem\".")}
#    }else{stop("Argument \"x\" is missing.")}

