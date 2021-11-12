#' copy number intra-tumor heterogeneity (CNH) from a single copy number measurement
#'
#' MATLAB CNH.m function translation into R to calculate copy number intra-tumor heterogeneity (CNH) from a single copy number measurement.
#' In this method, relative segmented copy numbers are transformed to absolute copy numbers and the distance
#' to integer values is measured for each segment. This procedure is
#' performed for either a range of ploidies/purities or for fixed
#' ploidy/purity, in case either or both quantities are known from previous
#' measurements.
#' Copy number heterogeneity (CNH) is defined as the minimum average distance of
#' segments to the closest integer value, weighted by segment length.
#' In addition to CNH, this function returns the ploidy and purity
#' corresponding to the inferred CNH.

#
#' @param seg_val numeric vector with values of relative copy numbers per segment.
#' @param seg_len numeric vector with segment lengths.
#' @param ploidy tumour ploidy
#' @param purity sample purity, default = NULL for grid search

#' @return
#' 1st output argument: \code{CNH_out}      --> inferred CNH
#'
#' 2nd output argument: \code{ploidy_out}  -->  inferred ploidy for empty input ploidy, otherwise same as input ploidy.
#'
#' 3th output argument: \code{purity_out}   -->  inferred purity for empty input purity, otherwise same as input purity.
#'


#' @usage  CNH(seg_val,seg_len,ploidy,purity)
#' @examples
#' file = paste0(path.package("CNH"), "/SJBALL247_D.bed")
#' SJBALL247_D <- read.delim(file, header=FALSE)
#' seg_len <- SJBALL247_D[,3]-SJBALL247_D[,2]
#' seg_val<- SJBALL247_D[,4]

#' CNH(seg_val,seg_len)

#' @export


CNH <-function(seg_val ,seg_len, ploidy = seq(1.5,5,0.01), purity = seq(0.2,1,0.01)){
  # check if input seg_val and seg_len are vectors of equal size
  if(!( is.vector(seg_val) && is.vector(seg_len) && length(seg_val) == length(seg_len) && NCOL(seg_val) == 1)){
    stop('Segment values (1st input argument) and segment lengths (second input argument) appear not to be column vectors of equal length')
  }
  # specify default range of ploidy for grid search, if input ploidy is empty

    if (!( is.numeric(ploidy) & NCOL(ploidy) == 1 & all(ploidy > 0))){
      stop('Ploidy is not a positive scalar or empty')
    }
  # specify default range of purity purity for grid search, if input purity is empty
    # check if purity is a scalar between 0 and 1 argument
    if (!( is.numeric(purity) & NCOL(purity) == 1 & all(purity <= 1)) ){
      stop('Purity is not a positive scalar between 0 and 1 or empty')
    }


  Npurity<-length(purity)
  Nploidy<-length(ploidy)

  # make grid with all combinations of ploidy and purity, for the transformation of
  # measured relative copy number profile (seg_val) to absolute values (q) using
  # q = seg_val*a1+a2.
  purity_all<-expand.grid(purity,ploidy)[,1]
  ploidy_all<-expand.grid(purity,ploidy)[,2]
  a1<-a2<-NULL
  for (i in 1:length(ploidy)){

    a1[(1+(i-1)*Npurity):(Npurity*i)] <-
      (purity*ploidy[i]+2*(1-purity))/purity
    a2[(1+(i-1)*Npurity):(Npurity*i)] <-
      -2*(1-purity)/purity

  }

  # iniatilize output: CNH_out, ploidy_out and purity_out
  CNH_out = 1;
  purity_out = 0;
  ploidy_out = 0;

  # grid search over all ploidies and purities to infer CNH
  for (i in 1:(Nploidy*Npurity)){
    # transform relative copy numbers to absolute numbers
    q = a1[i]*seg_val+a2[i]

    # measure distance to closest integer value of each segment
    q_dist_down = q %% 1
    q_dist_up = 1-(q %% 1)
    q_dist_min = pmin(q_dist_down, q_dist_up)

    # calculate the mean distance of segments to integer values,
    # weighted by segment length
    CNHnew <- sum(q_dist_min * seg_len) / (sum(seg_len))

    # if the new weighted distance to integers is smaller than any
    # previously calculated, replace CNH_out, ploidy_out and purity_out with the new values.
    if (CNHnew < CNH_out){
      CNH_out <- CNHnew
      purity_out = purity_all[i]
      ploidy_out = ploidy_all[i]

    }
  }
  return(list(CNH=round(CNH_out,4),ploidy=ploidy_out, purity=purity_out))
}
# Usage:
# SJBALL247_D <- read.delim("~/Documents/Projects/20211111CIN_B_ALL/SJBALL247_D.bed", header=FALSE)
# seg_len <- SJBALL247_D[,3]-SJBALL247_D[,2]
# seg_val<- SJBALL247_D[,4]
#
# CNH(seg_val,seg_len)
# $CNH
# [1] 0.0376
#
# $ploidy
# [1] 1.5
#
# $purity
# [1] 0.99

# SJBALL264_D <- read.delim("~/Documents/Projects/20211111CIN_B_ALL/SJBALL264_D.bed", header=FALSE)
# seg_len <- SJBALL264_D[,3]-SJBALL264_D[,2]
# seg_val<- SJBALL264_D[,4]
#
# CNH(seg_val,seg_len)

# $CNH
# [1] 0.0554
#
# $ploidy
# [1] 2.21
#
# $purity
# [1] 0.98
