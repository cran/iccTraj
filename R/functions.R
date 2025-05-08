

#' Preprocessing of data
#'
#' @keywords internal
#' @param data A data frame with the steps of trajectories. It is assumed the time between steps is uniform. It must contain contain at least five columns: subject identifier, trip identifier, latitude, longitude, and time of the reading.
#' @param ID Character string indicating the name of the subjects column in the dataset.
#' @param trip Character string indicating the trip column in the dataset.
#' @param LON Numeric. Longitude readings.
#' @param LAT Numeric. Latitude readings.
#' @param time Numeric. Time of the readings.
#' @return A data frame

get_data<-function(data,ID,trip,LON,LAT,time){

  data <- as.data.frame(data)
  vid <- data[,ID]
  vtrip <- data[,trip]
  vLON <- data[,LON]
  vLAT <- data[,LAT]
  vtime <- data[,time]
  idtrip = paste(vid,vtrip,sep="_")
  x<-data.frame(ID=vid,trip=vtrip,LON=vLON,LAT=vLAT,time=vtime,idtrip=idtrip)
  return(x)
}

#' Creates a list containing the trajectories information.
#' @keywords internal
#' @param projection Projection string of class CRS-class.
#' @param origin Optional. Origin of the date-time. Only needed in the internal process to create an object of type POSIXct.
#' @return Every component of the list has the following subcomponents:
#' #' \itemize{
#'   \item An object of class Track.
#'   \item An object of class SpatialPoints.
#'   \item Subject identifier
#'   \item Trip identifier.
#' }

tr_gen<-function(data,projection=CRS("+proj=longlat"), origin="1970-01-01 UTC"){

  id_trip<-unique(data$idtrip)

  out<-1:length(id_trip) %>% map(function(i) {
    tr_dat <- data %>% filter(idtrip == id_trip[i])
    pts = SpatialPoints(tr_dat[c("LON","LAT")], proj4string=projection)
    time<-as.POSIXct(tr_dat$time,origin=origin)
    tr<-STIDF(pts, time, data.frame(idtrip=tr_dat$idtrip))
    trs = Track(tr)
    id<-unique(tr_dat$ID)
    trip<-unique(tr_dat$trip)
    return(list(trs,pts,id,trip))
  }

  )
  out
}


#' Computes extended Hausdorff distance between two trajectories.
#' @param pp1 Set of spatial points for the first trajectory. It can be a matrix of 2D points, first column x/longitude, second column y/latitude, or a SpatialPoints or SpatialPointsDataFrame object.
#' @param pp2 Set of spatial points for the second trajectory. It can be a matrix of 2D points, first column x/longitude, second column y/latitude, or a SpatialPoints or SpatialPointsDataFrame object.
#' @param q Quantile for the extended Hausdorff distance. Default value q=1 uses the maximum that leads to classical Hausdorff distance.
#' @return A numerical value with the distance.
#' @export
#' @references{
#' Magdy, N., Sakr, M., Abdelkader, T., Elbahnasy, K. (2015). Review on trajectory similarity measures. 10.1109/IntelCIS.2015.7397286.
#'
#' Min, D., Zhilin, L., Xiaoyong, C. (2007) Extended Hausdorff distance for spatial objects in GIS. International Journal of Geographical Information Science, 21:4, 459–475
#'
#' }
#' @examples
#' # Take two trajectories
#' library(dplyr)
#' library(sp)
#' sample_data<-gull_data %>% filter(ID %in% c(5107912,5107913), trip %in% c("V02","V01"))
#' tr1<-gull_data %>% filter((ID == 5107912) & (trip=="V02"))
#' tr2<-gull_data %>% filter((ID == 5107913) & (trip=="V01"))
#' pts1 = SpatialPoints(tr1[c("LONG","LAT")], proj4string=CRS("+proj=longlat"))
#' pts2 = SpatialPoints(tr2[c("LONG","LAT")], proj4string=CRS("+proj=longlat"))
#' # Hausdorff distance
#' HD(pts1,pts2,q=1)
#' # Median Hausdorff distance
#' HD(pts1,pts2,q=0.5)

HD<-function(pp1,pp2,q=1){
  ds<-spDists(pp1, pp2, longlat = TRUE)
  m1<-quantile(apply(ds,1,min),q)
  m2<-quantile(apply(ds,2,min),q)
  max(m1,m2)
}

#' Creates a data frame with the Hausdorff distance between two trajectories
#' @keywords internal
#' @param data A list created with \code{"tr_gen"} function
#' @param i Position in the data of the first trajectory
#' @param j Position in the data of the second trajectory
#' @param q Quantile for the extended Hausdorff distance. Default value q=1 uses the maximum that leads to classical Hausdorff distance.
#' @export
#' @return A data frame with the subjects and trips identifiers and the extended Hausdorff distance
HD_data<-function(data,i,j,q){
  d = HD(data[[i]][[2]],data[[j]][[2]],q=q)
  ID1<-data[[i]][[3]]
  ID2<-data[[j]][[3]]
  trip1<-data[[i]][[4]]
  trip2<-data[[j]][[4]]
  data.frame(ID1,trip1,ID2,trip2,d)
}


#' Creates a data frame with the discrete Fréchet distance between two trajectories
#' @keywords internal
#' @param data A list created with \code{tr_gen} function
#' @param i Position in the data of the first trajectory
#' @param j Position in the data of the second trajectory
#' @export
#' @return A data frame with the subjects and trips identifiers and the discrete Fréchet distance

FD_data<-function(data,i,j){
  d =frechetDist(data[[i]][[1]],data[[j]][[1]])
  ID1<-data[[i]][[3]]
  ID2<-data[[j]][[3]]
  trip1<-data[[i]][[4]]
  trip2<-data[[j]][[4]]
  data.frame(ID1,trip1,ID2,trip2,d)
}



#' Generates a data frame with the pairwise distances between trajectories
#' @keywords internal
#' @param data_list A list created with \code{"tr_gen"} function
#' @param parallel TRUE/FALSE value. Use parallel computation? Default value is TRUE.
#' @param distance Metric used to compute the distances between trajectories. Options are **H** for median Hausforff distance, and **F** for discrete Fréchet distance.
#' @param q Quantile for the extended Hausdorff distance. Default value q=0.5 leads to median Hausdorff distance.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{furrr_options}.
#' @return A data frame with the subjects and trips identifiers and their distances



genHD<-function(data_list,parallel=TRUE,distance=c("H","F"),q=0.5,future_seed=123){

  distance<-match.arg(distance)

  k<-length(data_list)

  progressr::handlers(global = FALSE)

  if (parallel == TRUE) {
    ncores <- parallelly::availableCores(omit = 1)

    oplan <- future::plan("multisession", workers = ncores)
    on.exit(future::plan(oplan))


    progressr::with_progress({
      p <- progressr::progressor(along = 1:(k-1))
      D <- furrr::future_map_dfr(1:(k-1), function(i) {
        p()
        Sys.sleep(.2)
        furrr::future_map_dfr((i+1):k, function(j) {
          if (distance == "H") {
            iccTraj::HD_data(data_list, i, j, q = q)
          } else if (distance == "F") {
            iccTraj::FD_data(data_list, i, j)
          }
        })
      }, .options = furrr::furrr_options(seed = future_seed))
    })

  }


  if (parallel==FALSE){

  progressr::with_progress({
    p <- progressr::progressor(along = 1:(k-1))


    D<- 1:(k-1) %>% map_df(function(i){
      p()
      Sys.sleep(.2)
      j=(i+1):k %>% map_df(function(j){
        if (distance=="H") {
          HD_data(data_list,i,j,q=0.5)
        } else if (distance=="F") {
          FD_data(data_list,i,j)
        }
      }
      )
    }
    )
  })
  }


  return(D)
}



#' Computes the number of trips by subject
#' @keywords internal
#' @param data A data frame with the steps of trajectories. It is assumed the time between steps is uniform. It must contain contain at least five columns: subject identifier, trip identifier, latitude, longitude, and time of the reading.
#' @return A data frame with the number of trips by subject.

ntrips<-function(data){
  id<-unique(data$ID)
  ntrips<-1:length(id) %>% map_dfr(function(i){
    data %>% filter(ID == id[i]) %>% group_by(ID,trip) %>%
      summarize(n=n(), .groups="drop") %>%
      ungroup() %>% group_by(ID) %>% tally()
  }
  )

  return(ntrips)

}

#' Creates a matrix with the pairwise distances
#' @keywords internal
#' @param dataD Data frame created with \code{"genHD"} function.
#' @param nt Data frame with the number of trips by subject
#' @return A list with two components:
#' \itemize{
#'   \item *D*. A matrix with the pairwise distances
#'   \item *data*. A data frame with the pairwise distances
#' }


dist_mat<-function(dataD,nt){
  s1<-dataD %>% dplyr::select(ID1,trip1,ID2,trip2,d)

  s2<-dataD %>% mutate(ID2b=ID1,trip2b=trip1,ID1=ID2,trip1=trip2,
                       ID2=ID2b, trip2=trip2b) %>%
    dplyr::select(ID1,trip1,ID2,trip2,d)

  s3<-dataD %>% group_by(ID1,trip1) %>% summarise(n=n(),.groups = 'drop') %>%
    mutate(ID2=ID1,trip2=trip1,d=0) %>% dplyr::select(-n)

  s4<-dataD %>% group_by(ID2,trip2) %>% summarise(n=n(),.groups = 'drop') %>%
    mutate(ID1=ID2,trip1=trip2,d=0) %>% dplyr::select(-n)


  S<-rbind(s1,s2,s3,s4[nrow(s4),]
  ) %>% arrange(ID1,trip1,ID2,trip2)

  mat<-matrix(S$d,nrow=sum(nt$n),byrow=TRUE)
  return(list(D=mat,data=S))

}




#' Computes the intraclass correlation coefficient (ICC) using a matrix of distances.
#' @param X Matrix with the pairwise distances.
#' @param nt Data frame with the number of trips by subject
#' @return Data frame with the estimates of the ICC (r), the subjects' mean sum-of-squares (MSA), the between-subjects variance (sb), the total variance (st), and the within-subjects variance (se).
#' @export
#' @details
#'The intraclass correlation coefficient is estimated using the distance matrix among trajectories.


ICC<-function(X,nt){

  ns<-nrow(nt)
  n<-sum(nt$n)

  n0<-(n-sum(nt$n^2)/n)/(ns-1)


  A<-(X^2)*(-0.5)
  J<-matrix(rep(1,n^2),nrow=n)
  I<-diag(n)

  P<-(I-J/n)
  G<-P%*%A%*%P

  xx<-1:length(nt$n) %>% map(function(i)
    matrix(rep(1,nt$n[i]^2),nrow=nt$n[i])/nt$n[i]
  )

  H<-Reduce(adiag,xx)-J/n

  zz<-H%*%G%*%H

  a1<-sum(diag(zz))
  a2<-sum(diag(G))
  se<-(a2-a1)/(n-ns)
  MSA=a1/(ns-1)
  sb<-(MSA-se)/n0
  st<-sb+se
  sb2=MSA/n0
  r=min(sb/st,1)
  out<-data.frame(MSA,sb,st,se,r)
  return(out)
}


#' Generates a bootstrap sample and estimates the ICC
#' @keywords internal
#' @param X Data frame with the pairwise distances
#' @param nt Data frame with the number of trips by subject
#' @param Bmat Matrix with subjects identifiers in the resamples
#' @param indB Column in Bmat that geenrates the resample
#' @export
#' @return Data frame with the estimates of the ICC (r), the subjects' mean sum-of-squares (MSA), the between-subjects variance (sb), the total variance (st), and the within-subjects variance (se).

boot_ICC<-function(X,nt,Bmat,indB){
  mos.id<-Bmat[,indB]
  # All pairs
  bb<-rbind(cbind(mos.id,mos.id),t(combn(mos.id,2)),t(combn(mos.id,2))[,2:1])
  # All pairs new IDs
  idb<-rbind(cbind(1:length(mos.id),1:length(mos.id)),
             t(combn(1:length(mos.id),2)),t(combn(1:length(mos.id),2))[,2:1])
  ord<-order(idb[,1],idb[,2])
  bb<-bb[ord,]
  idb<-idb[ord,]
  dd<-data.frame(ID=mos.id)
  # number of trips
  nt_boot<-dd %>% left_join(nt,by="ID")
  # All pairs in data frame
  bbdat<-data.frame(ID1=bb[,1],ID2=bb[,2],ID1b=idb[,1],ID2b=idb[,2])
  # Add distances
  temp<-bbdat %>% left_join(X,by=c("ID1","ID2")) %>% arrange(ID1b,trip1,ID2b,trip2)
  mat<-matrix(temp$d,nrow=sum(nt_boot$n),byrow=TRUE)
  ICC(mat,nt_boot)
}


#' Estimates the intraclass correlation coefficient (ICC) for trajectory data
#' @export
#' @param data A data frame with the locations and times of trajectories. It is assumed the time between locations is uniform. It must contain at least five columns: subject identifier, trip identifier, latitude, longitude, and time of the reading.
#' @param ID Character string indicating the name of the subjects column in the dataset.
#' @param trip Character string indicating the trip column in the dataset.
#' @param LON Numeric. Longitude readings.
#' @param LAT Numeric. Latitude readings.
#' @param time Numeric. Time of the readings.
#' @param projection Projection string of class CRS-class.
#' @param origin Optional. Origin of the date-time. Only needed in the internal process to create an object of type POSIXct.
#' @param parallel TRUE/FALSE value. Use parallel computation? Default value is TRUE.
#' @param individual TRUE/FALSE value. Compute individual within-subjects variances? Default value is TRUE.
#' @param distance Metric used to compute the distances between trajectories. Options are "H" for median Hausforff distance, and "F" for discrete Fréchet distance.
#' @param bootCI TRUE/FALSE value. If TRUE it will generate boostrap resamples. Default value is TRUE.
#' @param nBoot Numeric. Number of bootstrap resamples. Ignored if \code{bootCI} is FALSE. Default value is 100.
#' @param q Quantile for the extended Hausdorff distance. Default value q=0.5 leads to median Hausdorff distance.
#' @param future_seed Logical/Integer. The seed to be used for parallellization. Further details in \code{\link[furrr]{furrr_options}}.
#' @return An object of class \code{iccTraj}.The output is a list with the following components:
#' \itemize{
#'   \item \code{est}. Data frame with the following estimates: the ICC (r), the subjects' mean sum-of-squares (MSA), the between-subjects variance (sb), the total variance (st), and the within-subjects variance (se).
#'   \item \code{boot}. If bootCI argument is set to TRUE, data frame with the bootstrap estimates.
#'   \item \code{D}. Data frame with the pairwise distances among trajectories.
#'   \item \code{indW}. Data frame with the following columns: the subject's identifier (ID), the individual within-subjects variances (w), the individual ICC (r), and the number of trips (n).
#' }
#' @details
#' The intraclass correlation coefficient is estimated using the distance matrix among trajectories.
#'
#' Bootstrap resamples are obtained using balanced randomized cluster bootstrap approach (Davison and Hinkley, 1997; Field and Welsh, 2007)
#'
#' @references{
#'
#' Davison A.C., Hinkley D.V. (1997). Bootstrap Methods and Their Application. Cambridge: Cambridge University Press.
#'
#' Field, C.A., Welsh, A.H. (2007). Bootstrapping Clustered Data. Journal of the Royal Statistical Society. Series B (Statistical Methodology). 69(3), 369-390.
#'
#' }
#' @examples
#'\donttest{
#' # Using median Hausdorff distance.
#'  Hd<-iccTraj(gull_data,"ID","trip","LONG","LAT","triptime")
#'  Hd$est
#' # Using discrete Fréchet distance.
#' Fd<-iccTraj(gull_data,"ID","trip","LONG","LAT","triptime", distance="F")
#' Fd$est
#'}

iccTraj<-function(data,ID,trip,LON,LAT,time,projection=CRS("+proj=longlat"),
                  origin="1970-01-01 UTC",parallel=TRUE, individual=TRUE, distance=c("H","F"),bootCI=TRUE,
                  nBoot=100,q=0.5,future_seed=123){

  t1<-Sys.time()


  data_tr<-get_data(data,ID,trip,LON,LAT,time)


  # Creates a list with the trajectories as spatial points objects
  # for subject and trip
  message("Processing data...")
  data_list<-tr_gen(data_tr,projection,origin)

  #Number of trips by subject
  nt<-data_tr %>% group_by(ID) %>% summarise(n=length(unique(trip)))

  message("Generating distances...")
  # Generates dataframe with pairwise distances
  D1<-genHD(data_list,parallel=parallel,distance,q=q,future_seed = future_seed)

  # Generates distance matrix
  Dmat<-dist_mat(D1,nt)
  D2<-Dmat$D
  X<-Dmat$data
  est<-ICC(D2,nt)

  #Bootstrap
  id<-unique(data_tr$ID)

  Bmat<-matrix(sample(rep(id,nBoot)),nrow=length(id),byrow=FALSE)

  if (bootCI==TRUE){
    message("Bootstrapping...")


    ncores <- parallelly::availableCores(omit = 1)

    oplan <- future::plan("multisession", workers = ncores)
    on.exit(future::plan(oplan))

    #progressr::handlers(global = FALSE)

    progressr::with_progress({
      p <- progressr::progressor(along = 1:nBoot)


      ICC_boot <- furrr::future_map_dfr(1:nBoot, ~{
        p()
        Sys.sleep(.2)
        boot_ICC(X,nt,Bmat,.x)
      }, .options = furrr::furrr_options(seed = future_seed),p=p)
    })


    message(paste(nBoot,"bootstrap samples generated",sep=" "))

  } else if (bootCI == FALSE) {
    ICC_boot <- NULL
  }


  if (individual == TRUE){
    message("Computing individual within-subjects variances...")

    wind_data<-1:length(id) %>% map_df(function(i){

      wind<-within_ind_var(D1,id[i])
      data.frame(ID=id[i],w=wind$w,r=est$sb/(est$sb+wind$w),n=wind$n)
    })

  } else if (individual == FALSE) {
    wind_data <- NULL
  }


  t2<-Sys.time()


  dt<-round(as.numeric(difftime(t2,t1)),2)

  message(paste("Process took",dt,attr(difftime(t2,t1), "units")
                ,sep=" "))
  out<-list(est=est,boot=ICC_boot,D=D1,indW=wind_data)
  class(out)<-c("iccTraj","list")
  return(out)
}

#' Prints the ICC
#' @keywords internal
#' @param x An object of class \code{"iccTraj"}
#' @return The ICC estimate.
#' @export
print.iccTraj<-function(x,...){
  print(x$est$r)
}


#' Computes the confidence interval for the ICC
#' @export
#' @param x An object of class \code{iccTraj}
#' @param conf Numeric. Level of confidence. Default is set to 0.95.
#' @param method String. Method used to estimate the confidence interval. Accepted values are: "BCa" for bias-corrected and accelerated bootstrap,  "EB" for empirical bootstrap, "Perc" for percentile bootstrap, "AN" for asymptotic Normal, and "ZT" for asymptotic Normal using the Z-transformation.
#' @return A vector with the two boundaries of the confidence interval.
#' @details
#' Let \eqn{\hat{\theta}} denote the ICC sample estimate and \eqn{\theta_i^{B}} denote the ICC bootstrap estimates with \eqn{i=1,\ldots,B}. Let \eqn{\delta_{\alpha/2}^{B}} and \eqn{\delta_{1-\alpha/2}^{B}} be the \eqn{\frac{\alpha}{2}} and \eqn{1-\frac{\alpha}{2}} percentiles of \eqn{\delta_{i}^{B}=\theta_i^{B}-\hat{\theta}}.
#'
#' The percentile bootstrap confidence interval is computed as \eqn{\hat{\theta}+\delta_{\alpha/2}^{B},\hat{\theta}+\delta_{1-\alpha/2}^{B}}.
#'
#' The empirical bootstrap confidence interval is estimated as \eqn{\hat{\theta}-\delta_{1-\alpha/2}^{B},\hat{\theta}-\delta_{\alpha/2}^{B}}
#'
#' Asymptotic Normal (AN) interval is obtained as \eqn{\hat{\theta} \pm Z_{1-\alpha/2}*SE_B} where \eqn{SE_B} denotes the standard deviation of \eqn{\theta_i^{B}}, and \eqn{Z_{1-\alpha/2}} stands for the \eqn{1-\alpha/2} quantile of the standard Normal distribution.
#'
#' In the ZT approach, the ICC is transformed using Fisher's Z-transformation. Then,  the AN approach is applied to the transformed ICC.
#'
#' @examples
#'\donttest{
#' # Using median Hausdorff distance
#' Hd<-iccTraj(gull_data,"ID","trip","LONG","LAT","triptime", parallel=FALSE, distance="H")
#' Hd$est
# Empirical bootstrap interval
#' interval(Hd)
#'}

interval<-function(x,conf=0.95,method=c("BCa","Perc","EB","AN","ZT")){

  method<-match.arg(method)
  q <- (1-conf)/2


  if (class(x)[1]=="iccTraj"){
    if (is.null(x$boot)){
      message("No bootstrap samples. Confidence interval is not computed")
      ci<-NULL
      }
    else if (!is.null(x$boot)){

      # Empirical bootstrap
      if (method=="Perc"){
        ci<-x$est$r+quantile(x$boot$r-x$est$r,probs=c(q,1-q))
        met<-"Percentile"
      }
      # Asymptotic Normal
      if (method=="AN"){
        se_r<-sd(x$boot$r)
        ci<-x$est$r+c(-1,1)*qnorm(1-q)*se_r
        met<-"Asymptotic Normal"
      }
      # Z transformation
      if (method=="ZT"){
        se_z<-sd(atanh(x$boot$r))
        ci<-tanh(atanh(x$est$r)+c(-1,1)*qnorm(1-q)*se_z)
        met<-"Z-transformation"
      }
      if (method=="EB"){
        ci<-2*x$est$r-quantile(x$boot$r,probs=c(q,1-q))
        ci<-ci[2:1]
        names(ci)<-names(ci)[2:1]
        met<-"Empirical bootstrap"
      }
      if (method=="BCa"){
        ci<-bca(x$boot$r, x$est$r, cl = conf)
        met<-"BCa"
      }
      message("Method:",met,sep="")
    }
  }
  else if (class(x)[1]!="iccTraj") {
    message("Object is not iccTraj class")
    ci<-NULL
   }


  return(ci)

}


within_ind_var<-function(data,id){

  dist_data<-data %>% filter(ID1==id, ID2==id)

  n<-(sqrt(nrow(dist_data)*8+1)+1)/2

  if (n>1){

  X<-matrix(0,nrow=n,ncol=n,byrow=T)
  X[lower.tri(X)]<-dist_data$d
  X<-t(X)
  X[lower.tri(X)]<-dist_data$d

  A<-(X^2)*(-0.5)
  J<-matrix(rep(1,n^2),nrow=n)
  I<-diag(n)

  P<-(I-J/n)
  G<-P%*%A%*%P
  w=sum(diag(G))/(n-1)
  } else if (n==1) w<-0
  data.frame(w,n=n)

}


bca<-function (theta, t0, cl = 0.95)
{

  theta<-theta[!is.na(theta)]
  cl_low <- (1 - cl)/2
  cl_hi <- 1 - cl_low
  nsims <- length(theta)
  z.inv <- length(theta[theta < t0])/nsims
  z <- qnorm(z.inv)
  U <- (nsims - 1) * (t0 - theta)
  A1 <- sum(U^3)
  A2 <- 6 * (sum(U^2))^{
    3/2
  }
  a <- A1/A2
  ll.inv <- pnorm(z + (z + qnorm(cl_low))/(1 - a * (z + qnorm(cl_low))))
  ll <- quantile(theta, ll.inv, names = FALSE)
  ul.inv <- pnorm(z + (z + qnorm(cl_hi))/(1 - a * (z + qnorm(cl_hi))))
  ul <- quantile(theta, ul.inv, names = FALSE)
  return(c(ll, ul))
}
