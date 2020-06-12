# Some functions for peak picking and smoothing

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

## 3. closest()
## Finds the above and below closest value
closest <- function(x,value){
  b <- x[x < value]
  if(length(b) == 0){
    below <- NA
  } else {
    below <- b[which.min(abs(b-value))]
  }
  
  a <- x[x > value]
  if(length(a) == 0){
    above <- NA
  } else {
    above <- a[which.min(abs(a-value))]
  }
  return(c(below, above))
}



peak_detect <- function(x, y, ...) {
  
  model <- loess(y ~ x, ...) # compute model
  #model <- mgcv::gam(y ~ s(x), family = Gamma(link = log))
  
  
  frq.x <- seq(min(x),max(x),by = 1)
  # make a x-axis with a higher frequency of points to make rollapply work
  #y.pred <- predict(model, newdata = frq.x)
  # predict y-values along high frequency x axis using the loess model
  newd <- data.frame(x = frq.x)
  y.pred <- as.numeric(predict(model, newdata = newd)) # backtransform log with exp()

  i.max <- localMaxima(y.pred)
  # find peaks
  if(y.pred[i.max] < mean(y)){
    t.max <- try(i.max[y.pred[i.max] > mean(y)])
  } else {
    t.max <- i.max
  }
  
  # decision tree, depending on the maximum y of the identified peaks
  # because rare ASVs with few observations tend to get many peaks that are actually not peaks
  # because ASVs with an extreme peak tend to get a small peak where there is none
  # all other ASVs will take peaks that are higher than 1i
  if(length(t.max) == 0L){
    #error handling code, maybe just skip this iteration using
    #output
    list(
      model = model,
      data = data.frame(
        x.pred = frq.x,
        y.pred = y.pred,
        stringsAsFactors = F
      ),
      pk.df = data.frame(no.peak = 0,
                         peak = NA,
                         peak.x = NA,
                         peak.y = NA,
                         closest = NA,
                         infp.x = NA,
                         infp.y = NA,
                         infp.slope = NA,
                         auc = NA,
                         stringsAsFactors = F)
    )
  } else {
    #####################
    # Find inflection points
    deriv3rd <- diff(diff(diff(y.pred)))
    d <- sign(deriv3rd)
    d <- cumsum(rle(d)$lengths)
    inf.pts <- d[seq.int(1L, length(d), 2L)]
    
    ####################
    # Take inflection points before and after each peak
    df.ls <- list()
    
    for(i in 1:length(t.max)){
      clos.inf <- closest(inf.pts, t.max[i])
      
      temp <- data.frame(no.peak = rep(length(t.max), times = length(clos.inf)),
                         peak = rep(i, times = length(clos.inf)),
                         peak.x = rep(frq.x[t.max[i]], times = length(clos.inf)),
                         peak.y = rep(y.pred[t.max[i]], times = length(clos.inf)),
                         closest = c("above", "below"),
                         infp.x = frq.x[clos.inf],
                         infp.y = y.pred[clos.inf],
                         infp.slope = diff(y.pred)[clos.inf],
                         auc = flux::auc(frq.x, y.pred, thres = min(y)),
                         stringsAsFactors = F)
      df.ls[[i]] <- temp
    }
    
    df <- do.call(rbind, df.ls)
    
    #output
    list(
      model = model,
      data = data.frame(
        x.pred = frq.x,
        y.pred = y.pred,
        stringsAsFactors = F
      ),
      pk.df = df
    )
    
    # Export values
    # model = model
    # x.pred = x-axis used for prediction
    # y.pred = predicted y values by model
    # x.pks = x value of detected peaks
    # y.pks = y value of detected peaks
    
  }
}


plotPks <- function(x, y, span, ID) {
  peaks <- peak_detect(x, y, span = span)
  #saveRDS(peaks, paste0("./Output/ASV_loess/",rids[j],".rds"))
  
  p <- ggplot() +
    theme_bw() +
    geom_point(mapping = aes(x = x, y = y), alpha = .8, colour = "Grey40", shape = 21, size = 2) +
    geom_line(data = peaks$data, mapping = aes(x = x.pred, y = y.pred)) +
    geom_point(data = peaks$pk.df, mapping = aes(peak.x, peak.y), colour = "red", size = 4) +
    geom_point(data = peaks$pk.df, mapping = aes(infp.x, infp.y, colour = closest), size = 2.5) +
    labs(x = "Catchment Area", y = "Mean read count (100 rarefactions)",
         title = paste("ID = ", ID,", span = ", span, sep="")) +
    guides(colour=guide_legend(title="Closest \nInflection Point"))
  
  ggsave(paste0("./Figures/Chapter1/ASV_smooth/LOESS/",rids[j],"_span",span,".png"), p, height = 10, width = 13, units = "cm")
}




peak_detect_gam <- function(x, y, total, ...) {
  #model <- loess(y ~ x, ...) # compute model
  #model <-
  #  try(mgcv::gam(y ~ s(x, bs = "ts"), weights = total,
   #               family = nb()), silent = T)
  
  model <-
    try(mgcv::gam(y ~ s(x, bs = "ts"), weights = total,
                  family = betar(link = "logit")), silent = T)
  # not so sure quasi is very curvy too
  # negative binomial is too sensitive
  # quasi is a good middle curvature
  # "ds" follows better the general trend than "ts"
  
  if (inherits(model, "try-error")) {
    frq.x <- seq(min(x), max(x), by = 1)
    y.pred <- rep(0, length(frq.x))
    
    list(
      model = NA,
      data = data.frame(
        x.pred = frq.x,
        y.pred = y.pred,
        stringsAsFactors = F
      ),
      pk.df = data.frame(
        no.peak = 0,
        peak = NA,
        peak.x = NA,
        peak.y = NA,
        closest = NA,
        infp.x = NA,
        infp.y = NA,
        infp.slope = NA,
        auc = NA,
        stringsAsFactors = F
      ))
  } else {
    frq.x <- seq(min(x), max(x), by = 1)
    # make a x-axis with a higher frequency of points to make rollapply work
    #y.pred <- predict(model, newdata = frq.x)
    # predict y-values along high frequency x axis using the loess model
    newd <- data.frame(x = frq.x)
    
    if(model$family[2] == "log" | model$family[2] == "logit"){
      y.pred <-
        exp(as.numeric(predict(model, newdata = newd))) # backtransform log with exp() for nb  
    } else {
      y.pred <-
        as.numeric(predict(model, newdata = newd)) # no backtransformation for quasi
    }
    
    i.max <- localMaxima(y.pred)
    # find peaks
    t.max <- i.max
    #if (y.pred[i.max] < quantile(y.pred)[2]) {
    #  t.max <- try(i.max[y.pred[i.max] > quantile(y.pred)[2]])
    #} else {
    #  t.max <- i.max
    #}
    # decision tree, depending on the maximum y of the identified peaks
    # because rare ASVs with few observations tend to get many peaks that are actually not peaks
    # because ASVs with an extreme peak tend to get a small peak where there is none
    # all other ASVs will take peaks that are higher than 1i
    if (length(t.max) == 0L) {
      #error handling code, maybe just skip this iteration using
      #output
      list(
        model = model,
        data = data.frame(
          x.pred = frq.x,
          y.pred = y.pred,
          stringsAsFactors = F
        ),
        pk.df = data.frame(
          no.peak = 0,
          peak = NA,
          peak.x = NA,
          peak.y = NA,
          closest = NA,
          infp.x = NA,
          infp.y = NA,
          infp.slope = NA,
          auc = NA,
          stringsAsFactors = F
        )
      )
    } else {
      #####################
      # Find inflection points
      deriv3rd <- diff(diff(diff(y.pred)))
      d <- sign(deriv3rd)
      d <- cumsum(rle(d)$lengths)
      inf.pts <- d[seq.int(1L, length(d), 2L)]
      
      ####################
      # Take inflection points before and after each peak
      df.ls <- list()
      
      for (i in 1:length(t.max)) {
        clos.inf <- closest(inf.pts, t.max[i])
        if (length(na.omit(clos.inf)) == 0L) {
          temp <-
            data.frame(
              no.peak = rep(length(t.max), times = length(clos.inf)),
              peak = rep(i, times = length(clos.inf)),
              peak.x = rep(frq.x[t.max[i]], times = length(clos.inf)),
              peak.y = rep(y.pred[t.max[i]], times = length(clos.inf)),
              closest = c("above", "below"),
              infp.x = NA,
              infp.y = NA,
              infp.slope = NA,
              #auc = flux::auc(frq.x, y.pred, thres = min(y)),
              stringsAsFactors = F
            )
        } else {
          temp <-
            data.frame(
              no.peak = rep(length(t.max), times = length(clos.inf)),
              peak = rep(i, times = length(clos.inf)),
              peak.x = rep(frq.x[t.max[i]], times = length(clos.inf)),
              peak.y = rep(y.pred[t.max[i]], times = length(clos.inf)),
              closest = c("above", "below"),
              infp.x = frq.x[clos.inf],
              infp.y = y.pred[clos.inf],
              infp.slope = diff(y.pred)[clos.inf],
              #auc = flux::auc(frq.x, y.pred, thres = min(y)),
              stringsAsFactors = F
            )
          
        }
        
        df.ls[[i]] <- temp
      }
      df <- do.call(rbind, df.ls)
      
      #output
      t<- list(
        model = model,
        data = data.frame(
          x.pred = frq.x,
          y.pred = y.pred,
          stringsAsFactors = F
        ),
        pk.df = df
      )
      
      # Export values
      # model = model
      # x.pred = x-axis used for prediction
      # y.pred = predicted y values by model
      # x.pks = x value of detected peaks
      # y.pks = y value of detected peaks
    }
  }
}


plotPks.gam <- function(x, y, total, ID) {
  peaks <- peak_detect_gam(x, y, total)
  #saveRDS(peaks, paste0("./Output/ASV_loess/",rids[j],".rds"))
  
  p <- ggplot() +
    theme_bw() +
    geom_line(data = peaks$data, mapping = aes(x = x.pred, y = y.pred)) +
    geom_point(data = peaks$pk.df, mapping = aes(peak.x, peak.y), colour = "red", size = 4) +
    geom_point(data = peaks$pk.df, mapping = aes(infp.x, infp.y, colour = closest), size = 2.5) +
    geom_point(mapping = aes(x = x, y = y), alpha = .8, colour = "Grey40", shape = 21, size = 2) +
    labs(x = "Catchment Area", y = "Mean read count (100 rarefactions)",
         title = paste("ID = ", ID, sep="")) +
    guides(colour=guide_legend(title="Closest \nInflection Point"))
  
  ggsave(paste0("./Figures/Chapter1/ASV_smooth/GAMs/Rel.abund/",rids[j],"_gam.png"), p, height = 10, width = 13, units = "cm")
}


#if(max(y.pred[i.max]) < 1) {
#  t.max <- try(i.max[y.pred[i.max] == max(y.pred)])
#} else if(max(y.pred[i.max]) > 100){
#  t.max <- try(i.max[y.pred[i.max] > 30])
#} else {
#  t.max <- try(i.max[y.pred[i.max] > 1])
#}

model <-
  try(mgcv::gam(y ~ s(x, bs = "ts"), weights = total,
                family = mgcv::nb(link = "log")), silent = T)
model2 <-
  try(mgcv::gam(y ~ s(x, bs = "ts"), weights = total,
                family = quasi()))


frq.x <- seq(min(x), max(x), by = 1)
newd <- data.frame(x = frq.x)
y.pred <-
  exp(as.numeric(predict(model, newdata = newd))) # backtransform log with exp()
y.pred2 <-
  as.numeric(predict(model2, newdata = newd)) # backtransform log with exp()
y.pred3 <-
  exp(as.numeric(predict(model3, newdata = newd))) # backtransform log with exp()
# is better than tw, bs= "ts" better than "ps"
# quasi

ggplot(data.frame(x,y), aes(x,y)) +
  theme_bw() +
  geom_point(size = 3, colour = "grey40") +
  geom_line(data = data.frame(x =frq.x, y =y.pred), aes(x,y), size = 1) +
  geom_line(data = data.frame(x =frq.x, y =y.pred2), aes(x,y), size = 1, colour = "red") +
  geom_line(data = data.frame(x =frq.x, y =y.pred3), aes(x,y), size = 1, colour = "royalblue")
