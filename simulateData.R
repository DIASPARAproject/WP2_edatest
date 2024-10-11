
ddatawd=modelwd="~/OneDrive_1_23-02-2024/"
load(file=file.path(ddatawd,"rn_rna.Rdata")) #rn_rna, rn_rnag

load(paste0(ddatawd,"calibration_dataset_delta_gamma_2020.Rdata"))
library(mgcv)
library(dplyr)
library(gratia)

years=1990:2018

network <- rn_rna %>%
  filter(startsWith(as.character(emu),"FR") & startsWith(as.character(idsegment),"FR") &
           !as.character(emu) %in% c("FR_Rhin", "FR_Meus") &
           !is.na(riverwidthm) & !is.na(cs_height_10_n) & !is.na(cs_height_10_p) &
           !is.na(distanceseakm.) & !is.na(altitudem.) & !is.na(lengthriverm) &
           !is.na(lddws) & !is.na(hydraulicdensityperm2.) & !is.na(dist_from_gibraltar_km)) %>%
  select(-all_of(c("gamma","delta","density")))
rm("rn_rna")
ddd <- ddd %>%
  filter(startsWith(as.character(emu), "FR") & year >=1990 &
           !as.character(emu) %in% c("FR_Rhin", "FR_Meus"))

ddg  <- ddg %>%
  filter(startsWith(as.character(emu), "FR") & year >=1990 &
           !as.character(emu) %in% c("FR_Rhin", "FR_Meus"))
model_delta2 <- mgcv::gam(densCS>0 ~ s(year, by=as.factor(emu), k = 10)  + 
                            emu +
                            s(ef_wetted_area, k=5) +
                            s(altitudem., k=5) +		
                            s(dist_from_gibraltar_km, k=15) +
                            s(hydraulicdensityperm2., k=6) +
                            s(lddws, k=6)+
                            ef_fishingmethod+
                            s(cs_height_10_p., bs='tp', k=6)+
                            s(distanceseakm.,bs='tp', k=6),
                          family = binomial,
                          data=ddd,
                          ctrl=list(nthreads=4),
                          method="REML")


model_gamma2 <- mgcv::gam(densCS ~ s(year, by=emu, k = 6)  + 
                            altitudem. * emu + 
                            riverwidthm + emu + 
                            te(cs_height_10_n., 
                               distanceseakm., by = emu, k = c(5, 5)),	
                          data=ddg,family = Gamma(link="log"),
                          ctrl=list(nthreads=3),samfrac=0.25,
                          method="REML")





create_sampling_random <- function(nb){
  observe_ef_wetted <- ddd %>%
    select(idsegment, ef_wetted_area) %>%
    left_join(network %>%
                mutate(wetted = lengthm * riverwidthm) %>%
                select(idsegment, wetted), by = "idsegment") %>%
    mutate(prop_ef_wetted=ef_wetted_area / wetted) %>%
    select(prop_ef_wetted) %>%
    pull()
  network %>%
    group_by(emu) %>%
    slice_sample(n = nb) %>%
    ungroup()%>%
    mutate(prop_ef_wetted = sample(observe_ef_wetted, nrow(.))) %>%
    mutate(ef_wetted_area = prop_ef_wetted * lengthm * riverwidthm) %>%
    mutate(ef_wetted_area=pmax(pmin(ef_wetted_area, 3000, na.rm=TRUE),
                               1, na.rm=TRUE)) %>%
    mutate(ef_fishingmethod = "com")
  
}



create_sampling_strahler <- function(nb){
  observe_ef_wetted <- ddd %>%
    select(idsegment, ef_wetted_area) %>%
    left_join(network %>%
                mutate(wetted = lengthm * riverwidthm) %>%
                select(idsegment, wetted), by = "idsegment") %>%
    mutate(prop_ef_wetted=ef_wetted_area / wetted) %>%
    select(prop_ef_wetted) %>%
    pull()
  strahler_weight <- network %>%
    group_by(emu, strahler) %>%
    summarize(habitats =sum(lengthm*riverwidthm, na.rm=TRUE), .groups =  "keep") %>%
    ungroup() %>%
    group_by(emu) %>%
    mutate(weights = habitats/sum(habitats)) %>%
    ungroup() %>%
    select(emu, strahler, weights)
  network %>%
    left_join(strahler_weight, by= c("emu","strahler")) %>%
    group_by(emu) %>%
    slice_sample(n = nb, weight_by = weights) %>%
    ungroup() %>%
    mutate(prop_ef_wetted = sample(observe_ef_wetted, nrow(.))) %>%
    mutate(ef_wetted_area = prop_ef_wetted * lengthm * riverwidthm) %>%
    mutate(ef_wetted_area=pmax(pmin(ef_wetted_area, 3000, na.rm=TRUE),
                               1, na.rm=TRUE)) %>%
    mutate(ef_fishingmethod = "com")
}


simulate_data <- function(nb=20, sampling_function=create_sampling_strahler,years = 1990:2018){
  dataset <- do.call(bind_rows,lapply(years, function(y){
    stations <- sampling_function(nb)
    stations$year <- y
    
    deltadata <- runif(nrow(stations)) <= predict(model_delta2, newdata=stations,
                                                  type = "response") 
    
    sim_dist <- gratia:::get_family_rd(model_gamma2)
    scale <- model_gamma2[["scale"]]
    
    gammadata <- sim_dist(mu = predict(model_gamma2, newdata=stations, type = "response"),
                          scale = scale)
    
    simulated_data <- cbind.data.frame(stations,
                                       data.frame(densCS = gammadata * deltadata))
    
  }))
  dataset <- dataset %>%
    mutate(emu = factor(emu, levels=levels(model_delta2$model$emu)))
  
  model_delta_bis <- mgcv::gam(densCS>0 ~ s(year, by=as.factor(emu), k = 10)  + 
                                 emu +
                                 s(ef_wetted_area, k=5) +
                                 s(altitudem., k=5) +		
                                 s(dist_from_gibraltar_km, k=15) +
                                 s(hydraulicdensityperm2., k=6) +
                                 s(lddws, k=6)+
                                 s(cs_height_10_p., bs='tp', k=6)+
                                 s(distanceseakm.,bs='tp', k=6),
                               family = binomial,
                               data=dataset,
                               ctrl=list(nthreads=4),
                               method="REML", start=coef(model_delta2)[!startsWith(names(coef(model_delta2)),"ef_fishing")])
  
  
  model_gamma_bis <- mgcv::gam(densCS ~ s(year, by=as.factor(emu), k = 6)  + 
                                 altitudem. * emu + 
                                 riverwidthm + emu + 
                                 te(cs_height_10_n., 
                                    distanceseakm., by = as.factor(emu), k = c(5, 5)),	
                               data=dataset %>% filter(densCS>0),family = Gamma(link="log"),
                               ctrl=list(nthreads=3),
                               method="REML", start=coef(model_gamma2))
  
  nety0 <- network %>%
    mutate(year = min(years),
           ef_wetted_area = 600) 
  
  ref <- nety0 %>%
    group_by(emu) %>%
    filter(idsegment == min(idsegment, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(reflogit = predict(model_delta_bis, newdata=., type="link"),
           reflog = predict(model_gamma_bis, newdata=., type="link")) 
  
  nety0$gammalog=predict(model_gamma_bis, newdata=nety0,type="link")
  nety0$deltalogit=predict(model_delta_bis,newdata=nety0,type="link")
  nety0$gamma=model_gamma_bis$family$linkinv(nety0$gammalog)
  nety0$delta=model_delta_bis$family$linkinv(nety0$deltalogit)
  
  
  
  biomass = do.call(bind_rows,lapply(years, function(y){
    if (y == min(years)){
      nety <- nety0
    } else{
      ref2 <- ref %>%
        mutate(year = y) %>%
        mutate(ylogit = predict(model_delta_bis, newdata=., type="link"),
               ylog = predict(model_gamma_bis, newdata=., type="link")) %>%
        mutate(offsetlogit=ylogit-reflogit,
               offsetlog=ylog-reflog) %>%
        select(emu,offsetlog,offsetlogit)
      nety <- nety0 %>%
        mutate(year=y) %>%
        left_join(ref2, by=c("emu")) %>%
        mutate(gamma = model_gamma_bis$family$linkinv(gammalog+offsetlog),
               delta = model_delta_bis$family$linkinv(deltalogit+offsetlogit))
    }
    nety %>%
      group_by(emu,year) %>%
      summarize(B=sum(delta*gamma*lengthriverm * riverwidthm,na.rm = TRUE), .groups =  "keep") %>%
      ungroup()
  }))
  rm("dataset")
  biomass
}


debut=Sys.time()
comb=expand.grid(nb=c(20,40,60,80,100),
                 f=c("create_sampling_strahler", "create_sampling_random"))
results = mapply(function(nb,fs) {
  print(paste("comb", nb, fs))
  do.call(bind_rows,lapply(1:100,
                               function(i) {
                                 print(paste("comb", nb, fs,i))
                                 tryCatch({
                                   res <- simulate_data(nb, get(as.character(fs)), 1990:2018)
                                   res <- res %>%
                                            mutate(nb = nb, 
                                                   iter = i,
                                                   f = fs)},
                                   error=function(e){
                                   },
                                   finally = {res})
                                 return(res)
                               }))
}, comb$nb, comb$f,SIMPLIFY = FALSE)
end=Sys.time()


save.image("eda_simulate.rdata")


