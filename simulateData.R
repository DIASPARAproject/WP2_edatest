
load("rn_rna.Rdata") #rn_rna, rn_rnag

load("calibration_dataset_delta_gamma_2020.Rdata")
library(mgcv)
library(dplyr)
library(gratia)

library(pbapply)

years=1990:2018

network <- rn_rna %>%
  filter(startsWith(as.character(emu),"FR") & startsWith(as.character(idsegment),"FR") &
           !as.character(emu) %in% c("FR_Rhin", "FR_Meus") &
           !is.na(riverwidthm) & !is.na(cs_height_10_n) & !is.na(cs_height_10_p) &
           !is.na(distanceseakm.) & !is.na(altitudem.) & !is.na(lengthriverm) &
           !is.na(lddws) & !is.na(hydraulicdensityperm2.) & !is.na(dist_from_gibraltar_km)) %>%
  select(-all_of(c("gamma","delta","density")))
rm("rn_rna")

library(ggplot2)
library(dplyr)
ddd |>
  group_by(emu, year) |>
  summarise(n = n()) |>
  mutate(cat = cut(n, breaks = c(0, 10, 20, 50, 100, 200, 500, 1000))) |>
  ggplot(aes(x = year, y = emu, fill = cat)) +
  geom_tile() + 
  theme_bw() +
  xlab("") + ylab("EMU") +
  scale_fill_viridis_d("number of operations") +
  xlim(1990, 2018)

ggsave("nb_elecrofishing.png", dpi = 300, width = 16 / 2.54, height = 16/2.54)

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
                          ctrl=list(nthreads=8),
                          method="REML")


model_gamma2 <- mgcv::gam(densCS ~ s(year, by=emu, k = 6)  + 
                            altitudem. * emu + 
                            riverwidthm + emu + 
                            te(cs_height_10_n., 
                               distanceseakm., by = emu, k = c(5, 5)),	
                          data=ddg,family = Gamma(link="log"),
                          ctrl=list(nthreads=8),samfrac=0.25,
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
                               ctrl=list(nthreads=1),
                               method="REML", start=coef(model_delta2)[!startsWith(names(coef(model_delta2)),"ef_fishing")])
  
  
  model_gamma_bis <- mgcv::gam(densCS ~ s(year, by=as.factor(emu), k = 6)  + 
                                 altitudem. * emu + 
                                 riverwidthm + emu + 
                                 te(cs_height_10_n., 
                                    distanceseakm., by = as.factor(emu), k = c(5, 5)),	
                               data=dataset %>% filter(densCS>0),family = Gamma(link="log"),
                               ctrl=list(nthreads=1),
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
library(parallel)
comb=expand.grid(nb=c(20,40,60,80,100),
                 f=c("create_sampling_strahler", "create_sampling_random"),
                 iter=1:100)

debut=Sys.time()

results = pblapply(seq_len(nrow(comb)),function(r) {
  nb <- comb$nb[r]
  fs <- comb$f[r]
  i <- comb$iter[r]
  print(paste("comb", nb, fs))
  
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
  write.table(res,paste("eda_res",
                        nb,
                        fs,
                        i,
                        ".csv",
                        sep ="_"),
              sep = ";", 
              row.names = FALSE,
              col.names = TRUE)
  
  return(res)
  
})
end=Sys.time()


save.image("eda_simulate.rdata")




##### graph to illsutrate


load("eda_simulate.rdata")


library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
allsegments <- read.table("segment_centroids.csv", header = TRUE, sep = ";")


get_maps <- function(stations){
  deltadata <- runif(nrow(stations)) <= predict(model_delta2, newdata=stations,
                                                type = "response") 
  
  sim_dist <- gratia:::get_family_rd(model_gamma2)
  scale <- model_gamma2[["scale"]]
  
  gammadata <- sim_dist(mu = predict(model_gamma2, newdata=stations, type = "response"),
                        scale = scale)
  
  simulated_data <- cbind.data.frame(stations,
                                     data.frame(densCS = gammadata * deltadata))
  
  
  
  segments <- simulated_data$idsegment
  seg <- allsegments |>
    dplyr::filter(idsegment %in% segments)
  
  france <- ne_countries(country = "France", scale = 10)
  ggplot(france) + geom_sf() + geom_point(data = simulated_data |>
                                            left_join(seg, by = "idsegment"),
                                          aes(x = x,
                                              y = y,
                                              col = densCS)) +
    scale_colour_viridis_c("simulated density") +
    coord_sf(xlim = c(-8, 10),
                    ylim = c(40, 52)) +
    xlab("") +
    ylab("") +
    theme_bw()
  
  
}
stations20 <- create_sampling_random(20)
y = 1998
stations60 <- create_sampling_strahler(60)


library(patchwork)
get_maps(stations20)/ get_maps(stations60)
ggsave("simulated_maps.png", width = 16/2.54, height = 16/2.54, dpi = 300)


