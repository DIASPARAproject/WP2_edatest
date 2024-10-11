
library(mgcv)
library(dplyr)
library(gratia)

years=1990:2018

load("eda_simulate.rdata")

nety0 <- network %>%
  mutate(year = min(years),
         ef_wetted_area = 600) 

true_state_ref <- nety0 %>%
  group_by(emu) %>%
  filter(idsegment == min(idsegment, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(reflogit = predict(model_delta2, newdata=., type="link"),
         reflog = predict(model_gamma2, newdata=., type="link"))


nety0$gammalog=predict(model_gamma2, newdata=nety0,type="link")
nety0$deltalogit=predict(model_delta2,newdata=nety0,type="link")
nety0$gamma=model_gamma2$family$linkinv(nety0$gammalog)
nety0$delta=model_delta2$family$linkinv(nety0$deltalogit)


biomass_true_state = do.call(bind_rows,lapply(years, function(y){
  if (y == min(years)){
    nety <- nety0
  } else{
    ref2 <- true_state_ref %>%
      mutate(year=y) %>%
      mutate(ylogit = predict(model_delta2, newdata=., type="link"),
             ylog = predict(model_gamma2, newdata=., type="link")) %>%
      mutate(offsetlogit=ylogit-reflogit,
             offsetlog=ylog-reflog) %>%
      select(emu,offsetlog,offsetlogit)
    nety <- nety0 %>%
      mutate(year=y) %>%
      left_join(ref2, by = "emu") %>%
      mutate(gamma = model_gamma2$family$linkinv(gammalog+offsetlog),
             delta = model_delta2$family$linkinv(deltalogit+offsetlogit))
  }
  nety %>%
    group_by(emu,year ) %>%
    summarize(B=sum(delta*gamma*lengthriverm * riverwidthm,na.rm = TRUE), .groups =  "keep") %>%
    ungroup()
}))



correlations <- do.call(bind_rows,lapply(results, function(res){
  res %>% left_join(biomass_true_state, by = c("emu","year"), suffix(".iter",".true")) %>%
    group_by(emu,iter,f,nb) %>%
    summarize(cor = cor(B.x,B.y), .groups = "keep") %>%
    ungroup()
}))
save.image("eda_simulate2.rdata")
