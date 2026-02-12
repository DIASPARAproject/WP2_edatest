

ddatawd=modelwd="~/OneDrive_1_23-02-2024/"
library(mgcv)
library(dplyr)
library(gratia)

years=1990:2018

#load("eda_simulate.rdata")
myfiles <- list.files() 
myfiles <- myfiles[startsWith(myfiles, "eda_res_")]
results = do.call(bind_rows,lapply(myfiles, function(f) read.table(f, sep = ";", header=TRUE)))



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



results |>
  mutate(strategy  = gsub("create_sampling_", "", f)) |>
  filter(emu %in% c("FR_Loir", "FR_Garo")) |>
  group_by(emu, nb, year, strategy) %>%
  summarize(Binf = quantile(B, 0.025),
            Bsup = quantile(B, 0.975),
            .groups = "keep") %>%
  ungroup() |>
  left_join(biomass_true_state) |>
  ggplot(aes(x = year,
             y = B,
             fill = strategy)) + 
  geom_ribbon(aes(ymin = Binf, ymax = Bsup), alpha = .21) +
  geom_point() +
  facet_grid(emu~nb, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  xlab("") + 
  ylab("Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("trends.png", width = 16/2.54, height = 16/2.54, dpi = 300)

save(correlations, file = "eda_simulate2.rdata")




results %>% left_join(biomass_true_state, by = c("emu","year"), suffix(".iter",".true"))  |>
  filter(emu == "FR_Loir") 
  