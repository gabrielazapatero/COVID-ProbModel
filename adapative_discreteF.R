library("lubridate")
library("rstan")
library("readxl")
library("tidyverse")
library("zoo")
library("ggplot2")

# stan setup --------------------------------------------------------------

STAN_MODEL = "adaptive_discreteSF.stan"
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
set.seed(3) # for reproductibility

# -------------------------------------------------------------------------

setwd("~/Desktop/4 BIOMÉDICA/SEM 1/PROY3")
data_file = "datMad.xlsx"


# preproc -----------------------------------------------------------------

covid = read_excel(data_file) %>%
  # eliminate columns 1, and 4 to 7
  select(-1, -(4:7)) %>%
  rename(date = "fecha", cases = "num_casos", deaths = "Muertes_diarias",
         recovered = "Altas")

plot_covid = function(covid) {
  p = ggplot(
    covid %>% select(-starts_with("Cumul")) %>% gather(type, measure, cases:deaths),
    aes(x = date, y = measure, col = type)
  ) + geom_point() +
    geom_line() +
    facet_wrap(. ~ type, scales = "free_y")
  print(p)
}

# inspection --------------------------------------------------------------
###############################################################################

plot_covid(covid)

# Model will be focused on the second wave of the pandemic from 15-7 to 21-11 to be more efficient

## FIXME
covid = covid %>% filter(date > '2020-07-15')
#plot_covid(covid)
################################################################################


# # En los datos se ven oscilaciones asociadas a los días de la semana. Para evitarlos,
# # vamos a suavizar la señal con una media por ventanas.

plot(deaths ~ date, data = covid, type = 'l')
covid$deaths = as.integer(
  rollapply(covid$deaths, width = 5, mean, fill = "extend")
)
lines(deaths ~ date, data = covid, col = 2)

# Ídem para recovered
plot(recovered ~ date, data = covid, type = 'l')
covid$recovered = as.integer(
  rollapply(covid$recovered, width = 7, mean, fill = "extend")
)
lines(recovered ~ date, data = covid, col = 2)

# ídem para cases
plot(cases ~ date, data = covid, type = 'l')
covid$cases = as.integer(
  rollapply(covid$cases, width = 7, mean, fill = "extend")
)
lines(cases ~ date, data = covid, col = 2)



covid = covid %>%
  mutate(Cumulative_cases = cumsum(cases),
         Cumulative_deaths = cumsum(deaths), 
         Cumulative_recovered = cumsum(recovered), 
         totalRec = (cumsum(cases)-cumsum(deaths)))
plot_covid(covid)


# -------------------------------------------------------------------------

N =  6747425; # Total population of Madrid
n_days = nrow(covid)
covid$nday = 1:nrow(covid)
d0 = 1


# data for Stan
data_sir = with(covid, {
  list(
    n_days = n_days,
    t0 = 0,
    ts = nday,
    N = N,
    #------------> Solo usamos cumulative deaths para inferir resultados.
    cum_deaths = Cumulative_deaths,
    infected = cases
  )
})

###############################################################################
niter = 2000
chains = 2
###############################################################################
model = stan_model(STAN_MODEL)

fitted = sampling(model,
                  data = data_sir,
                  iter = niter,
                  chains = chains,
                  #------------------> Limit control so samples are not rejected
                  control = list(adapt_delta = 0.99))



# check results ------------------------------------------------------------

# 1) Predicted deaths
smr_pred = cbind(
  as.data.frame(
    summary(fitted, pars = "pred_deaths", 
            probs = c(0.05, 0.5, 0.95))$summary
  ),
  covid$nday,
  covid$Cumulative_deaths
)
colnames(smr_pred) = make.names(colnames(smr_pred)) # to remove % in the col names

print(
  ggplot(smr_pred, mapping = aes(x = covid.nday)) +
    geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
    geom_line(mapping = aes(x = covid.nday, y = X50.)) +
    geom_point(mapping = aes(y = covid.Cumulative_deaths), col = 2) +
    labs(x = "Day", y = "Deaths")
)



# 2) inferred predictions

smr_pred2 = cbind(
  as.data.frame(
    summary(fitted, pars = c("real_infected"), probs = c(0.05, 0.5, 0.95))$summary
  ),
  #-------> Add + 18 to take into account time to death delay
  'nday' = covid$nday - 18
)

smr_pred2 = merge(
  smr_pred2, 
  covid %>% select(nday, cases), 
  by = 'nday'
)

colnames(smr_pred2) = make.names(colnames(smr_pred2)) # to remove % in the col names

print(
  ggplot(smr_pred2, mapping = aes(x = nday)) +
    #geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
    geom_line(mapping = aes(x = nday, y = X50.)) +
    geom_point(mapping = aes(y = cases), col = 2) +
    geom_line(mapping = aes(y = cases), col = 2) +
    labs(x = "Day", y = "(real) Infected") 
)


# 3) R0
smr_pred3 = cbind(
  as.data.frame(
    summary(fitted, pars = "R0", probs = c(0.05, 0.5, 0.95))$summary
  ), 
  covid$nday
)
colnames(smr_pred3) = make.names(colnames(smr_pred3)) # to remove % in the col names

print(
  ggplot(smr_pred3, mapping = aes(x = covid.nday)) +
    geom_ribbon(aes(ymin = X5., ymax = X95.), alpha = 0.6) +
    geom_line(mapping = aes(x = covid.nday, y = X50.)) + 
    labs(x = "Day", y = "R0") + 
    geom_hline(yintercept=1, col = 2)
)

