mu = 0.8
var = 0.01^2

alpha = ((1-mu)/var - 1/mu)*mu^2
beta = alpha*(1/mu - 1)

alpha = 11#1300
beta = 91#320

tibble(x = seq(0, 1, length.out=1000),
       prior = dbeta(x, alpha, beta)) %>%
  ggplot(aes(x = x, y = prior)) +
  geom_line()
