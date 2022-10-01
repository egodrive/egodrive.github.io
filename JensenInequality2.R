# https://en.wikipedia.org/wiki/Gompertz_function
# library
library(ggplot2)
library(ggExtra)

a1 = rnorm(n = SampSize, mean = 1.5, sd = .1)
a2 = rnorm(n = SampSize, mean = 1.5, sd = .4)
a3 = rnorm(n = SampSize, mean = 1.5, sd = .8)

b1 = rnorm(n = SampSize, mean = 4.5, sd = .1)
b2 = rnorm(n = SampSize, mean = 4.5, sd = .4)
b3 = rnorm(n = SampSize, mean = 4.5, sd = .8)


K = A = 2;x0 = B = 0.001; a = C = .75
x = seq(0,6,.1)
y4 = A*exp(-B*exp(-C*x))
respon = function(x){K*exp(log(x0/K)*exp(-a*x))}
curve(respon(x), from = 0, to = 6)
points(respon(a1)~a1)

df = data.table(xs = c(a3, a2, a1, b3, b2, b1))
df[,y := respon(xs)]
df[,Varians2 := rep(c("a3","a2", "a1", "a3","a2", "a1"), each = SampSize)]
df[,Varians := rep(c("a3","a2", "a1", "b3","b2", "b1"), each = SampSize)]

piris <- ggplot(mapping = aes(col = Varians2)) +
  geom_segment(data = df[,.(y = mean(y), x0 = -3, x1 = mean(xs), Var = var(xs)), c("Varians","Varians2")], 
               aes(x = x0, y = y, xend = x1, yend = y), size = 1)+
  geom_point(data = df, aes(x = xs, y = y))+ 
  coord_cartesian(xlim = c(0, 6.5))+
  xlab("Observert variabel")+
  ylab("Målvariabel")
piris

