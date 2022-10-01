# library
library(ggplot2)
library(ggExtra)
# The mtcars dataset is proposed in R
head(mtcars)

# classic plot :
p <- ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, size=cyl)) +
  geom_point() +
  theme(legend.position="none")

  

# Set relative size of marginal plots (main plot 10x bigger than marginals)
p1 <- ggMarginal(p, type="histogram", size=10)

# Custom marginal plots:
p2 <- ggMarginal(p, type="histogram", fill = "slateblue", xparams = list(  bins=10))

# Show only marginal plot for x axis
p3 <- ggMarginal(p,  color="purple", size=4)
p3


curve(2.5*exp(3.2*x), from = 0, to = 1)
SampSize = 1000
BetaB = 4.2
BetaA = 1.5
x1 =inv.logit(rnorm(n = SampSize, sd = .125))
x2 = inv.logit(rnorm(n = SampSize, sd = .35))
x3 =inv.logit(rnorm(n = SampSize, sd = .5))
y0 = BetaA*exp(BetaB*x)
y1 = BetaA*exp(BetaB*x1)
y2 = BetaA*exp(BetaB*x2)
y3 = BetaA*exp(BetaB*x3)
B0 = 1
B1 = -4

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
               aes(x = x0, y = y, xend = x1, yend = y))+
  geom_point(data = df, aes(x = xs, y = y))+ 
  coord_cartesian(xlim = c(0, 8))
piris

