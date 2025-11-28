# Example data
x <- seq(-0.1, 0.1, length.out = 200)

y1 <- 0.8 + 0.03*sin(20*x)       # black
y2 <- exp(-(x/0.03)^2) * 0.45    # magenta
y3 <- exp(-(x/0.05)^2) * 0.35    # orange
y4 <- exp(-(x/0.06)^2) * 0.20    # brown
y5 <- rep(0.35, length(x))       # grey

# Plot
plot(x, y1, type="l", lwd=2, ylim=c(0,1),
     xlab="", ylab="", main="S1", bty="o", xaxt="n")

lines(x, y2, col="magenta", lwd=2)
lines(x, y3, col="orange",  lwd=2)
lines(x, y4, col="brown",   lwd=2)
lines(x, y5, col="grey50",  lwd=2)

axis(1, at=c(-0.1, 0, 0.1), labels=c("-0.1", "CEN", "+0.1"))