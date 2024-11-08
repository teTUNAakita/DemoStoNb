# MakeFigs.R developed by Tetsuya Akita
#
# This code is provided for the purpose of reproducing the results presented in a manuscript that will be published,
# titled "Strong dependence of demographic stochasticity on a small effective number of breeders in iteroparous species,"
# for figure 2 and figure 3. See details in a manuscript.

### Figure 2a ###
# Relationship between Nb and exp[sigma_D^2]

pdf("~/burst.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
curve( 1 + 1/(x/2-1) , yaxs = "i",
       xlab = "",
       ylab = "",
       from=4, to=200, xlim = c(4,100), ylim = c(1,2.1), log = "xy", lwd=2, col = grey(0.5), xaxt = "n"
)
axis(1, at = c(4,10,20,50,100))
Nc = 200
curve( 1 + (1-x/Nc)/(x/2-1) ,xlab = "", ylab = "",
       from=4, to=Nc, ylim = c(1,2.1), log = "xy", lwd=2, add = TRUE, lty = 1
)
Nc = 100
curve( 1 + (1-x/Nc)/(x/2-1) ,xlab = "", ylab = "",
       from=4, to=Nc, ylim = c(1,2.1), log = "xy", lwd=2, add = TRUE, lty = 2
)
Nc = 50
curve( 1 + (1-x/Nc)/(x/2-1) ,xlab = "", ylab = "",
       from=4, to=Nc, ylim = c(1,2.1), log = "xy", lwd=2, add = TRUE, lty = 3
)
legend("topright", 
       legend = c("Approximation (Eq. 14)", 
                  paste0("Exact (Eq. 13, N = 200)"),
                  paste0("Exact (Eq. 13, N = 100)"),
                  paste0("Exact (Eq. 13, N = 50)")
       ),
       col = c(grey(0.5),"black","black","black"), 
       lty = c(1,1,2,3)
)
dev.off()

### Figure 2b ###
# Contour plot of the approximation accuracy in Eq.14 as a function of Nb and N

pdf("~/accuracy.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
n = 100 #
Nb = seq( 4, 1004, length = n )
N = seq( 4, 1004, length = n )
z <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    v_true = (1-Nb[i]/N[j])/(Nb[i]/2-1)
    v_approx = 1/(Nb[i]/2-1)
    if (Nb[i] > N[j]){
      z[i, j] = 100
    }else{
      z[i, j] = (v_approx - v_true)/v_true 
    }
  }
}
contour(Nb, N, z, drawlabels = F, levels=c(0.05, 0.1, 0.25), xlim = c(4,104), ylim = c(20,1004), log="y", lwd=2)
curve( 10*x ,xlab = "", ylab = "",
       from=4, to=100, log = "y", lwd=2, add = TRUE, lty = 2
)
legend("bottomright", 
       legend = c(expression(N[b]/N==0.1)),
       col = c("black"), 
       lty = c(3)
)
dev.off()

### Figure 3a ###
# Relationship between Nb and P_{E/ED} with approximation (Nb/N <<1)

pdf("~/missing.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
sigmaE = 0.2
curve( ( exp(sigmaE^2) - 1) / ( exp(sigmaE^2) * (1 + 1/(x/2-1)) - 1), yaxs = "i",
       xlab = "",#bquote("Effective breeding size,"~N[b]), 
       ylab = "",#bquote("Non-Poisson demographic factor, exp["~sigma[D]^2~"]"),
       from=4, to=200, xlim = c(4,200), ylim = c(0.01,1), log = "x", lty=1, xaxt = "n"
)
axis(1, at = c(4,10,20,50,100,200))
sigmaE = 0.4
curve( ( exp(sigmaE^2) - 1)  / ( (exp(sigmaE^2) * (1 + 1/(x/2-1))) - 1), xlab = "", ylab = "",
       from=4, to=200, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 2
)
sigmaE = 0.6
curve( ( exp(sigmaE^2) - 1)  / ( (exp(sigmaE^2) * (1 + 1/(x/2-1))) - 1), xlab = "", ylab = "",
       from=4, to=200, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 3
)
legend("bottomright", 
       legend = c(expression(sigma[E]==0.2),
                  expression(sigma[E]==0.4),
                  expression(sigma[E]==0.6)
                  ),
       col = c("black","black","black"), 
       lty = c(1,2,3)
)
dev.off()

### Figure 3b ###
# Relationship between Nb and P_{E/ED} with N = 50

pdf("~/missing50.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
sigmaE = 0.2
Nc=50
curve( ( exp(sigmaE^2) - 1) / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), yaxs = "i",
       xlab = "",#bquote("Effective breeding size,"~N[b]), 
       ylab = "",#bquote("Non-Poisson demographic factor, exp["~sigma[D]^2~"]"),
       from=4, to=50, xlim = c(4,50), ylim = c(0.01,1), log = "x", lty=1, xaxt = "n"
)
axis(1, at = c(4,10,20,50,100,200))
sigmaE = 0.4
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=50, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 2
)
sigmaE = 0.6
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=50, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 3
)
legend("bottomright", 
       legend = c(expression(sigma[E]==0.2),
                  expression(sigma[E]==0.4),
                  expression(sigma[E]==0.6)
       ),
       col = c("black","black","black"), 
       lty = c(1,2,3)
)
dev.off()

### Figure 3c ###
# Relationship between Nb and P_{E/ED} with N = 100

pdf("~/missing100.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
sigmaE = 0.2
Nc=100
curve( ( exp(sigmaE^2) - 1) / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), yaxs = "i",
       xlab = "",
       ylab = "",
       from=4, to=100, xlim = c(4,100), ylim = c(0.01,1), log = "x", lty=1, xaxt = "n"
)
axis(1, at = c(4,10,20,50,100,200))
sigmaE = 0.4
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=100, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 2
)
sigmaE = 0.6
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=100, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 3
)
legend("bottomright", 
       legend = c(expression(sigma[E]==0.2),
                  expression(sigma[E]==0.4),
                  expression(sigma[E]==0.6)
       ),
       col = c("black","black","black"), 
       lty = c(1,2,3)
)
dev.off()

### Figure 3d ###
# Relationship between Nb and P_{E/ED} with N = 200

pdf("~/missing200.pdf", width = 5, height = 5, pointsize = 12, bg = "white")
sigmaE = 0.2
Nc=200
curve( ( exp(sigmaE^2) - 1) / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), yaxs = "i",
       xlab = "",
       ylab = "",
       from=4, to=200, xlim = c(4,200), ylim = c(0.01,1), log = "x", lty=1, xaxt = "n"
)
axis(1, at = c(4,10,20,50,100,200))
sigmaE = 0.4
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=200, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 2
)
sigmaE = 0.6
curve( ( exp(sigmaE^2) - 1)  / ( exp(sigmaE^2) * ( 1+(1-x/Nc)/(x/2-1) ) - 1 ), xlab = "", ylab = "",
       from=4, to=200, ylim = c(0.01,1), log = "x", lwd=2, add = TRUE, lty = 3
)
legend("bottomright", 
       legend = c(expression(sigma[E]==0.2),
                  expression(sigma[E]==0.4),
                  expression(sigma[E]==0.6)
       ),
       col = c("black","black","black"), 
       lty = c(1,2,3)
)
dev.off()