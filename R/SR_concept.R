
library(gsl)

## params

ra <- 3
rb <- 1.2e-4

ba <- 3
bb <- 3/1.4e4

## REF PTS

rmr <- ra/rb*exp(-1)
rsy <- (1 - lambert_W0(exp(1)/ra)) / rb
ruy <- 1 - lambert_W0(exp(1)/ra)

bmr <- ba/bb
bsy <- (ba/bb)*sqrt(1/ba)-(1/bb)
bsy <- (sqrt(ba)-1)/bb
buy <- 1 - sqrt(1/ba)

## plots

ss <- seq(0,1.2e4,10)

rr <- ra*ss/exp(rb*ss)

br <- ba*ss/(1 + bb*ss)

# dev.new(height=17/2.54, width=8/2.54)
pdf("SR_schematics.pdf", height=17/2.54+0.2, width=8/2.54+0.2)

## plots
layout(matrix(c(1,0,2),3,1), heights=lcm(c(8,1,8)), widths=lcm(8))
par(mai=c(0.3,0.3,0.1,0.1), omi=c(0,0,0,0.1))
# par(mfrow=c(2,1), mai=c(0.5,0.5,0.5,0.1), omi=c(0.1,0.1,0,0.1))

## Ricker
plot(ss, rr, type="n", xlim=range(ss), ylim=range(ss), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="L")
mtext(expression(italic(S[t])), 1, line=1, cex=1.2, at=12000)
mtext(expression(italic(R[t])), 2, line=0.5, cex=1.2, at=12000, las=1)
rmod <- expression((a)~Ricker:~italic(R[t])==frac(italic(alpha * S[t]),italic(e^{beta * S[t]})))
#mtext(rmod, 3, line=-1, cex=1, adj=0.1, las=1)
text(100, 11500, rmod, cex=1.1, pos=4, xpd=NA)
## 1:1
abline(a=0, b=1, col="gray")
#text(1.2e4, 1.2e4, "1:1", adj=c(1,0))
## R-S
lines(ss, rr, lwd=2)
#rmod <- expression(frac(italic(alpha * S[t]),italic(e^{beta * S[t]})))
#text(12000, 7500, rmod, pos=2)
## alpha
segments(0, 0, 1900, ra*1900, lty="dashed")
text(2000, ra*2000, expression(alpha), adj=c(0.5,0.5))
## MSY
segments(rsy,0,rsy,ra*rsy/exp(rb*rsy), lty="dashed")
text(rsy, -100, expression(frac(1-italic(W)(frac(italic(e),alpha)),beta)), pos=1, xpd=NA)
segments(par()$usr[1],ra*rsy/exp(rb*rsy),rsy,ra*rsy/exp(rb*rsy), lty="dashed")
text(0, ra*rsy/exp(rb*rsy), expression(italic(R)[MSY]), pos=2, xpd=NA)
## Rmax
#segments(0, rmr, 1/rb, rmr, lty="dashed")
#abline(h=rmr, lty="dashed")
#text(-50, rmr, expression(frac(alpha,beta)~italic(e)^{-1}), adj=c(1,0.5), xpd=NA)
## K
segments(0, log(ra)/rb, log(ra)/rb, log(ra)/rb, lty="dashed")
segments(log(ra)/rb, 0, log(ra)/rb, log(ra)/rb, lty="dashed")
text(log(ra)/rb, -100, expression(frac(log(alpha),beta)), pos=1, xpd=NA)
text(0, log(ra)/rb, expression(italic(K)), pos=2, xpd=NA)


## B-H
plot(ss, br, type="n", xlim=range(ss), ylim=range(ss), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="L")
mtext(expression(italic(S[t])), 1, line=1, cex=1.2, at=12000)
mtext(expression(italic(R[t])), 2, line=0.5, cex=1.2, at=12000, las=1)
bmod <- expression((b)~Beverton-Holt:~italic(R[t])==frac(italic(alpha * S[t]),1+italic(beta * S[t])))
#mtext(bmod, 3, line=-1, cex=1, adj=0.1, las=1)
text(100, 11500, bmod, cex=1.1, pos=4, xpd=NA)
## 1:1
abline(a=0, b=1, col="gray")
#text(1.2e4, 1.2e4, "1:1", adj=c(1,0))
## R-S
lines(ss, br, lwd=2)
#bmod <- expression(frac(italic(alpha * S[t]),1+italic(beta * S[t])))
#text(12000, 7500, bmod, pos=2)
## alpha
segments(0, 0, 1500, ba*1500, lty="dashed")
text(1600, ba*1600, expression(alpha), adj=c(0.5,0.5))
## MSY
segments(bsy,0,bsy,ba*bsy/(1 + bb*bsy), lty="dashed")
text(bsy, -100, expression(frac(root(alpha)-1,beta)), pos=1, xpd=NA)
segments(par()$usr[1],ba*bsy/(1 + bb*bsy),bsy,ba*bsy/(1 + bb*bsy), lty="dashed")
text(0, ba*bsy/(1 + bb*bsy), expression(italic(R)[MSY]), pos=2, xpd=NA)
## K
segments(0, (ba-1)/bb, (ba-1)/bb, (ba-1)/bb, lty="dashed")
segments((ba-1)/bb, 0, (ba-1)/bb, (ba-1)/bb, lty="dashed")
text((ba-1)/bb, -100, expression(frac(alpha-1,beta)), pos=1, xpd=NA)
text(0, (ba-1)/bb, expression(italic(K)), pos=2, xpd=NA)


dev.off()
