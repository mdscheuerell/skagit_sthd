
library(gsl)

## params
## Ricker
ra <- 3
rb <- 1.2e-4
## B-H
ba <- 3
bb <- 3/1.4e4

## ref pts
## Ricker
rmr <- ra/rb*exp(-1)
rsy <- (1 - lambert_W0(exp(1)/ra)) / rb
ruy <- 1 - lambert_W0(exp(1)/ra)
## B-H
bmr <- ba/bb
bsy <- (ba/bb)*sqrt(1/ba)-(1/bb)
bsy <- (sqrt(ba)-1)/bb
buy <- 1 - sqrt(1/ba)

## S-R curves
## spawners
ss <- seq(0,1.2e4,10)
## Ricker
rr <- ra*ss/exp(rb*ss)
## B-H
br <- ba*ss/(1 + bb*ss)

## plots
dev.new(height=6.4, width=3.25)
# pdf("SR_schematics_2.pdf", height=6.4, width=3.25)

layout(matrix(c(1,0,2),3,1), heights=lcm(c(3,0.3,3)*2.54), widths=lcm(3*2.54))
par(mai=c(0.4,0.4,0.2,0.2), omi=c(0,0,0,0.25))

## Ricker
plot(ss, rr, type="n", xlim=range(ss), ylim=range(ss), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="L")
mtext(expression(italic(S[t])), 1, line=1, cex=1.1, at=max(ss))
mtext(expression(italic(R[t])), 2, line=0.5, cex=1.1, at=max(ss), las=1)
rttl <- "(a) Ricker"
text(400, max(ss), rttl, cex=1.1, adj=c(0,1), xpd=NA)
## 1:1
abline(a=0, b=1, col="gray")
#text(1.2e4, 1.2e4, "1:1", adj=c(1,0))
## R-S
lines(ss, rr, lwd=2)
rmod <- expression(frac(italic(alpha * S[t]),italic(e^{beta * S[t]})))
text(12300, ra*max(ss)/exp(rb*max(ss)), rmod, adj=c(0,0.5), xpd=NA)
## alpha
segments(0, 0, 1900, ra*1900, lty="dashed")
text(2000, ra*2000, expression(alpha), adj=c(0.5,0.5))
## MSY
segments(rsy,0,rsy,ra*rsy/exp(rb*rsy), lty="dashed")
text(rsy, 0, expression(frac(1-italic(W)~bgroup("(",frac(italic(e),alpha),")"),beta)), adj=c(0.5,1.1), xpd=NA)
segments(par()$usr[1],ra*rsy/exp(rb*rsy),rsy,ra*rsy/exp(rb*rsy), lty="dashed")
text(0, ra*rsy/exp(rb*rsy), expression(italic(R)[MSY]), pos=2, xpd=NA)
## K
segments(0, log(ra)/rb, log(ra)/rb, log(ra)/rb, lty="dashed")
segments(log(ra)/rb, 0, log(ra)/rb, log(ra)/rb, lty="dashed")
text(log(ra)/rb, 0, expression(frac(log(alpha),beta)), adj=c(0.5,1.2), xpd=NA)
text(0, log(ra)/rb, expression(italic(K)), pos=2, xpd=NA)

## B-H
plot(ss, br, type="n", xlim=range(ss), ylim=range(ss), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n", bty="L")
mtext(expression(italic(S[t])), 1, line=1, cex=1.1, at=max(ss))
mtext(expression(italic(R[t])), 2, line=0.5, cex=1.1, at=max(ss), las=1)
bttl <- "(b) Beverton-Holt"
text(400, max(ss), bttl, cex=1.1, adj=c(0,1), xpd=NA)
## 1:1
abline(a=0, b=1, col="gray")
## R-S
lines(ss, br, lwd=2)
bmod <- expression(frac(italic(alpha * S[t]),1+italic(beta * S[t])))
text(max(ss)+300, ba*max(ss)/(1 + bb*max(ss)), bmod, adj=c(0,0.5), xpd=NA)
## alpha
segments(0, 0, 1500, ba*1500, lty="dashed")
text(1600, ba*1600, expression(alpha), adj=c(0.5,0.5))
## MSY
segments(bsy,0,bsy,ba*bsy/(1 + bb*bsy), lty="dashed")
text(bsy, 0, expression(frac(root(alpha)-1,beta)), adj=c(0.5,1.2), xpd=NA)
segments(par()$usr[1],ba*bsy/(1 + bb*bsy),bsy,ba*bsy/(1 + bb*bsy), lty="dashed")
text(0, ba*bsy/(1 + bb*bsy), expression(italic(R)[MSY]), pos=2, xpd=NA)
## K
segments(0, (ba-1)/bb, (ba-1)/bb, (ba-1)/bb, lty="dashed")
segments((ba-1)/bb, 0, (ba-1)/bb, (ba-1)/bb, lty="dashed")
text((ba-1)/bb, 0, expression(frac(alpha-1,beta)), adj=c(0.5,1.2), xpd=NA)
text(0, (ba-1)/bb, expression(italic(K)), pos=2, xpd=NA)

# dev.off()
