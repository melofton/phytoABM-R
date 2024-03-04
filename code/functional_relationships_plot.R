
f <- function(par, I_K){
  y = 1 - exp(-(par/I_K))
}
I_K = 500

f1 <- function(par, I_S){
  y = (par/I_S) * exp(1 - (par/I_S))
}
I_S = 300

jpeg("./HABs-ABM/light_eqs_EFI.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f(x,I_K = I_K),from=0, to=1500,ylab='Photosynthesis rate',xlab = expression(paste("Light intensity"," (",mu,"mol ",m^-2," ",s^-1,")")),ylim = c(0,1), yaxt="n")
axis(2, las = 2)
curve(f1(x,I_S = I_S),from=0, to=1500, add = TRUE, lty = 2)
legend("bottomright", lty = c(1,2), legend = c("Monod: saturating light response","Steele: photoinhibition"),bty = "n")
dev.off()

jpeg("./plot_output/light_eq.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f1(x,I_S = I_S),from=0, to=1500,ylab='Photosynthesis rate',xlab = expression(paste("Light intensity"," (",mu,"mol ",m^-2," ",s^-1,")")),ylim = c(0,1), yaxt="n")
axis(2, las = 2)
legend("topright", lty = c(1), legend = c("Steele: photoinhibition"),bty = "n")
dev.off()

f2 <- function(x,Tmin,Topt,Tmax){
  y = ((x - Tmin) / (Topt - Tmin)) *((Tmax - x) / (Tmax - Topt)) ^((Tmax - Topt) / (Topt - Tmin))
}
Tmin = 4
Topt = 23
Tmax = 30

jpeg("./plot_output/temp_eq_EFI.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f2(x,Tmin=Tmin,Topt=Topt,Tmax=Tmax),from=4, to=30,ylab='Growth rate',xlab = 'Water temperature (°C)')
dev.off()

