
f <- function(par, I_K){
  y = 1 - exp(-(par/I_K))
}
I_K = 500

f1 <- function(par, I_S){
  y = (par/I_S) * exp(1 - (par/I_S))
}
I_S = 250

jpeg("./HABs-ABM/light_eqs_EFI.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f(x,I_K = I_K),from=0, to=1500,ylab='Photosynthesis rate',xlab = expression(paste("Light intensity"," (",mu,"mol ",m^-2," ",s^-1,")")),ylim = c(0,1), yaxt="n")
axis(2, las = 2)
curve(f1(x,I_S = I_S),from=0, to=1500, add = TRUE, lty = 2)
legend("bottomright", lty = c(1,2), legend = c("Monod: saturating light response","Steele: photoinhibition"),bty = "n")
dev.off()

jpeg("./plot_output/light_eq.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f1(x,I_S = I_S),from=0, to=1500,ylab='fI',xlab = expression(paste("Light intensity"," (",mu,"mol ",m^-2," ",s^-1,")")),ylim = c(0,1), yaxt="n")
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

f3 <- function(x, T_0, q){
  y = exp(-( (x - T_0) / q )^2)
}
T_0 = 22
q = 5

jpeg("./plot_output/temp_eq_Hellweger.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f3(x,T_0 = T_0, q = q),from=4, to=30,ylab='fT',xlab = 'Water temperature (°C)')
dev.off()

f4 <- function(x, R_resp, theta_resp){
  y = R_resp*theta_resp^(x-20)
}
R_resp = 0.08/1440
theta_resp = 1.08

jpeg("./plot_output/resp_eq.jpeg", res = 300, width = 5, height = 3.5, units = "in")
par(cex.lab = 1.5, mgp = c(2.7,1,0))
curve(f4(x,R_resp = R_resp, theta_resp = theta_resp),from=4, to=30,ylab='Probability of death',xlab = 'Water temperature (°C)')
dev.off()




