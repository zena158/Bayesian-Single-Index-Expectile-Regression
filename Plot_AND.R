# Define the Laplace PDF and Lasso penalty
AND_pdf <- function(x, xx, w, sigma2=1, tau) {
ss=2/(sqrt(sigma2*pi))*(sqrt(tau*(1-tau))/(sqrt(tau)+sqrt(1-tau))) 
dd= exp(-w*(x-xx)^2/sigma2)
pdf=ss*dd

}

ws <- function(tau, x, xx) {
w=rep(0, length(x))
for(i in 1:length(x)){
 if( x[i]>= xx) w1=tau 
 if( x[i] < xx) w1=1-tau 
w[i]=w1}
w
}

# Parameters
sigma2 <- 1
tau1 <- 0.5
tau2 <- 0.7
tau3 <- 0.9

xx=0
# Generate x values
x <- seq(-2, 2, length.out = 1000)
# Compute w values
w=ws(tau=tau1, x, xx)
pdf_values <- AND_pdf(x, xx=0, w, sigma2=1, tau=tau1)

# Plot AND PDF
pdf(file = "C:\\Users\\Al-Noor Center\\Desktop\\Expectile Regression\\Figures\\AND.pdf",   width = 8, height = 4) 
plot(x, pdf_values, type = "l", col = "red", lwd = 2,
     main = "", xlab = "y", ylab = "f(y)",
     ylim = c(0, max(pdf_values)))

w=ws(tau=tau2, x, xx)
pdf_values1 <- AND_pdf(x, xx=0, w, sigma2=1, tau=tau2)
lines(x, pdf_values1, col = "blue", lwd = 2, lty = 2)

w=ws(tau=tau3, x, xx)
pdf_values2 <- AND_pdf(x, xx=0, w, sigma2=1, tau=tau3)
lines(x, pdf_values2, col = "green", lwd = 2, lty = 4)


# Add a legend
legend("topright", legend = c(expression(paste(tau==0.50)), 
       expression(paste(tau==0.70)), expression(paste(tau==0.90))), 
       col = c("red", "blue","green"), lwd = 2,lty=c(1, 2, 4))
dev.off()
