N = 180
tempo = 1:N
erro = rnorm(N, sd= 3)

Xt = 5 + 3 * cos(2 * pi * 1/12 * tempo) + 4 * sin(2 * pi * 1/12 * tempo) + 
         2.5 * cos(2 * pi * 1/18 * tempo) + 1.5 * sin(2 * pi * 1/18 * tempo) + erro

plot(y = Xt, x = tempo, type = "l")
# plot(Xt ~ tempo)

a0 = mean(Xt)
Xt_2 = Xt - a0


# Fw = fft(Xt_2)
# Ij = Mod(Fw)^2
# 
# Ij_2 = Ij[1:(N/2)]
# 
# plot(Ij_2, type = "l")
# which(Ij_2 > 2000)


omega = 0:(N-1)/N
dj = c()
tempo = 1:N
dcos = c()
dsin = c()

for (j in 0:(N-1)) {
  dj[j+1] = sum(Xt *exp(-2 * pi * 1i * omega[j+1] * tempo)) / sqrt(N)
  dcos[j+1] = sum(Xt * cos(2 * pi *  omega[j+1] * tempo)) / sqrt(N)
  dsin[j+1] = sum(Xt * sin(2 * pi *  omega[j+1] * tempo)) / sqrt(N)
}

Ij = tempo*abs(dj)^2

round(Ij,4) == round(dcos^2 + dsin^2, 4)

Ij_2 = Ij[1:(N/2)]
plot(Ij_2, type = "l")
idx = which(Ij_2>30)

1 / omega[idx]


a12 = (2/sqrt(N)) * dcos[11]
b12 = (2/sqrt(N)) * dsin[11] 
a18 = (2/sqrt(N)) * dcos[16] 
b18 = (2/sqrt(N)) * dsin[16] 


Xt_est = a0 + a12*cos(2*pi*omega[11]*tempo) + b12*sin(2*pi*omega[11]*tempo) +
               a18*cos(2*pi*omega[16]*tempo) + b18*sin(2*pi*omega[16]*tempo)

plot(Xt, type = "l")
lines(Xt_est, col = 2)

xcos12 = cos(2*pi*omega[11]*tempo)
xsin12 = sin(2*pi*omega[11]*tempo)
xcos18 = cos(2*pi*omega[16]*tempo)
xsin18 = sin(2*pi*omega[16]*tempo)

xcos20 = cos(2*pi*omega[10]*tempo)
xsin20 = sin(2*pi*omega[10]*tempo)


modelo = lm(Xt ~ xcos12 + xsin12 + xcos18 + xsin18 + xcos20 + xsin20)
summary(modelo)

resid = Xt- Xt_est

plot(resid, type = "l")
acf(resid)






## Analise de série real




dados = read.table("https://www.ime.usp.br/~pam/m-cananeia.76.85")
serie = ts(dados$V1)
N = length(serie)

plot(serie)
Zt = serie-mean(serie)


dft = fft(Zt)/sqrt(N)
dj =  abs(dft)^2

plot(dj)
which(dj>100)

temp = 1:N
omega = 0:(N-1)/N
dcos = c()
dsin = c()
Ij = c()

for (j in 1:N) {
  Ij[j] = sum(Zt*exp(-2 * pi * 1i *omega[j]*temp)) / sqrt(N)
  dcos[j] = sum(Zt*cos(2*pi*omega[j]*temp)) / sqrt(N)
  dsin[j] = sum(Zt*sin(2*pi*omega[j]*temp)) / sqrt(N)
}


round(dcos^2+dsin^2, 4) == round(dj, 4)
round(dcos^2+dsin^2, 4) == round(abs(Ij)^2, 4)

a_0 = mean(serie)
a12 = (2/sqrt(N)) * dcos[11]
b12 = (2/sqrt(N)) * dsin[11]

temp = 1:N
Zt_est = a_0 + a12*cos(2*pi*omega[11]*temp) + b12*sin(2*pi*omega[11]*temp)


par(mfrow = c(1,1))
plot(serie)
lines(Zt_est, col = 2)


xcos = cos(2*pi*omega[11]*temp)
xsin = sin(2*pi*omega[11]*temp)

modelo = lm(serie ~ xcos + xsin)
summary(modelo)

plot(residuals(modelo), type = "l")
astsa::acf2(residuals(modelo))
shapiro.test(residuals(modelo)) # teste para verificar se os residuos são normais


