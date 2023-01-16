library(fda)

scriptsDir <- system.file('scripts', package='fda')
Rscripts <- dir(scriptsDir, full.names=TRUE, pattern='R$')
fdarm <- grep('fdarm', Rscripts, value=TRUE)
chapters <- length(fdarm)
# NOTE: If R fails in any of these scripts,
# this for loop will not end normally,
# and the abnormal termination will be displayed:
for(ch in 1:chapters){
  cat('Running', fdarm[ch], '\n')
  invisible(source(fdarm[ch]))
}

Rscripts
Rscripts[15]

viewR

f <- fd(c(-1,2),create.bspline.basis(norder=2))

plot(f)

g <- f^2

fminus1 <- f^(-1)

plot(fminus1)
plot(g)


g
