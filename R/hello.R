# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

#aaa

#bbb
a=1:5
b=1:6
c=1:4
Reduce(intersect,list(a,b,c))

m <- matrix(1:4, nrow = 2, byrow = T)
rownames(m) <- c("a","b")
colnames(m) <- c("c","d")

f <- function(eee) {
  assign("a",1, envir = eee)
  b <<- 1
  #print(a)
}

ff <- function() {
  eee <- environment()
  f(eee)
  print(a)
  ff0()
}

ff0 <- function() {
  print(a)
}

fff <- function(b) {
  print(b)
  rm(b)
  print(b)
}

ffff <- function(b) {
  fff(b)
  print(b)
}

gg <- function() {
  data("lrdb", envir=environment())
}

gg <- function() {
  data("lrdb")
}

ggg <- function() {
  print(environment())
}


h <- function(d=r) {
  print(d)
}

hh <- function() {
  d=1
  h(d)
}


hhh <- function() {
  r=1
  h()
}

hhhh <- function() {
  tt <- 1
  assign("w",0)
}

h5 <- function() {
  hhhh()
}
h6 <- function() {
  x=1
  x
}

f <- function(x) {
  print(x)
  ff(x)
}

ff <- function(y) {
  print(y)
}

suppressWarnings(h6())

h7 <- function() {
  assign("w",10)
  return(NULL)
}

A <- 1:5
B <- 3:6
setequal(intersect(A,B), A)

ee <- environment()

j <- function(ee) {
  rm(A, envir = ee)
}

mm <- function() {
  message("a \n")
  message("b")
}

a <- combn(1:10,2)
apply(a, 2, rev)


for (i in c
     ("a")) {
  print(i)
}


X=matrix(c(rep(0:1, each=10), rep(1:0, each=10)), ncol = 2)
X
y=X%*%c(1,2)+rnorm(20,0,.6)
y

fit <- lm(y~0+X)
fit
X2=matrix(c(rep(1, each=20), rep(1:0, each=10)), ncol = 2)
X2

fit2 <- lm(y~0+X2)
fit2

vcov(fit)
vcov(fit2)

f <- function() {
  for (j in 1:200) {
  Sys.sleep(0.01)
  cat("")
  progress(j, max.value = 200)
}
}

dat <- data.frame(gene=c("a","b","c"),c1=1:3, c2=2:4)
dat
expr.l <- t(dat[dat$gene=="c", -1])
expr.l
