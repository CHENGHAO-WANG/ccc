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
