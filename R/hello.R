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


f1 <- function(x) {
  env <- environment()

  f2(y=x, env=env)

  f3()
}

f2 <- function(y, env) {
  z <- y + 1

  assign("z",z, envir = env)
}

f3 <- function() {
  print(z)
}

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

ff <- function(x) {
  print(x)
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

f <- function() {
  x <- 1
  g()
  g <- function() {
    print(x)
  }
  
  
}

library(data.table) 

dt1 <- data.table(x=1:3, y = 2:4)
dt2 <- data.table(z=5)

l <- list(dt1, dt2)

rbindlist(l,fill = TRUE, use.names = TRUE)

rbindlist(l,fill = TRUE, use.names = TRUE, idcol = "notice")

dt <- data.table(x=1:4)


m <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("a","b")))
m

preprocess_dt_center_all <- function(dt) {
  dt_copy <- copy(dt) # Avoid modifying original data.table
  categorical_cols <- names(dt_copy)[sapply(dt_copy, function(col) is.factor(col) || is.character(col))]
  
  # Convert categorical to numeric (one-hot encoding)
  for (col_name in categorical_cols) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL] # Remove original categorical column
  }
  
  # Center all columns
  for (col_name in names(dt_copy)) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}
preprocess_dt_center_covar <- function(dt, covar) {
  dt_copy <- copy(dt) # Avoid modifying original data.table
  covar_exists <- covar %in% names(dt_copy)
  if(!all(covar_exists)){
    stop(paste0("The following columns specified in covar do not exist in the data.table: ", paste(covar[!covar_exists], collapse = ", ")))
  }
  
  categorical_covar <- covar[sapply(dt_copy[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  # Convert categorical covar columns to numeric (one-hot encoding)
  for (col_name in categorical_covar) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL] # Remove original categorical column
    #Add new one hot encoded columns to the covar list, and remove the original column.
    covar <- c(covar[covar != col_name], names(dt_copy)[startsWith(names(dt_copy), paste0(col_name, "_"))])
  }
  
  # Center specified columns
  for (col_name in covar) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}
# Example Usage:
dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 4:7
)

dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 5:7
)

preprocessed_dt <- preprocess_dt_center_all(dt)
print(preprocessed_dt)

pdt <- preprocess_dt_center_covar(dt, covar = c("A","B"))

dt2 <- data.table(
  E = c(1,2,3,4),
  F = factor(c("red", "blue", "red", "green"))
)

preprocessed_dt2 <- preprocess_dt(dt2)
print(preprocessed_dt2)

preprocess_dt_no_copy <- function(dt, covar) { # No copy here
  covar_exists <- covar %in% names(dt)
  if(!all(covar_exists)){
    stop(paste0("The following columns specified in covar do not exist in the data.table: ", paste(covar[!covar_exists], collapse = ", ")))
  }
  categorical_covar <- covar[sapply(dt[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  for (col_name in categorical_covar) {
    col <- dt[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt[, (new_col_name) := as.integer(col == level)]
    }
    dt[, (col_name) := NULL]
    covar <- c(covar[covar != col_name], names(dt)[startsWith(names(dt), paste0(col_name, "_"))])
  }
  
  for (col_name in covar) {
    dt[, (col_name) := scale(dt[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt)
}

dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 4:7,
  E = 8:11
)

covar_cols <- c("A", "B", "E")
preprocessed_dt <- preprocess_dt_no_copy(dt, covar_cols)

print(preprocessed_dt)
print(dt) # The original 'dt' is now modified!

# Create a sample data.table
dt <- data.table(
  A = 1:5,
  B = 6:10,
  C = 11:15
)

# Take a subset
subset_dt <- dt[A > 2]

# Modify the subset
subset_dt[, B := B * 2]

# Check the original data.table
print(dt)
print(subset_dt)

# Example 2: using := with with=FALSE
dt2 <- data.table(
  A = 1:5,
  B = 6:10,
  C = 11:15
)

subset_dt2 <- dt2[A > 2]

subset_dt2[, (B) := B*2, with=FALSE]

print(dt2)
print(subset_dt2)

DT <- data.table(A = 1:5, B = 6:10)
DT2 <- DT[1:3]
DT3 <- DT[A>=3]
DT4 <- DT[,.(A,B)]
DT5 <- DT


df <- data.frame(A= 1:5, B= 7:11)
f <- function(d) {
  setDT(d)
  return(1)
}
f(df)

f <- function(d) {
  d <- as.data.table(d)
  return(class(d))
}

dt0 <- as.data.table(df)



# what if I have a column of constants?
# what if I have two column having the same columns?

df <- data.frame(y=rnorm(10), x1=2, x2 = rnorm(10))
fit <- lm(y~x1, data = df)

df <- data.frame(y=rnorm(10), x1=2, x2 = rnorm(10),x2=rnorm(10))
names(df) <- c("y","x1","x2","x2")
fit <- lm(y~x2, data = df)


num_ids <- DT[, uniqueN(get("A"))]

DT[1,1] <- 2

df <- data.frame(y=rnorm(10),x=rep(c("a","b"), each = 5))
df
fit <- lm(y~0+x, data = df)
summary(fit)
fit1 <- lm(y~x, data = df)
summary(fit1)


tryCatch({lm(y~x, data = dat)}, error = function(e) {
  return(e$message)
})
k <- tryCatch({lm(y~x, data = dat)
  lm(y~x, data = dat2)
  }, error = function(e) {return(list(e$message))})
for (i in 1:2) {
  tryCatch({lm(y~x, data = dat)}, error = function(e) {
    #return(e$message)
  })tryCatch({lm(y~x, data = dat)}, error = function(e) {})
}

counts <- c(10,15,3,19,6,9)
x <- c(0,0,0,1,1,1)
glm(counts ~ x, family = poisson(), control = list())

y=9
f <- function(x,y=NULL) {
  x+y
}
f(1)

f <- function(x){
  x+y
}

g <- function(x,y) {
  f <- function(x) {
    x+y
  }
  f(x)
}
g(1,2)
g <- function(x,y) {
  f(x)
  f <- function(x) {
    x+y
  }

}
g(1,2)



a <- if(F) 1
a


f <- function(x) {
  x+5
}
g <- function(f) {
  f(f)
}

library(lme4)
df <- data.frame(y=rnorm(10),x=rep(c("a","b"), each = 5))
df
fit <- lm(y~0+x, data = df)
summary(fit)
fixef(fit)
coef(fit)
vcov(fit)

f <- function(x, ...) {
  x+1
}


m1 <- matrix(c(1,-1), nrow = 1)
m2 <- matrix(1:2,nrow = 2)
v <- m1%*%m2

library(sandwich)
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## model fit and HC3 covariance
fm <- lm(y ~ x)
vcovHC(fm)
## usual covariance matrix
vcovHC(fm, type = "const")
vcov(fm)
sigma2 <- sum(residuals(lm(y ~ x))^2)/98
sigma2 * solve(crossprod(cbind(1, x)))

library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = F)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
(v1 <- vcov(fm1))
v2 <- vcov(fm1, correlation = TRUE)
v2
vcov(fm1, correlation = FALSE)
lme4::vcov.merMod(fm1)
v3 <- lme4::vcov.merMod(fm1, correlation = FALSE)

m <- fm1
sand <- sandwich(m, bread. = bread.lmerMod(m, full = TRUE),
                 meat. = meat(m, level = 2))

dat <- data.frame(y=rbinom(10,1,.5),x=rep(0:1,each=5))
dat <- data.frame(y=rbinom(10,1,.5),x=rnorm(10))
gm <- glm(y~x,data=dat, family = "binomial")

dat <- data.frame(y=rbinom(100,1,.5),x=rnorm(100))
gm <- glm(y~x,data=dat, family = "binomial")

# Number of observations
n <- 100  

# Generate categorical variables
x <- factor(sample(c("A", "B"), n, replace = TRUE))  
sample <- factor(sample(c("S1", "S2", "S3", "S4"), n, replace = TRUE))  

# Generate normally distributed response variable
y <- rnorm(n, mean = 0, sd = 1)  
x1 <- runif(n, min = 0, max = 1) # Uniformly distributed between 0 and 1
x2 <- runif(n, min = 0, max = 1) # Uniformly distributed between 0 and 1
x3 <- runif(n, min = 0, max = 1)
# Create data table
data <- data.table(y = y, x = x, sample = sample, x1 = x1, x2 = x2,x3=x3)
# Create data frame

m. <- lmer(y~x+x1+x2+x3+(1|sample), data = data)
m2 <- lm(y~x+x1+x2+x3, data = data)

dt1 <- data.table(x=1)
dt2 <- data.table(y=3)
dt3 <- NULL
l <- list(dt1, dt2, dt3)
l <- list(dt1, dt2, dt3, rbindlist(list()))
rbindlist(l[!is.na(l)], fill = T)
rbindlist(l, fill = T)

dt3 <- data.table(x=1, u=1:3)
dt3 <- data.table(x=1:2, u=1:3)




results <- list()
for (i in 1:10) {
  temp_dt <- data.table(x=1:5, y=i)
  temp_dt[,z:=x*y]
  results[[i]] <- temp_dt
}
result1 <- rbindlist(results)
result1

library(foreach)
library(doParallel)
registerDoParallel(cores = 2)
results_p <- foreach(i=1:10, .packages = "data.table") %dopar% {
  temp_dt <- data.table(x=1:5,y=i)
  temp_dt[, z:=x*y]
  temp_dt
}
stopImplicitCluster()
result2 <- rbindlist(results_p)
result2

h <- function(n) {
  setDTthreads(n)
}

dt <- data.table()
dt[, new_col := 1]
dt[, vec_col := list(list(c(1, 2, 3)))]
dt[, col2 := list(c(1:3))]

dt <- data.table(new_col = 1, vec_col = list(c(1, 2, 3)))

# Extract the vector into separate columns
vec_dt <- data.table(t(unlist(dt$vec_col)))

# Optionally rename the new columns
setnames(vec_dt, paste0("vec_col_", seq_len(ncol(vec_dt))))

# Combine with the original data.table (without vec_col)
dt_final <- cbind(dt[, .(new_col)], vec_dt)

print(dt_final)

centre <- function(x, type) {
  switch(type,
         mean = mean(x),
         median = median(x),
         trimmed = mean(x, trim = .1),
  stop("fku"))
}
x <- rcauchy(10)

f <- function(x) {
  x
  sys.function()
  }
f(1)
formals(f(1))

f <- function(x,y=1,z) {
  mc <- match.call()
  fmls <- formals(sys.function())
  required_args <- names(fmls)[sapply(fmls, is.symbol)]
  user_args <- names(mc)[-1]
  missing_required <- setdiff(required_args, user_args)
  
  missing_required
  #user_args
}
f()
f(1,2)
f(1,2,4)

g <- function(x) {
  1
}
g()


f2 <- function(f=sys.function()) {
  f
}
f2()
f3 <- function(){
  f2()
  
}
f3()

message("a",
        "b")

f <- function(x) {
  assertthat::assert_that(assertthat::is.flag(x))
  invisible(x)
}


f <- function(x) {
  if (x %% 3 == 0) {
    warning("hey")
  }
  x
}

g <- function(xs) {
  for (x in xs) {
    f(x)
  }
} 

xs <- 1:7
g(xs)

suppressWarnings(g(xs))

do_something <- function(x) {
  if (x %% 3 == 0) {
    warning("ahhh")
  }
  if (x %% 2 == 0) {
    warning("shit")
  }
  x
}
things <- 1:7

warnings_list <- list()

for (i in seq_along(things)) {
  this_warning <- NULL
  
  result <- withCallingHandlers({
    # Your code that might produce a warning
    risky_result <- do_something(things[[i]])
  }, warning = function(w) {
    this_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")  # suppress the warning
  })
  
  warnings_list[[i]] <- this_warning
}
warnings_list

results_list <- list()
warnings_list <- list()
errors_list <- list()

for (i in seq_along(things)) {
  #thing_name <- names(things)[i]  # or paste0("item_", i)
  thing_name <- paste0("things",i)
  
  # Initialize an empty list to collect warnings for this iteration
  these_warnings <- list()
  this_error <- NULL
  this_result <- NULL
  
  # Wrap the function call with tryCatch and warning handler
  this_result <- tryCatch(
    withCallingHandlers({
      # Call the potentially warning-raising function
      do_something(things[[i]])
    }, warning = function(w) {
      # Each warning gets added to the list
      these_warnings[[length(these_warnings) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")  # Muffles the warning output
    }),
    error = function(e) {
      # Capture any errors and set the result to NA
      this_error <<- conditionMessage(e)
      return(NA)
    }
  )
  
  # Store the result of the function call
  results_list[[thing_name]] <- this_result
  
  # If there were any warnings, store them in the warnings list
  if (length(these_warnings) > 0) {
    warnings_list[[thing_name]] <- these_warnings
  }
  
  # If there was an error, store the error message
  if (!is.null(this_error)) {
    errors_list[[thing_name]] <- this_error
  }
}
warnings_list




outer <- function() {
  x <- 1
  inner <- function() {
    x <<- 2
  }
  inner()
  x
}

outer()  # returns 1
x         # this will be 2 in the global environment now!

g <- function() {
  x <<- 2
}
h <- function() {
  g()
}
h()


func <- function(x) {
  if (x %% 3 == 0) {
    stop("oh no")
  }
  x
}

xs <- 1:7
re <- list()
e.list <- list()
for (x in xs) {
  re[[x]] <- tryCatch(func(x), error = function(e) {
    e.list[[x]] <<- e$message
    })
}
re
e.list

f1 <- function() {
  return()
}
f1()
f2 <- function() {
  NULL
}
f2()

tryCatch({
  warning("wwww")
  stop("eeee")
}, warning = function(w) {
  w$message
}, error = function(e) {
  e$message
})

f <- function() {
  x <- 1
}


f <- function(x) {
  l[[1]] <- x
}
l <- list()
f(1)

f <- function(x) {
  l[[1]] <<- x
}
l <- list()
f(1)

g <- function(x) {
  f <- function(x) {
    l[[1]] <<- x
  }
  l <- list()
  f(1)
}
g(1)


library(future.apply)
plan(multisession)

library(progressr)
handlers(global = TRUE)
oldh <- handlers("cli")

my_fcn <- function(xs) {
  p <- progressor(along = xs)
  y <- 1
  future_lapply(xs, function(x, ...) {
    Sys.sleep(6.0-x)
    #p(sprintf("x=%g", x))
    p(class="sticky")
    sqrt(x+y)
  })
}

my_fcn(1:5)

library(data.table)

# Create a data.table with double vectors
dt <- data.table(
  id = 1:3, 
  values = list(c(1/3, 2/3), c(pi, exp(1)), c(sqrt(2), log(2)))
)

# Save as CSV
fwrite(dt, "output.csv")

# Read the CSV file
dt_read <- fread("output.csv")

# Check structure
str(dt_read)

dt_read[, values := lapply(values, function(x) as.numeric(strsplit(x, "\\|")[[1]]))]

# Check the updated structure
str(dt_read)

dt_read

dt <- data.table(
  id = 1:3, 
  values = list(c(1/3, 2/3), c(pi, exp(1)), c(sqrt(2), log(2), 5))
)
# Find the maximum length of vectors
max_len <- max(lengths(dt$values))

# Expand list-column into separate columns
dt_wide <- dt[, paste0("value_", seq_len(max_len)) := transpose(values)]

# Remove the original 'values' column
dt_wide[, values := NULL]

# Check the result
print(dt_wide)

library(lme4)
library(detectseparation)
# Generate example data
set.seed(123)
data <- data.frame(
  y = rbinom(100, 1, 0.5),    # Binary outcome (0 or 1)
  x = factor(rep(1:2, each = 50)),  # Categorical predictor
  group = factor(rep(1:4, each = 25))  # Random effect grouping (10 groups)
)

# Fit logistic mixed-effects model
model <- glmer(y ~ x + (1 | group), data = data, family = binomial)

# Print model summary
summary(model)

glm(y ~ x, data = data, family = binomial, method = 'detect_separation')



set.seed(123)

# Generate example data
data <- data.frame(
  y = rbinom(100, 1, 0.5),  # Binary outcome
  x = factor(rep(1:2, each = 50)),  # Categorical predictor
  group = factor(rep(1:10, each = 10))  # Random effect grouping (10 groups)
)

# Force three clusters (e.g., groups 1, 3, and 5) to have y = 0
data$y[data$group %in% c(1, 3, 5)] <- 0

# Fit logistic mixed-effects model
model <- glmer(y ~ x + (1 | group), data = data, family = binomial)

# Print model summary
summary(model)

library(data.table)


dt <- data.table(x=1:3)
f <- function(d) {
  d[,y := 1:3]
  d
}
f(dt)

g <- function(d) {
  d[,y := 1:3]
  setDF(d)
}
g(dt)

iris

dd <- iris

h <- function(dd) {
  setDT(dd)
  dd[, xx := Sepal.Length * Sepal.Width]
  setDF(dd)
  setDF(dd)
  list(a=dd)
}
res <- h(dd = dd)
h <- function(dd) {
  setDT(dd)
  dd[, xx := Sepal.Length * Sepal.Width]
  setDF(dd)
}

res2 <- h(dd)

dd <- as.data.frame(iris)
class(dd)
h2 <- function(dd) {
  setDT(dd)
  dd[, xx := Sepal.Length * Sepal.Width]
  
}

res3 <- h2(dd)
head(dd)
class(dd)
