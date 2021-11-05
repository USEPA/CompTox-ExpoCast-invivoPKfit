# function to make line
myline <- function(x,m=1,b=0)
{
  y <- m*x + b
  
  return(y)
}

# a series of regularly spaced points
xpoints <- seq(0,10,by=0.5)
# set the random number generator seed so that we can get the same "random"
# numbers next time:
set.seed(1234)
# set the amount of error we want, the smaller this is the easier it will be to
# find the true values of m and b
ERROR <- 0.1
# a line with slope 1 and intercept 0 with random noise added:
ypoints <- (1-ERROR+runif(length(xpoints),-ERROR,ERROR))*myline(xpoints)

# turn it into a data.frame:
mydata <- data.frame(x=xpoints,y=ypoints)
head(mydata)

# a function to calculated the squared difference (this is similar to the 
# likelihood function):
mysquareddiff <- function(this.data,this.m,this.b)
{
  pred.y <- myline(this.data$x, m=this.m, b=this.b)
  squarediff <- (pred.y - this.data$y)^2
  return(sum(squarediff))
}
# Test it out for slope 0 and intecrept 0 (the wrong answer):
mysquareddiff(mydata,0,0)

# Now use the optimizer to see if we can find our original values of m and b
# despite adding ERROR
library(optimx)

# these values are where we start, if the initial valuees are too bad we won't
# find the answer:
initial.values <- c(0,0)
names(initial.values) <- c("m","b")

# Running optimx gets us close (m = 0.86, b = 0.129) but not perfect
# Less error should get us closer
optimx(initial.values, 
  function(x) mysquareddiff(mydata, this.m=x["m"], this.b=x["b"]))

# If we set the limits incorrectly (in this case, an upper limit of 0.5 on both 
# values when we know the slope is > 0.5) we get a bad answer:
optimx(initial.values, 
  function(x) mysquareddiff(mydata, this.m=x["m"], this.b=x["b"]),
  upper = c(0.5,0.5))
  

  