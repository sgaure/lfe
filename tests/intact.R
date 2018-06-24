library(lfe)
options(lfe.threads=1,digits=3,warn=1)
set.seed(42)
N <- 10000
x <- rnorm(N)
y <- rnorm(N)
time_id <- factor(sample(5,N,repl=TRUE))
group <- factor(sample(c('A','B','C','M','P','Q'),N,repl=TRUE))
data <- data.frame(x=x, y=y, time_id=time_id, group=group)

model_1<-felm(y~x|group:time_id,data=data)
model_2<-felm(y~x|time_id+group:time_id,data=data)
model_3<-felm(y~x|group:time_id+ time_id+group,data=data)

lm_1 <- lm(y ~ x + group:time_id, data=data)
lm_2 <- lm(y ~ x + time_id + group:time_id, data=data)
lm_3 <- lm(y ~ x + group:time_id + time_id + group,data=data)

message('felm 1'); print(model_1)
message('lm 1'); print(lm_1)
message('felm 2'); print(model_2)
message('lm 2'); print(lm_2)
message('felm 3'); print(model_3)
message('lm 3'); print(lm_3)

data[,'time+id'] <- data[,'time_id']
m <- felm(y~x|`time+id`+group + group:`time+id`,data=data)
print(getfe(model_3))
print(getfe(m))
