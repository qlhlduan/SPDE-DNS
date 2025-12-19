## ===============================
## Download and preprocess yield data
## ===============================

u <- url("https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt")
yield <- scan(u, skip = 14)
yield_mat <- NULL
for (i in 1:372)
  yield_mat <- rbind(yield_mat, yield[(19*(i-1)+1):(19*(i-1)+19)])
colname  <- c(1,3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)

colnames(yield_mat) <- c("Date", colname)

yield_mat <- as.data.frame(yield_mat)

yield_mat$Date <- as.Date(as.character(yield_mat$Date),format = "%Y%m%d")
