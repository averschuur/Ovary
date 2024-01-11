## merge betas UMCU and betas ovary

size <- as.data.frame(as.numeric(c("25", "14", "9", "10", "9", "15", "19", "10", "4", "9", "7", "11", "11")))
colnames(size) <- "size"

size %>%
  #group_by(size) %>%
  summarize(n = n(), avg = mean(size), sd = sd(size))

mean(size$size)
size %>%
  mean()
class(size)


pass <- as.data.frame(as.numeric(c("2", "2", "3", "3", "3", "3", "4", "4", "4", "5", "5")))
colnames(pass) <- "pass"

pass %>%
  #group_by(size) %>%
  summarize(n = n(), avg = mean(pass), sd = sd(pass))


pass <- as.data.frame(as.numeric(c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
                                   "2", "2", "2", 
                                   "3", "3",
                                   "4",
                                   "6")))
colnames(pass) <- "pass"


