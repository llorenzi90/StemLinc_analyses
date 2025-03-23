library(ggplot2)
library(plotly)
sampnames=c("StemLinc",
            "Klimmeck")
title="All assembled genes"
cc_dat=as.data.frame(table(tracking_GL$StemLinc_cc,
                           tracking_GL$Klimmeck_cc))
cc_dat <- cc_dat %>%
  mutate(text = paste0("StemLinc: ", Var1, "\n", "Kli: ", Var2, "\n", "Value: ",Freq, "\n", "Frac:", round(Freq/sum(Freq),2)))
myplots <- htmltools::tagList()

p=ggplot(cc_dat, aes(Var1, Var2, fill= Freq,text=text)) +
  geom_tile() + ggtitle(title) + xlab(sampnames[1]) +
  ylab(sampnames[2])

myplots[[1]]=plotly::as_widget(plotly::ggplotly(p, tooltip="text"))
myplots


# DF <- data.frame(A = c(1:10), B = c(1:10), C = (1:10))
# myplots <- htmltools::tagList()
#
# for (i in 1: dim(DF)[2]) {
#   p = ggplot(DF, aes(x = DF[, i], y = C)) +
#     geom_point()
#   myplots[[i]] = plotly::as_widget(plotly::ggplotly(p))
# }
#
# myplots
