library(networkD3)
library(xlsx)
d=read.xlsx("C:/Users/cxu2/Documents/Book2.xlsx", header=T, sheetIndex=9)

nodes=unique(c(d$Source, d$Target))
nodes=data.frame("name"=nodes)

links=data.frame(source= match(d$Source, nodes$name) - 1,
			target= match(d$Target, nodes$name) - 1,
			value=rep(10,nrow(d)))

sankeyNetwork(Links = links, Nodes = nodes,
  Source = 'source', Target = 'target', Value = 'value',
  NodeID = 'name', sinksRight = FALSE, fontSize = 16)
