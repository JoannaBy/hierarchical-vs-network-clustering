library(stylo)
setwd("German/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks directed
network_wurzburg = stylo(network = TRUE, network.type = "directed", mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

# preparing wurzburg networks usual style
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="German_networks.RData")

setwd("..")

setwd("French/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="French_networks.RData")

setwd("..")

setwd("English/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="English_networks.RData")

setwd("..")

setwd("Latin-Steffen/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="Latin-Steffen_networks.RData")

setwd("..")

setwd("TragedyComedy/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="TragedyComedy_networks.RData")

setwd("..")

setwd("Latin-NewOld/")

# loading the entire corpus, getting a frequency list
texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
frequency.list = make.frequency.list(texts)
my.frequencies = make.table.of.frequencies(texts, features=frequency.list, relative=FALSE)

# preparing wurzburg networks
network_wurzburg = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.wurzburg", frequencies = my.frequencies, gui = FALSE)

#preparing delta networks
network_delta = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT", distance.measure="dist.delta", frequencies = my.frequencies, gui = FALSE)

save(network_wurzburg, network_delta, file="Latin-NewOld_networks.RData")