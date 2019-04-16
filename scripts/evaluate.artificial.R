library(stylo)
source("clustering.indices.R")

# available measures
measures = c("delta", "wurzburg")

# available number of MFWs
options.of.MFWs = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

# available linkage clustering.methods
clustering.methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

# get the corpora from your directory
corpora = list.dirs('.', recursive=FALSE)

# settings for artificial texts
artificial.text.size = 1000
no.of.artificial.texts = 10
artificial.corpus = list()
list.cutoff = 1000

# run the experiment for each of the available corpora
for (current.corpus.name in corpora){
  #current.corpus.name=corpora[2]
	# store initial variables
	init.variables = ls()	
	
	# set the current.corpus.name
	current.corpus.name =  sub("[^[:alpha:]]+", "", current.corpus.name)
	
	# set the location of the current corpus
	setwd(current.corpus.name)
	
		# loading the entire corpus, getting a frequency list
	texts = load.corpus.and.parse(files = "all", corpus.dir = "corpus", language = "Other")
	common.wordlist = make.frequency.list(texts)
	actual.classes = factor(gsub("_.*", "", names(texts)))

	for(current.class in unique(actual.classes)) {

		current.class.samples = list()

		class.members = names(texts)[gsub("_.*", "", names(texts)) == current.class]
		all.from.the.class = unlist(texts[c(class.members)])
		for(k in 1:no.of.artificial.texts) {
			current.class.samples[[k]] = sample(all.from.the.class, artificial.text.size, replace = TRUE)
		}

		current.freqs = make.table.of.frequencies(current.class.samples, common.wordlist[1:list.cutoff], absent.sensitive = FALSE)
		rownames(current.freqs) = paste(current.class, rownames(current.freqs), sep = "_")
		current.freqs = unclass(current.freqs)
		assign(paste("table_of_freqs_", current.class, sep = ""), current.freqs)
	}

	# joining frequency tables for particular classes
	frequency.table = c()
	for(p in grep("table_of_freqs_", ls(), value = TRUE)) {
		y = as.data.frame(mget(p))
		colnames(y) = colnames(current.freqs)
		frequency.table = rbind(frequency.table, y)
	}
	
	# get cluster names
	class.names = gsub("_.*", "", rownames(my.frequencies))
	
	# determine the number of unique clusters
	unique.clusters = length(unique(class.names))
	
	# perform the analysis for a given measure
	for (current.measure in measures){

		# perform stylo analysis to evaluate and save distances for the whole consensus loop over MFWs
		stylo.results = stylo(frequencies = my.frequencies, mfw.min = min(options.of.MFWs), mfw.max = max(options.of.MFWs), distance.measure =  current.measure, gui = FALSE, write.png = FALSE, save.analyzed.freqs = TRUE, save.distance.tables = TRUE)
	  
	 	for(current.MFWs in options.of.MFWs){
	    	file.rename(paste0("distance_table_", current.MFWs, "mfw_0c.txt"), paste0("distance_table_", current.MFWs, "mfw_0c_", current.measure, ".txt"))

			# load distance table if you have it already
			distance.file=paste0("distance_table_", current.MFWs,"mfw_0c_", current.measure,".txt")
			distance.table = t(read.table(distance.file, encoding="ANSI"))
			
			# Alternatively, perform stylo analysis
			#stylo.results = stylo(frequencies = my.frequencies, mfw.min = current.MFWs, mfw.max = current.MFWs, distance.measure =  current.measure, gui = FALSE, write.png=FALSE, save.analyzed.freqs = TRUE,save.distance.tables = TRUE)
			
			# run the experiment for each of the linkage clustering.methods available
			for (current.method in clustering.methods) {
				# assign current method
				# current.method = clustering.methods[1]

				# assign index variable name - to be used for saving results
				test.name = paste0(current.corpus.name, '_', current.measure, '_', current.MFWs, '_', current.method)

				# assign indices' file name
				text.file.name = paste0('L-I_', test.name, '.txt')
				
				# write current method and number of unique clusters to output file
				write(paste0('# Clustering method: ', current.method), file = text.file.name, append = TRUE, ncolumns=1000)
				write(paste0('# No. of ground truth clusters: ', unique.clusters),  file = text.file.name, append = TRUE, ncolumns=1000)
				iteration.info = paste0("# Clusters\tAMI\tARI\tNMI")
				write(iteration.info, file = text.file.name, append = TRUE, ncolumns=1000)

				# create a dendrogram for current method
				dendrogram = hclust(d = as.dist(distance.table), method = current.method)
				# alternatively, use (due to numerical precision of the distance table): dendrogram1 = hclust(d = as.dist(stylo.results$distance.table), method = current.method)
				
				# set max number of clusters aka number of texts
				max.clusters = as.integer(min(75, length(texts)))

				# evaluate quality of clustering for various numbers of clusters
				for(no.of.clusters in 2:max.clusters){
					# assign number of clusters to a variable
					# no.of.clusters = 3

					# evaluate clustering
					clustered = cutree(dendrogram, k = no.of.clusters)
					ground.truth = lapply(unique(class.names), function(class.ID){which(class.names == class.ID)})
					result = clustering.indices(ground.truth, as.vector(clustered))

					# save the number of clusters and the results of evaluation to a text file
					write(c(no.of.clusters,as.numeric(result)), file = text.file.name, append = TRUE, ncolumns=1000,sep="\t")

					# save the result to a variable
					assign(paste(test.name, no.of.clusters, sep=''), result)
				}
			}
		
	  # create network
		networkU = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT",network.type="undirected", distance.measure = current.measure, frequencies = my.frequencies, gui = FALSE)
		networkD = stylo(mfw.min = 100, mfw.max = 1000, analysis.type="BCT",network.type="directed", distance.measure = current.measure, frequencies = my.frequencies, gui = FALSE)
			
		# store network results
		#net.settings = paste0(current.corpus.name, '_', current.measure)			
		#assign(paste(net.settings, "_network", sep=''), network)

		# set file name
		current.file.name = paste0(current.corpus.name,'_', current.measure, '.RData')

		# store variables created in the study
		variables.to.be.saved = c("ground.truth","current.corpus.name","current.measure","unique.clusters","class.names","networkU","networkD")
    
		# save the results
		save(list = variables.to.be.saved, file = current.file.name, version = NULL)		
		 
		# run the tests for another measure
	  }
	}		 
	# run the tests for another corpus	 
	# move one directory up
	setwd("..")			 

}