###23 october 2023
# Supplementary material of manuscript entitled: 
# "Modeling evolutionary tempo and mode using formal morphological spaces 
# and Markov-chain principles" 

# List of required packages
library(ape)
library(castor)
library(grid)
library(igraph)
library(lattice)
library(markovchain)
library(phylobase)
library(phytools)
library(vegan)

# Read the NEXUS file with Lepidosaurian pyhlogeny and phalangeal formulas of terminals. 
# Set your working directory accordingly
setwd(choose.dir()) # Select the main folder where this supplementary material has been saved 
getwd() # Confirm that the working directory has been set correctly
squamata <- read.nexus("Data/Lepi.NEX")
PhFla <- readRDS("Data/PFphylotips.rds") #Counting of phalanges across digits 

# The first two instances of the approach consist of enumerating phenotypical options 
# and reconstructing phenotypes throughout the phylogenetic tree.
# The dataframe PhFla represents a comprehensive collection of phalangela formulas (PFs). 
# To investigate the evolutionary paths of PFs, we treat the number of phalanges in each digit as a discrete character
# and trace its evolutionary history onto a super-tree available at 
# Fontanarrosa et al. (2021) - DOI: 10.1002/jez.b.23026 
# The assignment of states to inner nodes is performed using maximum parsimony on a digit-by-digit basis. 
# Transition costs are considered proportional to the number of phalanges gained or lost. 
# Combining the inferred phalangeal counts for each digit, we reconstruct the PFs at inner nodes.

d1 <- PhFla[,1]
d2 <- PhFla[,2]
d3 <- PhFla[,3]
d4 <- PhFla[,4]
d5 <- PhFla[,5]
 
d1[d1 == "?"] <- 2 # Whenever a dubious character state codification exists at digit I, 
                   # it was filled with the plesiomorphic condition because of overall parsimony    
d1 <- as.numeric(as.character(d1))
d5[d5 == "?"] <- 3 # The same as above
d5 <- as.numeric(as.character(d5))

aux <- asr_max_parsimony(squamata, d1 + 1,  transition_costs="proportional") 
d1r <- apply(aux$ancestral_likelihoods, 1, which.max) - 1

aux <- asr_max_parsimony(squamata, d2 + 1,  transition_costs="proportional") 
d2r <- apply(aux$ancestral_likelihoods, 1, which.max) - 1

aux <- asr_max_parsimony(squamata, d3 + 1,  transition_costs="proportional") 
d3r <- apply(aux$ancestral_likelihoods, 1, which.max) - 1

aux <- asr_max_parsimony(squamata, d4 + 1,  transition_costs="proportional") 
d4r <- apply(aux$ancestral_likelihoods, 1, which.max) - 1

aux <- asr_max_parsimony(squamata, d5 + 1,  transition_costs="proportional") 
d5r <- apply(aux$ancestral_likelihoods, 1, which.max) - 1

flatip <- apply(cbind(d1, d2, d3, d4, d5), 1, paste, collapse = "")
flanode <- apply(cbind(d1r, d2r, d3r, d4r, d5r), 1, paste, collapse = "")
vectorflas <- c(flatip, flanode)

# Following, we frame the entire phylogeny within an evolutionary pseudo-time window ranging from 0 (root) to 1 (terminals). 
# We assign branch lengths using the procedure based on geometric means. This allocation is determined by considering the reciprocal 
# of path lengths from the parent node to all relevant terminals. In other words, the location of an inner node within the pseudo-time window,
# is influenced by the extent of divergence among its descendants and the number of branches involved. 
# Our method of 'relative positioning' for ancestral nodes is useful whenever lineage dating is not available. 
# Each internal node can be seen as a point in a polyline (linked segments) joining the immediate ancestor (parent node) and the terminal (tip).
# Such a polyline consists of as many legs (segments) as there are vertices from the internal node in question to the terminal.

squamata$edge.length[] <- 0
sq4 <- as(squamata, "phylo4") 
df <- as(sq4, "data.frame") 
desdeancdesc <- df[,3:2]  
ancestro <- unlist(ancestors(sq4, c(1:nrow(df)), "parent"))
lineages <- lapply(ancestors(sq4, 1:Ntip(squamata), "all"), function(x) rev(x)) 
for(i in 1:length(lineages))
   lineages[[i]] <- c(lineages[[i]], i) 
names(lineages) <- paste(squamata$tip.label, 1:Ntip(squamata))

nodolargo <- cbind(unlist(lineages), unlist(lapply(lineages, function(x) length(x):1))) 

# The geometric mean can also be expressed as the exponential of the arithmetic mean of logarithms.
frac <- exp(tapply(log(1/nodolargo[,2]), nodolargo[,1], sum) / tapply(nodolargo[,2], nodolargo[,1], length))
ancestriadirecta <- unlist(ancestors(sq4, 1:length(frac), "parent")) #Root node (i.e. id 673) does not have parents, so it appears as NA                           

#Traverse the tree with a breadth first search
posnodo <- rep(0, length(frac)) 
raiz <- which(is.na(ancestriadirecta))
queue <- which(ancestriadirecta == raiz)
branchlength <- c()
while(length(queue) > 0) {
   cual <- queue[1]
   queue <- c(queue, which(ancestriadirecta == cual))
   queue <- queue[-1]
   locancestro <- posnodo[ancestriadirecta[cual]]
   posnodo[cual] <- (1 - locancestro)*frac[cual] + locancestro
   branchlength <- c(branchlength, posnodo[cual] - locancestro)
} 

# Visualize the statistical distribution of branch lengths. Note the exponential tail 
hist(branchlength, freq = FALSE, breaks = 100, xlim = c(0, quantile(branchlength, 0.99)))

# Construction of the flow network from the frequency of evolutionary transitions between unique phenotypic conditions across the trajectory 
# from root to tips. Networks are collections of nodes and links (edges). In the case of flow networks, the links are directed, meaning they 
# have a specific source and target nodes.  In our analysis, we have categorized the phenotypic conditions (character states) of the nodes present in each
# evolutionary trajectory (ET). We focus on the unique conditions and consider them as the nodes of our network. The links in the network represent evolutionary
# transitions and are directed from the source node to the target node. These links are weighted based on the frequency of a particular state change observed in
# the ETs. 
nodosred <- unique(vectorflas)
redflujo <- matrix(0, length(nodosred), length(nodosred)) #Initialize with zeroes the adjacency matrix between nodes (unique phalangeal formulations)
for(i in 1:nrow(df)){
 link <- match(vectorflas[unlist(desdeancdesc[i,])], nodosred)  #vector of two elements, namely the source and the target
 redflujo[link[1], link[2]] <- redflujo[link[1], link[2]] + 1
} 
rownames(redflujo) <- colnames(redflujo) <- nodosred
diag(redflujo)[1:53] <-  diag(redflujo)[1:53] + 1
redflujo <- prop.table(redflujo, 1) #Transform the matrix into an stochastic one
redfinal <- graph_from_adjacency_matrix(redflujo, mode = "directed", 
                                        weighted = TRUE, diag = FALSE)
tkplot(redfinal)
# The next sentence retrieves the number id of the root node useful to customize the layout of the graph
match(vectorflas[673], nodosred)


# In the next lines we produce aligned trajectories using dummy nodes interspersed at regular intervals
# Here, we obtain another transtion matrix taking into account the relative positioning of internal nodes
trayectorias <- lapply(ancestors(sq4, 1:Ntip(squamata), "all"), function(x) rev(x)) 
for(i in 1:length(trayectorias))
   trayectorias[[i]] <- c(trayectorias[[i]], i) 
names(trayectorias) <- paste(squamata$tip.label, 1:Ntip(squamata))
trajflas <- lapply(trayectorias, function(x) vectorflas[x])

posictraj <- lapply(trayectorias, function(x) posnodo[x])
cortestree <- seq(0, 1, length.out = 50)
evoltraj <- matrix(0, length(posictraj), length(cortestree))
for(i in 1:nrow(evoltraj)) {
  cualesflas <- apply(abs(outer(cortestree, posictraj[[i]], "-")), 1, which.min) # dummy node is tagged with the formula of closest tree node
  evoltraj[i,] <- trajflas[[i]][cualesflas]
}

redflujo2 <- matrix(0, length(nodosred), length(nodosred)) #Initialize with zeroes the adjacency matrix between nodes (unique phalangeal formulations)
ndummies <- ncol(evoltraj)
for(i in 1:nrow(evoltraj)){
 links <- match(evoltraj[i,], nodosred)  #vector of two elements, namely the source and the target
 redflujo2[cbind(links[1:(ndummies - 1)], links[2:ndummies])] <- redflujo2[cbind(links[1:(ndummies - 1)], links[2:ndummies])] + 1
} 
rownames(redflujo2) <- colnames(redflujo2) <- nodosred
diag(redflujo2)[1:53] <-  diag(redflujo2)[1:53] + 1 
redflujo2 <- prop.table(redflujo2, 1) #Transform the matrix into an stochastic one

# The transtion matrix (also present in the RData file saved) is submitted to a Markov chain analysis.
# Among other parameters, we study the hitting time which is the relative required time to achieve a determined state from another one. 
# To calculate the hitting time, we consider each of the phenotypic states that are mapped onto the tree, with the plesiomorphic condition as the initial state.
# It can be used as a measure of the relative difficulty or complexity of the evolutionary paths leading to different states

# load("Data/redflujo.RData")
# redflujo <- redflujo2
dedos <- matrix(as.integer(unlist(strsplit(rownames(redflujo), ""))), ncol = 5, byrow = T)
distdedo <- as.matrix(dist(dedos, "manhattan"))
xy <- cmdscale(as.dist(distdedo))
plot(xy) #Project each point (phalangeal formula) onto a simplified two-dimensional space that captures the overall structure of distances among data points
distancias <- c()
acumtrans <- c()

# Installing the version of markovchain package used during the work
# install.packages("markovchain_0.8.0.tar.gz", repos = NULL, type = "source")
# Be sure you have installed Rtools 
simpleMc <- new("markovchain", states = rownames(redflujo), 
                transitionMatrix = redflujo)
idx <- 1:nrow(redflujo)
names(idx) <- rownames(redflujo)
for(i in 1:1000) {
 outs <- markovchainSequence(n = 100, markovchain = simpleMc, t0 = "23453",
                             include.t0 = TRUE)
 filtro <- idx[rle(outs)$values]
 if(length(filtro) == 1) next
 distancias <- c(distancias, mean(distdedo[cbind(filtro[1:(length(filtro) -1)], filtro[-1])]))
 acumtrans <- rbind(acumtrans, cbind(filtro[1:(length(filtro) -1)], filtro[-1])) 	
}
hist(distancias) #Statistical distribution of phenotypic distance between pairs of states involved in a simulated transition 

#Function for performing a random shuffling
shufflemtx <- function(mtx){
  stopifnot(isSymmetric(mtx))
  n <- nrow(mtx)
  rndn <- sample(1:n)
  out <- mtx[rndn, rndn]                     
  return(out)                     
}

shuffletransicion <- function(mtx){
  stopifnot(all(apply(mtx, 1, sum) == 1))
  out <- mtx
  for(i in 1:nrow(mtx))
        out[i,] <- sample(mtx[i,])
  return(out)                     
}

barplot(tapply(shuffletransicion(redflujo), shufflemtx(distdedo), sum))

#####
distanciasr <- c()
acumtransr <- c()
distdedo <- as.matrix(dist(dedos, "manhattan"))
azar <- sample(1:nrow(distdedo))
distdedor <- distdedo[azar, azar]
idx <- 1:nrow(redflujo)
names(idx) <- rownames(redflujo)
for(i in 1:1000) {
 outs <- markovchainSequence(n = 100, markovchain = simpleMc, t0 = "23453",
                             include.t0 = TRUE)
 filtro <- idx[rle(outs)$values]
 if(length(filtro) == 1) next
 #lines(xy[filtro,] + rnorm(length(filtro), 0.01, 0.05)) 
 acumtransr <- rbind(acumtransr, cbind(filtro[1:(length(filtro) -1)], filtro[-1])) 	
 distanciasr <- c(distanciasr, mean(distdedor[cbind(filtro[1:(length(filtro) -1)], filtro[-1])]))
}
hist(distanciasr)
mean(distanciasr <= mean(distancias))

###mean number of visits
meanNumVisits(simpleMc) -> out
maxlink <- pmax(out, t(out)) 
ej <- ifelse(maxlink > 0.5, 1, 0)

red <- igraph::graph_from_adjacency_matrix(redflujo, "directed", weighted = T)
maxlink <- pmax(out, t(out)) 
E(red)$weight <- 1/E(red)$weight

fptime <- matrix(0, nrow(redflujo), ncol(redflujo))
for(i in 1:nrow(redflujo)){
 firstPassage(simpleMc, states(simpleMc)[i], 10000) -> out
 fptime[i,] <- sweep(1:10000%*%out, 2, apply(out, 2, sum), "/")[1,]
}
fptime[is.na(fptime)] <- Inf
fptimeconstrained <- fptime*ifelse(redflujo > 0, 1, Inf)

# Next, we relate the number of steps required to achieve a given state from pseudo-time = 0, with  
# the degree of phenotypic distance between plesiomorphic condition at root and the target state.
# So, the initiation state is the PF 23453.
out <- firstPassage(simpleMc, states(simpleMc)[1], 10000) 
# 'out' represents a matrix of size 10000 rows x 56 columns (number of state characters).
# It shows the probability of the first time of passage in states to be exactly the number in the row starting from the root.
x <- apply(out, 2, which.max)
y <- as.matrix(dist(dedos, "manhattan"))[1,]


plot(x, y, ylab = "Manhattan distance from root", xlab = "Mode first passage time")
#identify(x, y)
#for(i in 1:ncol(out)){Sys.sleep(1);  barplot(out[1:100,i], main = colnames(out)[i])}

# Towards the synthetic map of tempo and mode.
# The last milestone of the approach refers to the combination of the phylogenetic (flow network) 
# and phenetic (morphospace)  information into a synthetic map of evolutionary change.
# The next levelplot aims to capture the time necessary to reach any phenotypic option starting from the root condition. 
modefpt <- apply(out, 2, which.max)
distdedo <- as.matrix(dist(dedos, "manhattan"))
xy <- cmdscale(as.dist(distdedo))

salida <- ordisurf(xy ~ modefpt)
plot(envfit(xy ~ dedos), add = T)

# Just for checking the correctness of the way the argument newdata should be used
# predict(salida, newdata = data.frame(x1 = xy[,1], x2 = xy[,2])) == predict(salida)

coordxrg <- range(xy[,1]) + c(-0.5, 0.5)
coordyrg <- range(xy[,2]) + c(-0.5, 0.5)
newxy <- expand.grid(seq(coordxrg[1], coordxrg[2], length.out = 10),
                     seq(coordyrg[1], coordyrg[2], length.out = 10))
predichos <- predict(salida, newdata = data.frame(x1 = newxy[,1], x2 = newxy[,2]))

#pdf("prueba.pdf" height = 4, width = 4, units = 'in')
#bitmap("Plotlevel2.tiff", height = 4, width = 4, units = 'in', res=300)
levelplot(predichos ~ newxy[,1] * newxy[,2], col.regions = topo.colors(12), 
          at = seq(0, ceiling(max(predichos)), length.out = 12), xlab = "Morphospace Axis 1", 
          ylab = "Morphospace Axis 2")
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text(expression("Mode fpt"), 0.2, -0.1, hjust=0.5, vjust=1)
trellis.unfocus()
#dev.off() #iun case of bitmap
#add data points corresponding to hand configurations either observed or phylogenetically reconstructed
#vector of positions, 1 bottom, 4 right. Ad hoc location after considering overlap between labels
pos <- rep(1, 56)
pos[44] <- 3
pos[48] <- 3
pos[53] <- 4
pos[3] <- 2
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(xy[,1], xy[,2], pch="*", col = "black", cex=2)
ltext(xy[,1], xy[,2], colnames(out), col="magenta", pos = pos, cex=0.7)
trellis.unfocus()
#dev.off()

