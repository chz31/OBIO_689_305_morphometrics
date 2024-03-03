data("scallops")
A <- gpagen(scallops$coorddata, scallops$curvslide, surfaces = scallops$surfslide)
A <- A$coords
ind <- factor(scallops$ind)
n <- dim(A)[3]; k <- dim(A)[2]; p <- dim(A)[1]; nind <- nlevels(ind); spec.names <- dimnames(A)[[3]]
object.sym = TRUE

land.pairs = scallops$land.pairs
npairs <- nrow(land.pairs); nl <- p-2*npairs

#Create reflections and relabel lms
A2 <- A
A2[land.pairs[,1],,] <- A[land.pairs[,2],,]
A2[land.pairs[,2],,] <- A[land.pairs[,1],,]
A2[,1,] <- A2[,1,]*(-1)
A <- array(c(A,A2), c(p,k, 2*n)) #get a landmark set of original A and mirrored A2)
ind <- factor(rep(ind,2)); side <- gl(2,n); 
ind

gpa.res <- gpagen(A, print.progress = TRUE) 
A <- gpa.res$coords  

#Symmetries
X.ind <- model.matrix(~ind + 0, data = as.data.frame(dat.shape[-1]))
symm.component <- arrayspecs(coef(lm.fit(X.ind, Y)),p,k) #equivalent to get an average between original and mirror for each pair of shapes

ind.names <- substr(dimnames(symm.component)[[3]], start=4,
                    stop=nchar(dimnames(symm.component)[[3]]))
ind.names
dimnames(symm.component)[[3]] <- ind.names
X.side <- model.matrix(~(side:ind) + 0, data = as.data.frame(dat.shape[-1]))

avg.side.symm <- coef(lm.fit(X.side, Y)) #The averages of orignal and reflections

#Asymmetry
#Asymmetry = reflection - original
asymm.component <- avg.side.symm[indsq,] - avg.side.symm[-indsq,]

#Grand mean
mn.shape <- mshape(A) 

# #Asymm shape = (reflection - original) + grand mean mn.shape
asymm.component <- simplify2array(lapply(1:n.ind, function(j) {
  t(matrix(asymm.component[j,],k,p)) + mn.shape
}))


###########################################
#get the original shape of specimen 1:
#Directly using the asymm residual calculated from (reflection - original) + symm.component
specimen1_test <- (asymm.component[, , 1] - mn.shape) + symm.component[, , 1]

#Calculating asymm residual as (relfection _ original)/2 
specimen1_test2 <- (asymm.component[, , 1] - mn.shape)/2 + symm.component[, , 1]

#The reflection of the specimen 1's shape
#A[, , 6]
A[1:5, , 6]
specimen1_test[1:5, ]
specimen1_test2[1:5, ]

