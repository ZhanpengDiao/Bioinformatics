#
# Smith-Waterman algorithm
#
# DESCRIPTION: 
# This file basically is Needleman-Wunsch. What I do to achieve the goal
# of Q5 is by adding additional limitions to the traceback step.
# The endpoint of an aligned sequence is point (i, j) with the maximal
# score at the last column of sequence Y on the F-matrix, while the
# start point is always placed at the first column behind the “-”. By
# doing that all elements of sequence Y will be kept.
#
plotting <- TRUE
#	Input sequences are encoded 1, 2, 3, 4 for A, C, G, T
#
# Q5-(1)
# X sequence: GAATCT
# Y sequence: CATT
#
X <- c(3, 1, 1, 4, 2, 4)
Y <- c(2, 1, 4, 4)
#
tStart <- proc.time()
#
NuclAcids <- c("A", "C", "G", "T")
#
seqX <- NuclAcids[X]
seqY <- NuclAcids[Y]
#
nX <- length(X)
nY <- length(Y)
#
#	Set up scoring matrix
# match = 1, mismatch = -1, gap penalty = -2
d <- 2  # gap penalty
s <- array(c(1,-1,-1,-1,
             -1, 1,-1,-1,
             -1,-1, 1,-1,
             -1,-1,-1, 1), dim=c(4, 4))
#
#	Set up the matrix F and pointers for traceback
# Caution: the indices of F[i, j] differ by 1 from the convention in the lecture notes
#
#	Initialisation
#
Fmatrix <- array(dim=c(nX + 1, nY + 1))
Fmatrix[,1] <- -(0:nX)*d
Fmatrix[1,] <- -(0:nY)*d
pointers <- array(rep(FALSE,nX*nY), dim=c(nX, nY, 3))

#
# Recursion: this algorithm uses Needleman-Wunsch
#
for(i in 1:nX){
    for(j in 1:nY){
        Fmatrix[i + 1, j + 1] <- max(Fmatrix[i, j] + s[X[i], Y[j]],
                                     Fmatrix[i, j + 1] - d,
                                     Fmatrix[i + 1, j] - d
        )
        pointers[i, j, ] <-  Fmatrix[i + 1, j + 1] == c(Fmatrix[i, j] + s[X[i], Y[j]],
                                                        Fmatrix[i, j + 1] - d,
                                                        Fmatrix[i + 1, j] - d
        )
    }
}
#
tEnd <- proc.time()
cat("\n Time taken =", (tEnd - tStart)[1], "seconds\n")
#
#	(Traceback not included in this program)
#
#	Draw a matrix showing paths
#
if(plotting){
    #
    plot.new()
    plot.window(xlim=c(0,nY), ylim=c(0,nX))
    text(rep(0:nY, times=(nX + 1)), rep(nX:0, each=(nY + 1)), labels=paste(t(Fmatrix)))
    mtext(c("-", seqY), at = 0:nY, side=3, line=2, col="blue", cex=1.5)
    mtext(c("-", seqX), at = nX:0, side=2, line=2, col="blue", cex=1.5, las=2)
    #
    #	Arrows along first row and column
    #
    arrows(x0=0, x1=0, y0=(nX:1) -1 + 0.2, y1=(nX:1) - 1 + 0.8, length=0.1)
    arrows(x0=(1:nY)-0.2, x1=(1:nY)-0.8, y0=nX, y1=nX, length=0.1)
    #
    #	Arrows in body of diagram
    #
    xFrom <- rep(1:nY, times=nX)
    yFrom <- rep((nX:1) - 1, each=nY)
    #
    xDiagArrows <- xFrom[t(pointers[,,1])]
    yDiagArrows <- yFrom[t(pointers[,,1])]
    arrows(x0=xDiagArrows-0.2, x1=xDiagArrows-0.8, y0=yDiagArrows+0.2, y1=yDiagArrows+0.8, length=0.1)
    #
    xUpArrows <- xFrom[t(pointers[,,2])]
    yUpArrows <- yFrom[t(pointers[,,2])]
    arrows(x0=xUpArrows, x1=xUpArrows, y0=yUpArrows + 0.2, y1=yUpArrows + 0.8, length=0.1)
    #
    xLeftArrows <- xFrom[t(pointers[,,3])]
    yLeftArrows <- yFrom[t(pointers[,,3])]
    arrows(x0=xLeftArrows-0.2, x1=xLeftArrows-0.8, y0=yLeftArrows, y1=yLeftArrows, length=0.1)
    #
}
#
