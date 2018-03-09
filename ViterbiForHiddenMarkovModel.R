# 
#   Viterbi algorithm for HMMs
#
#   Zhanpeng Diao - u5788688
#
#	Markov transition matrix:
#		Pi must be square,
#		rows of Pi plus p_i0 must sum to 1
#		elements of p_0j must sum to 1
#	Warning: it looks like the transpose of what's needed in the following lines of code
#		because of the way R converts vectors to arrays
#
#   p++
#   matrix from "+" states to "+" states
Pin_in <- t(array(c(0.2, 0.2, 0.4, 0.1,
                  0.1, 0.3, 0.3, 0.2,
                  0.1, 0.3, 0.4, 0.1,
                  0.1, 0.3, 0.4, 0.1), dim=c(4, 4)))
#
#   p+-
#   matrix from "+" states to "-" states
Pin_out <- t(array(c(0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02), dim=c(4, 4)))
#
#   p--
#   matrix from "-" states to "-" states
Pout_out <- t(array(c(0.3, 0.2, 0.2, 0.2,
                      0.3, 0.2, 0.1, 0.3,
                      0.2, 0.2, 0.3, 0.2,
                      0.1, 0.2, 0.3, 0.3), dim=c(4, 4)))

#   p-+
#   matrix from "-" states to "+" states
Pout_in <- t(array(c(0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02,
                   0.02, 0.02, 0.02, 0.02), dim=c(4, 4)))
#
# combine four separate 4*4 matrices together to make a 8*8 matrix
Pi = cbind(rbind(Pin_in, Pout_in), rbind(Pin_out, Pout_out))
#
p_0j <- c(0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25)
p_i0 <- c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02)
states <- c("A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-")
#
#	emission probabilities: q[a, i] is the probability of emitting symbol a from state i
#
q <- array(c( 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1,
              1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1), dim=c(4,8))
#
#	Observed sequence of symbols
#
# test sequence 1
# ATATATATCGCGCGCATATATAT
# x <- c(1, 4, 1, 4, 1, 4, 1, 4, 2, 3, 2, 3, 2, 3, 2, 1, 4, 1, 4, 1, 4, 1, 4)

# test sequence 2
# ATATATATCGCGCGCGATATATAT
# x <- c(1, 4, 1, 4, 1, 4, 1, 4, 2, 3, 2, 3, 2, 3, 2, 3, 1, 4, 1, 4, 1, 4, 1, 4)

# test sequence 3
# ATATATATCGCGCGCGATATATAG
 x <- c(1, 4, 1, 4, 1, 4, 1, 4, 2, 3, 2, 3, 2, 3, 2, 3, 1, 4, 1, 4, 1, 4, 1, 3)
#
#
N <- dim(Pi)[1]			# Number of Markov states extracted from Pi
alphaSize <- dim(q)[1]	# alphabet size
L <- length(x)			# length of path
#
#	Consistency checks on input parameters
#
if(dim(Pi)[2] != N)	stop("Pi is not square")
if(length(p_0j) != N) stop("p_0j is the wrong length")
if(length(p_i0) != N) stop("p_i0 is the wrong length")
if(!isTRUE(all.equal(rowSums(Pi) + p_i0, rep(1,N))))
    stop("rows of (Pi, p_i0) don't sum to 1")
if(sum(p_0j) != 1)	stop("elements of p_0j don't sum to 1")
if(length(states)!=N) stop("wrong number of state labels")
if(dim(q)[2]!=N)		stop("size of array q doesn't match number of Markov states")
if(!isTRUE(all.equal(colSums(q), rep(1,N))))
    stop("q probabilities don't sum to 1")
#
#	Set up arrays for v and pointers:
#		v[l,i] is probability of Viterbi path ending in state l with timestep i
#		ptr[i, l] points back to optimum state at timestep i-1 from timestep i
#		piStar is the optimum Markovian path
#
v <- array(dim=c(N, L))
ptr <- array(dim=c(L, N))
piStar <- array(dim=L)
#
#	Initialisation
#
v[,1] <- p_0j*q[x[1],]
ptr[1,] <- rep(0,N)
#
#	Recursion
#
for(i in 2:L){
    for(l in 1:N){
        v[l,i] <- q[x[i],l]*max(v[,i - 1]*Pi[,l])
        ptr[i,l] <- which.max(v[,i - 1]*Pi[,l])
    }
}
#
#	Termination
#
Prob <- max(v[,L]*p_i0)
piStar[L] <- which.max(v[,L]*p_i0)
#
#	Traceback
#
for(i in L:2){
    piStar[i - 1] <- ptr[i, piStar[i]]
}
#
cat("\n Probability of the most probable path =", Prob,
    "\n\n Most probable path:  \n")
cat(c("B", states[piStar], "E"), "\n\n")
#
