#
#   Viterbi algorithm for HMMs
#
#   p++
#   matrix from "+" states to "+" states
Pin_in <- t(array(c(0.2, 0.2, 0.4, 0.1,
                    0.1, 0.3, 0.3, 0.2,
                    0.1, 0.3, 0.4, 0.1,
                    0.1, 0.3, 0.4, 0.1), dim=c(4, 4)))

#   p+-
#   matrix from "+" states to "-" states
Pin_out <- t(array(c(0.02, 0.02, 0.02, 0.02,
                     0.02, 0.02, 0.02, 0.02,
                     0.02, 0.02, 0.02, 0.02,
                     0.02, 0.02, 0.02, 0.02), dim=c(4, 4)))

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

Pi = cbind(rbind(Pin_in, Pout_in), rbind(Pin_out, Pout_out))

p_0j <- c(0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25)
p_i0 <- c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02)
states <- c("A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-")
emitted_symbols <- c("A", "C", "G", "T", "A", "C", "G", "T")
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
new_state = 0
new_symbol = 0
path <- c()
sequence <- c()
#
#	Initialisation
#
# let the system select the first state
new_state <- sample(seq(1, 8), size = 1, prob = p_0j)
# add the first emitted symbol
new_symbol <- sample(seq(1, 4), size = 1, prob = q[, new_state])
#
#	Recursion
#
pointer = 1
while(new_state != 9) {
    path <- c(path, new_state)
    sequence <- c(sequence, new_symbol)
    new_state <- sample(seq(1, 9), size = 1,  prob = c(Pi[path[pointer],], p_i0[path[pointer]]))
    new_symbol <- sample(seq(1, 4), size = 1, prob = q[, new_state])
    pointer <- pointer + 1
}
#
# COPY FROM Viterbi_algorithm.R
# used for check if the Viterbi path can match the true path
#
# let x be the new sequence generated above
x <- sequence
#
N <- dim(Pi)[1]			# Number of Markov states extracted from Pi
alphaSize <- dim(q)[1]	# alphabet size
L <- length(x)			# length of path
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
# compute the number of matches
#
matched <- path == piStar
match_rate = sum(matched) / length(matched)

# print out the result
print("Observed emitted symbols of the new chain:", quote = FALSE)
cat(c("B", emitted_symbols[path], "E"), "\n")
print("State path of the new chain:", quote = FALSE)
cat(c("B", states[path], "E"), "\n")
cat("\n Probability of the most probable path =", Prob,
    "\n Most probable path:  \n")
cat(c("B", states[piStar], "E"), "\n\n")
cat(sprintf("%.3f%% of states matched!", match_rate))
