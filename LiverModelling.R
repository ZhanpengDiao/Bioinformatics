#
#	Liver modelling
#
	V_max <- 1 		# maximal rate, mmol/min
	F <- 1  		# flow rate, l/min
	K <- 0.1 	    # paramter of reaction, mmol/l

#
	# set the boundary for v
	V <- seq(0, 2, length=1000)

	# set the boundary for ci
	ci <- seq(0, 3, length=1000)

	# define the function of the equation 66
	# V = F * ci * (1 - exp(-(V_max - V)/(F * K)))
	conc_elim_equation <- function(ci, V) {
	    F * ci * (1 - exp(-(V_max - V)/(F * K))) - V
	}

	# caluate the result matrix
	eq_result_set <- outer(ci, V, conc_elim_equation)

	# draw the ci-V line
	contour(ci, V, eq_result_set, levels = 0,
	        xlab = "input substrate concentration (ci)", ylab = "total rate of elimination of substrate (V)")

# draw the curve for the homogeneous
# x = y - 1 on the note

	# define the function for the homogeneous
	# (K / ci) = (V_max / V) - 1
	homogeneous <- function(ci, V) {
	    (V_max / V) - 1 - (K / ci)
	}

	homo_result_set <- outer(ci, V, homogeneous)

	# draw the homogeneous line
	contour(ci, V, homo_result_set, levels = 0, labels = "homogeneous",
	        method = "edge", add = TRUE, lty = "dashed")

# draw the curve for the flow-limited regime
# rx = y on the note

	# define the function for the flow limited regime
	# (V_max / (F * K)) * (K / ci) = V_max / V
	flow_limited <- function(ci, V) {
	    ((F * ci) / V) - 1
	}

	# draw the flow limited line
	contour(ci, V, fl_result_set, levels = 0, labels = "flow-limited",
	        method = "edge", add = TRUE, lty = "dashed")

# draw the line for the changeover
# let y - 1 - x = y - rx

	# draw a vertical line for c*
	lines(rep(0.9, length=length(V)), c, type = 'l', lty = 2)
	text(1, 1.1, labels = c("c* = 0.9"))
