"""
Cross-correlation exercise

Problem:

▽:1 ---------> ▽:2
wave speed: c [m/s]
distance: d [m]
"""

using Plots

#------------------------------------------------#
cwave = 20		# wave speed between sensors [m/s]
dist = 1000	# distance between sensors [m]
f = 0.1		# frequency of signal [Hz]
dt = 0.1	# sampling time [s]
T = 100 	# time length of input signal [s]

#input parameters
t_init = 20	# initial time of signal at sensor 1

#------------------------------------------------#

# Making time series (pseudo inputs)
N = round(Int, T/dt + 1) # number of data point 

t = dt .* collect(0:N-1) # Time [s]

u1 = zeros(N)
u2 = zeros(N)

period = 1/f

if period < dt
	println("dt is too large.")
	return
end

init_id = round(Int, t_init/dt) + 1
delay = dist/cwave #[s]


for i = 1:round(Int, period/dt)
	u1[i+init_id] = sin(2*pi*f*dt*i)
	u2[i+init_id + round(Int,delay/dt)] = sin(2*pi*f*dt*i)
end

p1 = plot(t, u1, line=(:black, 1, :solid),
	marker = (:cross, 2, :green),
    ylabel = "Signal", 
    xlabel = "Time [s]",
    title = "Station 1",
    xlim = (0, T),
    ylim = (-1.5, 1.5)
    )

p2 = plot(t, u2, line=(:blue, 1, :solid),
	marker = (:cross, 2, :green),
    ylabel = "Signal", 
    xlabel = "Time [s]",
    title = "Station 2",
    xlim = (0, T),
    ylim = (-1.5, 1.5)
    )

#plot(p1, p2, layout = (2,1), legend=false)

#savefig("test.png")

#Do cross-correlation in time domain
#add zero for the length of N-1 to both sides of u2
#normalized by the number of data points N-k (k is the data number in zero padding)

Rcorr	= zeros(Float64, 2*N-1)
u2corr	= vcat(zeros(N-1), u2, zeros(N-1)) 
k 		= vcat(collect(N-1:-1:0), collect(1:N-1))

for n = 1:2*N-1

	#scale = 1/(N-k[n]) #scale by 1/(N-k)
	scale = 1/N

	for m = 1:N

		Rcorr[n] += scale * u1[m] * u2corr[n+m-1]

	end
end

#time axis is determined by the sampling rate
tcorr = dt .* collect(-(N-1):(N-1))

p3 = plot(tcorr, Rcorr, line=(:red, 1, :solid),
	marker = (:cross, 2, :green),
    ylabel = "Correlation", 
    xlabel = "Time [s]",
    title = "Cross-correlation between u1 and u2",
    xlim = (-T, T),
    ylim = (-1.2*maximum(abs.(Rcorr)), 1.2*maximum(abs.(Rcorr)))
    )

plot(p1, p2, p3, layout = (3,1), 
				 size = (800, 1200),
				 legend=false)

