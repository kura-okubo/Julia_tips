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
N = round(Int, T/dt + 1) 

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

plot(p1, p2, layout = (2,1), legend=false)

#savefig("test.png")

#Do cross-correlation in time domain
