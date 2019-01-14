using SAC, Plots

datanameE="E.BKKM..HNE.2012.001.sac"
datanameN="E.BKKM..HNN.2012.001.sac"
datanameU="E.BKKM..HNU.2012.001.sac"

yE = SAC.read("../KANTO/sac/event_2012_001/$datanameE")
yN = SAC.read("../KANTO/sac/event_2012_001/$datanameN")
yU = SAC.read("../KANTO/sac/event_2012_001/$datanameU")

tempE = yE.kstnm * "_" * yE.kcmpnm
tempN = yN.kstnm * "_" * yN.kcmpnm
tempU = yU.kstnm * "_" * yU.kcmpnm
label = [tempE, tempN, tempU]


yplotE = yE.t[1:end]
yplotN = yN.t[1:end]
yplotU = yU.t[1:end]

xplot = yE.delta .* collect(1:length(yplotE)) ./ 60 ./ 60

p1 = plot(xplot, yplotE, xlabel = "Time [h]",
					ylabel = "Acceleration",
					label=label[1],
					linecolor = :black,
					ylim = (-1e6,1e6)
					)

p2 = plot(xplot, yplotN, xlabel = "Time [h]",
					ylabel = "Acceleration",
					label=label[2],
					linecolor = :black,
					ylim = (-1e6,1e6)
					)

p3 = plot(xplot, yplotU, xlabel = "Time [h]",
					ylabel = "Acceleration",
					label=label[3],
					linecolor = :black,
					ylim = (-1e6,1e6)
					)

plot(p1, p2, p3, layout = (3,1),
				 size = (800, 1200))
savefig("./test1.png")