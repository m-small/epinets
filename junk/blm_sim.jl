cd("/Users/michael/work/GitHub/epinets")
include("EpiSim.jl")
using CSV
using Plots
using LightGraphs
using JLD2, FileIO
cd("/Users/michael/work/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series")
file="time_series_covid19_confirmed_global.csv"
allstates=["Australian Capital Territory" "New South Wales"  "Northern Territory" "Queensland" "South Australia" "Tasmania" "Victoria" "Western Australia"]
allcities=["Canberra" "Sydney" "Darwin" "Brisbane" "Adelaide" "Hobart" "Melbourne" "Perth"]
z=Array{Any,1}
for state in allstates
    z=push!(z,EpiSim.getdata("Australia",state))
end
ddays=CSV.File(file)[1]
ddays=propertynames(ddays,4)[5:end]
ddays=String.(ddays)
ndays=length(ddays)
pops=[426709 8089526 245869 5095100 1751693 534281 6594804 2621680]
cpop=[447457 4741874 132708 2326656 1315346 208324 4677157 2004696]
cinf=[107 3110 29 1061 440 228 1681 596]
#'reasonable' parameters
epiparam=Dict()
epiparam["p0"]=0.2 #a guess - tuned to match observed data 
epiparam["p2"]=1/12 #revised infection rate with distancing measure
epiparam["q"]=1/7 #"up to" two weeks
epiparam["r0"]=1/14 #about two weeks for mild, 3-6 for severe
epiparam["r2"]=1/4 #revised removal rate (now due to testing and isolation)
epiparam["nseeds"]=5 #probably too many, consider dropping.

for i in [2, 4, 7]
	epiparam["pop"]=pops[i]
	epiparam["gridsize"]=Int(floor(sqrt(pops[i])))
	y=z[i];
	ddt(z,zt)=count(z->z>0, z[1:zt])+ count(z->z<0, z[zt+1:end])
	~,tpday=findmax([ddt(diff(diff(y)),nx) for nx in 1:(ndays-2)])
	#this is the turning point between exponential growth and decay. totItp total infections at day tpday
	totItp=y[tpday+1]
	plot(1:tpday+1,y[1:tpday+1],lw=4,label="growth phase",title=allstates[i])
	plot!(tpday+1:ndays,y[tpday+1:ndays],lw=4,label="plateau")
	cd("/Users/michael/work/GitHub/epinets")
	gridsize=epiparam["gridsize"]
	pop=epiparam["pop"]
	bamodel=barabasi_albert(gridsize^2, 3, 2)
	#lattice=LightGraphs.grid((gridsize,gridsize),periodic=true)
	#wattstrog95=watts_strogatz(gridsize^2, 4, 0.013)  #s=0.013 => 95% compliance
	#wattstrog90=watts_strogatz(gridsize^2, 4, 0.026)  #s=0.026 => 90% compliance
	wattstrog80=watts_strogatz(gridsize^2, 4, 0.053)  #s=0.053 => 80% compliance
	#St95,Et95,It95,Rt95=EpiSim.episim(bamodel,wattstrog95, epiparam, totItp, tpday+2, 90, 100)
	#St90,Et90,It90,Rt90=EpiSim.episim(bamodel,wattstrog90, epiparam, totItp, tpday+2, 90, 100)
	
	th=cinf[i]

	St80,Et80,It80,Rt80=EpiSim.episim(bamodel,wattstrog80, epiparam, totItp, tpday+2, 90, 100)
	(Sp,Ep,Ip,Rp) = EpiSim.epipred(th, St80,Et80,It80,Rt80)


	iso=0.8 #was 0.95, I think 0.8 should be enough ... actually this parameter *may* significantly affect the results as it determines the susceptible pool size
	pop=epiparam["pop"]
	grd=Int(floor(sqrt(Int(floor(pop*(1-iso))))))
	isosize=pop-grd^2
	sg=SimpleGraph(isosize)
	ws=watts_strogatz(grd^2, 4, 0.026)
	isograph80=union(sg,ws)
	St9,Et9,It9,Rt9=EpiSim.episim3(bamodel,wattstrog80, isograph80, epiparam, totItp, tpday+2, 1, 180, 400)
	(Sp,Ep,Ip,Rp) = EpiSim.epipred(th, St9,Et9,It9,Rt9)

	RR = Rp[2:end,:] - Rp[1:end-1,:] ### daily increase of removed
	II = RR+Ip[2:end,:]-Ip[1:(end-1),:] ### daily increase of infected
	EE = II+Ep[2:end,:]-Ep[1:(end-1),:] ### daily increase of exposed
	ET = Ep[2:end,:] ### current number of exposed

	@save allcities[i] RR II EE ET

end