cd("/Users/michael/work/GitHub/epinets")
include("EpiSim.jl")

using CSV
using Plots
using LightGraphs
using JLD2, FileIO

cd("/Users/michael/work/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series")

usstates=["New York","California","Michigan","Illinois","Connecticut","Pennsylvania","Louisiana","Florida"]
uslocale=["New York","Los Angeles","Wayne","Cook","Fairfield","Philadelphia","Orleans","Miami-Dade"]

pops=[8398748 3990000 1749343 943332 1584064 343829 2716940]

z=Array{Any,1}(undef,length(usstates))
for (i,state) in enumerate(usstates)
    z[i]=EpiSim.getdataus(uslocale[i],state)
end

ddays=CSV.File("time_series_covid19_confirmed_US.csv")[1]
ddays=propertynames(ddays,4)[12:end]
ddays=String.(ddays)
ndays=length(ddays)

#pops=[426709 8089526 245869 5095100 1751693 534281 6594804 2621680]
cd("/Users/michael/work/GitHub/epinets")

#'reasonable' parameters
epiparam=Dict()
epiparam["p0"]=0.2 #a guess - tuned to match observed data
epiparam["p2"]=1/12 #revised infection rate with distancing measure
epiparam["q"]=1/7 #"up to" two weeks
epiparam["r0"]=1/14 #about two weeks for mild, 3-6 for severe
epiparam["r2"]=1/4 #revised removal rate (now due to testing and isolation)
epiparam["nseeds"]=5 #probably too many, consider dropping.

allrs=Dict()
totItps=Dict()
tpdays=Dict()
for (i,locale) in enumerate(uslocale)
    println("Working on ",locale)
    #get relevant data
    gridsize=Int(floor(sqrt(pops[i])))
    pop=gridsize^2
    epiparam["gridsize"]=gridsize
    epiparam["pop"]=pop
    y=z[i];
    #compute turning point
    ddt(z,zt)=count(z->z>0, z[1:zt])+ count(z->z<0, z[zt+1:end])
    ~,tpday=findmax([ddt(diff(diff(y)),nx) for nx in 1:(ndays-2)])
    #this is the turning point between exponential growth and decay. totItp total infections at day tpday
    while tpday>(ndays-2)
        tpday -=1
    end
    totItp=y[tpday+1]
    #build networks
    bamodel=barabasi_albert(gridsize^2, 3, 2)
    lattice=LightGraphs.grid((gridsize,gridsize),periodic=true)
    wattstrog95=watts_strogatz(gridsize^2, 4, 0.013)  #s=0.013 => 95% compliance
    wattstrog90=watts_strogatz(gridsize^2, 4, 0.026)  #s=0.026 => 90% compliance
    wattstrog80=watts_strogatz(gridsize^2, 4, 0.053)  #s=0.053 => 80% compliance
    #simulations
    St95,Et95,It95,Rt95=EpiSim.episim(bamodel,wattstrog95, epiparam, totItp, tpday+2, 90, 100)
    St90,Et90,It90,Rt90=EpiSim.episim(bamodel,wattstrog90, epiparam, totItp, tpday+2, 90, 100)
    St80,Et80,It80,Rt80=EpiSim.episim(bamodel,wattstrog80, epiparam, totItp, tpday+2, 90, 100)
    #store it all away
    rs=Array{Any,2}(undef,3,4)
    rs[1,1]=St95
    rs[1,2]=Et95
    rs[1,3]=It95
    rs[1,4]=Rt95
    rs[2,1]=St90
    rs[2,2]=Et90
    rs[2,3]=It90
    rs[2,4]=Rt90
    rs[3,1]=St80
    rs[3,2]=Et80
    rs[3,3]=It80
    rs[3,4]=Rt80
    allrs[locale]=deepcopy(rs)
    totItps[locale]=totItp
    tpdays[locale]=tpday
    filetitle="done_"*join(split(locale))
    @save filetitle allrs rs z totItps tpdays statename epiparam i
end
@save "done_US" allrs rs z totItps tpdays statename epiparam i
