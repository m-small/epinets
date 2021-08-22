
#new calculations of regional transmission, at request of the health department. August 2021

# load libraries
using LightGraphs, Plots, GraphPlot, Compose, Random, StatsBase, CSV, DataFrames
using JLD2, FileIO
include("ContactModels.jl")
include("EpiSim.jl")

#load data
wapop=CSV.read("all-wa.csv", DataFrame)
waloc=CSV.read("all-wa-pos.csv", DataFrame)
#########################################
#the following parameter is hardcoded here, even though it is used within ContactModels.jl just to ensure it is not forgotten....
global majorlocalities=58    #global scope intended. 
global brk=majorlocalities*2
#########################################

#configure and parse data
(nr,nc)=size(wapop)
townpop1=parse.(Int, replace.(wapop[1:2:brk,2], r","=> ""))
townpop2=parse.(Int, replace.(wapop[brk+1:end,2], r","=> ""))
popl=[townpop1; townpop2]
rowidx=[]
for i in majorlocalities+1:length(popl)
    if popl[i]!=0
    	global rowidx
        rowidx=[rowidx; i]
    end
end
rowidx1=[1:brk; rowidx.+majorlocalities]
rowidx2=[1:majorlocalities; rowidx]
#########
wapop2=wapop[rowidx1,:] #remove all empty settlements.
wapop=wapop2; #we'll just keep the one without the empty towns
waloc2=waloc[rowidx2,:]
waloc=waloc2
popl=popl[rowidx2]
posn=Array(waloc[:,3:4]);

wapop=Array(wapop)
#build (presumed) conectivity structure between localities
(nr,nc)=size(wapop)
nr=nr-majorlocalities
######################################
dist=Array{Int,2}(undef,nr,nr)
for i in 1:nr
    for j in 1:nr
        dist[i,j]=similaritytown(wapop,i,j,posn)
    end
end

# extract information from dataframe
#
#townpop=townpop[1:end-1]
#region and assign colouring
nodecolor = [colorant"lightseagreen", colorant"orange", colorant"red", colorant"green", colorant"yellow", colorant"brown", colorant"pink", colorant"purple", colorant"blue", colorant"violet", colorant"black"]
townlabel=[wapop[1:2:brk,3]; wapop[brk+1:end,3]]
#townlabel=townlabel[1:end-1]
regions=unique(townlabel)
regioni=Array{Int,1}(undef,nr)
for i in 1:nr
    regioni[i]=findall(x->x==1,townlabel[i].==regions)[1]
end
nodefillc = nodecolor[regioni]
#town names
townlabel=[wapop[1:2:brk,1]; wapop[brk+1:end,1]]
#townlabel=townlabel[1:end-1]




#'reasonable' parameters - alpha/beta ################################
epiparam=Dict()
epiparam["p0"]=0.2 #a guess - tuned to match observed data 
epiparam["q"]=1/7 #"up to" two weeks
epiparam["r0"]=1/14 #about two weeks for mild, 3-6 for severe
#########################################################
epiparam["nseeds"]=5 #probably too many, consider dropping. 
epiparam["pop"]=sum(popl) # =size(net)[1]



#'reasonable' parameters - delta - to override the preceeding parameters ################################
epiparam=Dict()
epiparam["p0"]=0.3 #a guess - tuned to match observed data, need more data  
epiparam["q"]=1/3 # 3-5 days - fairly robust
epiparam["r0"]=1/7 #about two weeks for mild, 3-6 for severe - bit fudgy
#########################################################
epiparam["nseeds"]=5 #probably too many, consider dropping. 
epiparam["pop"]=sum(popl) # =size(net)[1]


for ppt in [0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.1]
	for vac in [0.3 0.5 0.6 0.7 0.8 0.9 0.95]
	#build the transition network 
	#
	# in here we have parameters for probability of movement and assumed local connectivity structure
	#
	#pt=0.01 #probability of travel/proportion of population connecting multiple regions
	#
	# Various possible connectivity structures, listed below.
	#
	#net, transit, locale, popl, nedges_added = buildstate(wapop, x -> covidsafe(nomassmix(x),0.5), pt, true, true)
	#net, transit, locale, popl, nedges_added = buildstate(wapop, x -> nomassmix(x), pt, true, true)
	#net, transit, locale, popl, nedges_added = buildstate(wapop, x -> limitmix(x,500), pt, true, true, posn)
	#net, transit, locale, popl, nedges_added = buildstate(wapop, x -> fullmixing(x), pt, true, false)
	# note, covidsafe assumptions with probability p match vaccination with proportion p^2.
	net, transit, locale, popl, nedges_added = buildstate(wapop, x -> covidsafe(nomassmix(x),sqrt(vac)), ppt, true, true) #50% vaccination

	#the two boolean arguments are:
	#smalltown==true will bias toward travel to/from remote communities, otherwise travel will be predominantly to the largest places
	#richclub==true preferences "popular" (highly connected) nodes as travellers
	  

	#now for the MC simulation
	ndays=180
	nsims=500
	locpop=cumsum([1; popl],dims=1)
	St,Et,It,Rt = EpiSim.episimcom(net, locpop, epiparam, ndays, nsims)  #no change-point nsims simulations for ndays days
;
	
	#St,Et,It,Rt = 0,0,0,0

	fname="new_remote_p$(Int(ppt*1000))_TT_nomassmix$(Int(vac*100)).jld"

	#save the epidemic data
	save(fname,"St",St,"Et",Et,"It",It,"Rt",Rt,"net",net,"transit",transit,"locale",locale,"popl",popl)
end
end


exit()

