#outbreak.jl

using LightGraphs, Plots, GraphPlot, Compose, Random
include("ContactModels.jl") #see  ContactModels.ipynb for a description of what this does
include("EpiSim.jl")
using CSV, DataFrames
using LightGraphs, SimpleWeightedGraphs, Statistics


#build the WA connectivity matrix
wapop=CSV.read("remote-wa.csv", DataFrame)

function similaritytown(df,i,j)
    # similaritytown(df,i,j)
    # Define a "similarity" between two locations. This code is an ad-hoc rule that is an
    # attempt to model the amount of movement between discrete regions. Dataframe df is
    # assumed to be in the Orlando format as illustrated above.
    bonus=100  #effectively, the larger the bonus, the more we include non-indigeneos travel.
    x1= Matrix(df[2*i-1:2*i,5:end])# convert(Matrix{Int}, df[2*i-1:2*i,5:end])
    x2= Matrix(df[2*j-1:2*j,5:end])#convert(Matrix{Int}, df[2*j-1:2*j,5:end])
#   languagesimilarity=sum(x1[:],dims=1)'*sum(x2[:],dims=1) #amalgamates the two different classifications of the boundary of a locale
    languagesimilarity=x1[:]'*x2[:]        #or not
    regionbonus= (df[2*i,3]==df[2*j,3])*bonus
end

dist=Array{Int,2}(undef,58,58) #distance, and eventually connetivity between towns
for i in 1:58
    for j in 1:58
        dist[i,j]=similaritytown(wapop,i,j)
    end
end

# extract information from dataframe
#
#population
townpop=wapopint=parse.(Int, replace.(wapop[1:2:end,2], r","=> ""))
townpop=townpop[1:end-1]
#region and assign colouring
nodecolor = [colorant"lightseagreen", colorant"orange", colorant"red", colorant"green", colorant"yellow", colorant"brown", colorant"pink", colorant"purple", colorant"blue", colorant"violet"]
townlabel=wapop[1:2:end,3]
townlabel=townlabel[1:end-1]
regions=unique(townlabel)
regioni=Array{Int,1}(undef,58)
for i in 1:58
    regioni[i]=findall(x->x==1,townlabel[i].==regions)[1]
end
nodefillc = nodecolor[regioni]
#town names
townlabel=wapop[1:2:end,1]
townlabel=townlabel[1:end-1]

function buildstate(statedata,netbuilder)
    # buildstate(statedata, netbuilder)
    ##############################################################################
    #statedata is a csv table loaded from the structure /format above (AKA the Orlando Format)
    #netbuilder is the intra-locale network construction rule (a function)
    ##############################################################################
    #
    #extract the data from CSV table
    locale=statedata[1:2:end, 1]
    popl=parse.(Int, replace.(statedata[1:2:end,2], r","=> ""))
    distrc=statedata[1:2:3, 3]
    npl=length(popl)
    ###############
  #  nppl=1000 # assume 1 in 1000 people moving, on average
    ###############
    net=SimpleGraph() #empty graph
    transit=Array{Int64,2}(undef,npl,npl) #number of transits beween locale[i] and locale[j]
    for (i,town) in enumerate(locale)
        println("Adding ",town," (population: ",popl[i],")")
        #add the intralocale links
        netadd=netbuilder(popl[i]) #use the preassigned method to add the new component
        net=blockdiag(net,netadd)
        #compute the amount of transit
        for j in 1:(i-1) #number of transit between here and every other previous part
        transit[i,j] = Int(floor(similaritytown(statedata,i,j)))
            transit[j,i]=transit[i,j]
        end
        transit[i,i]=0
    end
    #add the inter-locale links too --- could do it in the same (previous) loops, but this is clearer (I think)
    rr=[0; cumsum(popl)]
    nedges_add = 0
    for i in 1:npl
        for j in (i+1):npl
            edg1=rand(collect(rr[i]+1:rr[i+1]),transit[i,j])
            edg2=rand(collect(rr[j]+1:rr[j+1]),transit[i,j])
            for k in 1:transit[i,j]
                add_edge!(net,edg1[k],edg2[k]) #this is grossly inefficient
                nedges_add += 1
            end
        end
    end
    #all done
    println(" $nedges_add edges added")
    return net, transit, locale, popl, distrc, nedges_add
end


###############################3
# Now, this is the clever bit, we will build two connectivity matrices, one with mixing and the other with "lockdown"
locknet, transit, locale, popl, distrc = buildstate(wapop[1:end-1,:], x -> covidsafe(socialdist(x,0.6),0.8))
nomixnet, transit, locale, popl, distrc = buildstate(wapop[1:end-1,:], x -> covidsafe(fullmixing(x),0.5))
mixnet, transit, locale, popl, distrc = buildstate(wapop[1:end-1,:], x -> nomassmix(x))
fullmixnet, transit, locale, popl, distrc = buildstate(wapop[1:end-1,:], x -> fullmixing(x))

#'reasonable' parameters ################################
epiparam=Dict()
epiparam["p0"]=0.2 #a guess - tuned to match observed data
epiparam["p2"]=1/12 #revised infection rate with distancing measure
epiparam["q"]=1/7 #"up to" two weeks
epiparam["r0"]=1/14 #about two weeks for mild, 3-6 for severe
epiparam["r2"]=1/4 #revised removal rate (now due to testing and isolation)
#########################################################
epiparam["nseeds"]=1 #suppose there is just one person on the lose
epiparam["pop"]=size(locknet)[1]
#
nlock=5
ndays=90
nsims=200



#do something
#St,Et,It,Rt = EpiSim.episim(net, epiparam, ndays, nsims)  #no change-point nsims simulations for ndays days
infect=Array{Any,1}(undef,14)
for nlock in 1:14
    St,Et,It,Rt = EpiSim.episimdays(locknet, nomixnet, epiparam, nlock, ndays-nlock, nsims)  #no change-point nsims simulations for ndays days
    infect[nlock]=It+Rt
end


function plotqnt(x,y,col,labl,qnt::Float64=0.45)
    nt,ny=size(y)
    low=Array{Float64,1}(undef,nt)
    mid=Array{Float64,1}(undef,nt)
    hig=Array{Float64,1}(undef,nt)
    for i in 1:nt
        low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
    end
    plot!(x,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=labl)
end

plot()
usecol=[:red,:green,:orange,:yellow,:blue,:black,:red,:green,:orange,:yellow,:blue,:black,:red,:green,:orange,:yellow,:blue,:black,:red,:green,:orange,:yellow,:blue,:black,:red,:green,:orange,:yellow,:blue,:black]
for nlock in [1 4 7 14]
    plotqnt(1:ndays,infect[nlock],usecol[nlock],"$nlock-day lockdown")
end
