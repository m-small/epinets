#as before, because, why not....
using LightGraphs
using GraphPlot
using Plots
using ProgressBars
using Statistics
using LinearAlgebra
using JLD2

#epidemic simulation function
function episim(net1, net2, ndays::Int64=120, nsims::Int64=50, p2::Float64=1/12, r2::Float64=1/4)

    state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
    state[1:pop].=1;
    state[rand((1:pop),(1,nseeds))].=2; #seeds should be exposed cases

    st=Array{UInt64,2}(undef,ndays,nsims)
    ex=Array{UInt64,2}(undef,ndays,nsims)
    fe=Array{UInt64,2}(undef,ndays,nsims)
    rm=Array{UInt64,2}(undef,ndays,nsims)
    st[1,1:nsims] .= count(x->x==1,state)
    ex[1,1:nsims] .= count(x->x==2,state)
    fe[1,1:nsims] .= count(x->x==3,state)
    rm[1,1:nsims] .= count(x->x==4,state)

    print("p2=$p2; r2=$r2")

    iter = ProgressBar(1:nsims)
    for j in iter
        net=net1
        #reinitialise the state vector
        state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
        state[1:pop].=1;
        state[rand((1:pop),(1,nseeds))].=2; #seeds should be exposed cases
        r=r0
        p=p0

        for i in 1:ndays

            for v in vertices(net)
                if state[v]==1
                    for n in all_neighbors(net, v)
                        if state[n]==3 && rand(Float64).<p
                            state[v]=2   #susceptible becomes exposed
                        end
                    end
                elseif state[v]==2 && rand(Float64).<q
                    state[v]=3 #exposed becomes infectious
                elseif state[v]==3 && rand(Float64).<r
                    state[v]=4 #infectious becomes removed
                end
            end
            st[i,j] = count(x->x==1,state)
            ex[i,j] = count(x->x==2,state)
            fe[i,j] = count(x->x==3,state)
            rm[i,j] = count(x->x==4,state)

            #########################################
            #switch network after 300 infections
            if fe[i,j]>100
                net=net2 #change network structure
                r=r2     #change the removal rate
                p=p2     #and the infection rate
            end
            #########################################

        end
    end

    return st, ex, fe, rm

end


#requires global parameters
gridsize=1450
pop=gridsize^2
p0=0.2 #a guess - tuned to match observed data
p2=1/12 #revised infection rate with distancing measures
q=1/2 #"up to" two weeks
r0=1/14 #about two weeks for mild, 3-6 for severe
r2=1/4 #revised removal rate (now due to testing and isolation)
 #   ndays=120 #120 prediction (from patient(s) zere
 #   nsims=50 #50 runs
nseeds=5 #probably too many, consider dropping.

#and networks:
swp=0.01
bamodel=barabasi_albert(gridsize^2, 3, 2)
#lattice=LightGraphs.grid((gridsize,gridsize),periodic=true)
wattstrog=watts_strogatz(gridsize^2, 4, swp)


#pretty plot - I can't believe this is not built in
function plotquantiles(y,col,labl,qnt::Float64=0.3)
    nt,ny=size(y)
    low=Array{Float64,1}(undef,nt)
    mid=Array{Float64,1}(undef,nt)
    hig=Array{Float64,1}(undef,nt)
    for i in 1:nt
        low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
    end
    plot!(1:nt,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=labl)
end


#performance evaluation metrics
function recover50(A) #how many days does it take to eliminate infection in 50% of simulations?
    ndays,nsims=size(A)
    epigone=sum(x->x<1, A, dims=2) #number of zeros per column
    days=findall(x->x>(nsims/2),vec(epigone)) #days when epidemic is zero
    days=sort(days)
    if isempty(days)
        return -1   #never recovers
    end
    dayone=1
    while dayone==days[1] && length(days)>1
        dayone += 1
        deleteat!(days, 1)
    end
    if length(days)==1
        return 0 #never endemic
    else
        return days[1]
    end
end
#performance evaluation metrics

function meantotalinfected(A) #total infected
    ndays,nsims=size(A)
    return mean(A[end,:])
end

function meanmaxinfected(A) #maximumload
    return mean(maximum(A, dims=2)) #not sure about dims here
end

function parametersweep(swaps,p2s,r2s)
    ndays=300#0
    nsims=20#100

    ns=length(swaps)
    np=length(p2s)
    nr=length(r2s)
    mti=Array{Float64,3}(undef,ns,np,nr)
    r50=Array{Int64,3}(undef,ns,np,nr)
    mmi=Array{Float64,3}(undef,ns,np,nr)

    si=1
    bigiter=ProgressBar(swaps)
    for swp in bigiter
        wattstrog=watts_strogatz(gridsize^2, 4, swp)
        ri = 1
        for r2 in r2s
            pi=1
            for p2 in p2s
                St,Et,It,Rt = episim(bamodel,wattstrog,ndays,nsims,p2,r2) #switch between scale free and lattice (i.e. mixing to diffusion)
                r50[si,pi,ri] = recover50(Et+It) #irradication in 50% of simulations
                mti[si,pi,ri] = meantotalinfected(Rt) #total infections
                mmi[si,pi,ri] = meanmaxinfected(Rt) #highest total caseload
                pi += 1
            end
            ri += 1
        end
        print("Done swp=$(swp)")
        si += 1
    end

    return mti, mmi, r50
end

mdeg=4;
C=[1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]
swaps= 1 .- C.^(1 ./ mdeg)
p2s=[0.02:0.02:0.2;]
r2s=[0.2:0.05:0.5;]

mti,mmi,r50 = parametersweep(swaps,p2s,r2s)

@save("parameters_swept2020_th100")
