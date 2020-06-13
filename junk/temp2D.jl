using LightGraphs
include("EpiSim.jl")
using Plots
using GraphPlot, Compose


function getbiggestbit(net)
    cc=connected_components(net)
    bigi=0
    bign=0
    for i in 1:length(cc)
        if length(cc[i])>bign 
            bign=length(cc[i])
            bigi=i
        end
    end
    return cc[bigi]
end

net=dorogovtsev_mendes(1000)

net=net[getbiggestbit(net)]
ndsz=[length(neighbors(net,vert)) for vert in vertices(net)]
gplot(net, nodefillc=colorant"black",nodesize=ndsz,NODESIZE=0.05)

function rewire!(net,p)
    #rewire edges of a graph with probability p
    vert=vertices(net)
    for edg in edges(net) #just do it the slow and literal way - we've got time
        if rand(1)[1]<p
            #rewire
            rem_edge!(net,edg)
            add_edge!(net,edg.src,rand(vert[vert.!=edg.src]))
        end
    end
    return net
end


net=LightGraphs.grid([100, 10], periodic=true)
gplot(net)

rewire!(net,0.01)
gplot(net)