using LightGraphs
using GraphPlot
using Plots

bamodel=barabasi_albert(10000, 3, 2)
lattice=LightGraphs.grid((100,100),periodic=true)
net=lattice
gplot(net)


pop=100^2
p=0.1
q=1/7
r=1/5
ndays=50

state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
state[1]=3
state[2:pop].=1

st=ex=fe=rm=Array{UInt64,1}(undef,ndays)
st[1]=count(x->x==1,state)
ex[1]=count(x->x==2,state)
fe[1]=count(x->x==3,state)
rm[1]=count(x->x==4,state)

gen=1
while gen<ndays
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
    gen += 1
    st[gen]=count(x->x==1,state)
    ex[gen]=count(x->x==2,state)
    fe[gen]=count(x->x==3,state)
    rm[gen]=count(x->x==4,state)
end


plot(1:ndays,[st,ex,fe,rm])
