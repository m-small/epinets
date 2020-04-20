module EpiSim
    export episim, plotquantiles, getdata

    using LightGraphs, ProgressBars, CSV, Plots, Statistics

    #gridsize=1450
#pop=gridsize^2
#p0=0.2 #a guess - tuned to match observed data
#p2=1/12 #revised infection rate with distancing measures
#q=1/7 #"up to" two weeks
#r0=1/14 #about two weeks for mild, 3-6 for severe
#r2=1/4 #revised removal rate (now due to testing and isolation)
 #   ndays=120 #120 prediction (from patient(s) zere
 #   nsims=50 #50 runs
nseeds=5 #
    function episim(net1, net2, epiparam, feth::Int64=100, predays::Int64=40, postdays::Int64=80, nsims::Int64=50)
        #feth is the threhold to switch between growth and control
        pop=epiparam["pop"]
        p0=epiparam["p0"]
        p2=epiparam["p2"]
        q=epiparam["q"]
        r0=epiparam["r0"]
        r2=epiparam["r2"]
        nseeds=epiparam["nseeds"]

        state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
        state[1:pop].=1;
        state[rand((1:pop),(1,nseeds))].=2; #seeds should be exposed cases

        ndays=predays+postdays
        st=Array{UInt64,2}(undef,ndays,nsims)
        ex=Array{UInt64,2}(undef,ndays,nsims)
        fe=Array{UInt64,2}(undef,ndays,nsims)
        rm=Array{UInt64,2}(undef,ndays,nsims)
        st[1,1:nsims] .= count(x->x==1,state)
        ex[1,1:nsims] .= count(x->x==2,state)
        fe[1,1:nsims] .= count(x->x==3,state)
        rm[1,1:nsims] .= count(x->x==4,state)

        iter = ProgressBar(1:nsims)
        for j in iter
            net=net1
            #reinitialise the state vector
            state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
            state[1:pop].=1;
            state[rand((1:pop),(1,nseeds))].=2; #seeds should be exposed cases
            r=r0
            p=p0

            i=1
            notdoneyet=true
            while i<=(predays+postdays)

                for v in vertices(net)
                    if state[v]==1
                        for n in all_neighbors(net, v)
                            if state[n]==3 && rand(Float64).<p
                                state[v]=2   #susceptible becomes exposed
                            elseif state[n]==14 && rand(Float64).<p
                                state[v]=2   #susceptible becomes exposed
                            end
                        end
                    elseif state[v]==2 && rand(Float64).<q
                        state[v]=13 #exposed becomes infectious
                    elseif state[v]==3 && rand(Float64).<r
                        state[v]=14 #infectious becomes removed
                    end
                end
                state[state.>10]= state[state.>10] .- 10
                st[i,j] = count(x->x==1,state)
                ex[i,j] = count(x->x==2,state)
                fe[i,j] = count(x->x==3,state)
                rm[i,j] = count(x->x==4,state)

                #########################################
                #switch network after 300 infections
                inftot = pop-(st[i,j]+ex[i,j])
                if notdoneyet && inftot>feth
                    net=net2 #change network structure
                    r=r2     #change the removal rate
                    p=p2     #and the infection rate
                    #realign everything - this is clunky, should be a nicer way to achieve the same thing...
                    if i>predays
                            shft=i-predays
                            for k in 1:predays
                                st[k,j]=st[k+shft,j]
                                ex[k,j]=ex[k+shft,j]
                                fe[k,j]=fe[k+shft,j]
                                rm[k,j]=rm[k+shft,j]
                            end
                    elseif i<predays
                            shft=predays-i
                            for k in (predays-shft):-1:1
                                st[k+shft,j]=st[k,j]
                                ex[k+shft,j]=ex[k,j]
                                fe[k+shft,j]=fe[k,j]
                                rm[k+shft,j]=rm[k,j]
                            end
                            for k in 1:shft
                                st[k,j]=pop
                                ex[k,j]=0
                                fe[k,j]=0
                                rm[k,j]=0
                            end
                    end
                    i=predays
                    notdoneyet=false
                end
                #########################################
                i += 1
            end
        end

        return st, ex, fe, rm

    end


    function plotquantiles(y,col,labl,qnt::Float64=0.45)
        nt,ny=size(y)
        low=Array{Float64,1}(undef,nt)
        mid=Array{Float64,1}(undef,nt)
        hig=Array{Float64,1}(undef,nt)
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
        end
        plot!(1:nt,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=labl)
    end

    function getdata(country::String="Australia", state::String="Western Australia", file::String="time_series_covid19_confirmed_global.csv")

        x=[]
        for row in CSV.File(file)
            if row[2]==country
                if length(state)==0 && row[1]===missing
                    x=row
                elseif row[1]==state
                    x=row
                end
            end
        end

        if length(x)==0
            println("Error: "*country*", "*state*" not found.")
            return
        end

        z=Array{Int64,1}(undef,length(x)-4)

        for (index, value) in enumerate(x)
            if index>4
                z[index-4]=value
            end
        end

        return z
    end


    function getdataus(locality::String="Santa Fe", state::String="New Mexico", file::String="time_series_covid19_confirmed_US.csv")

        x=[]
        for row in CSV.File(file)
            if row[7]==state
                if length(locality)==0 && row[6]===missing
                    x=row
                elseif row[6]==locality
                    x=row
                end
            end
        end

        if length(x)==0
            println("Error: "*locality*", "*state*" not found.")
            return
        end

        z=Array{Int64,1}(undef,length(x)-11)

        for (index, value) in enumerate(x)
            if index>11
                z[index-11]=value
            end
        end

        return z
    end

end
