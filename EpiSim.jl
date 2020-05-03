module EpiSim

    export episim, plotquantiles, getdata
    using LightGraphs, ProgressBars, CSV, Plots, Statistics

    function epistep!(state, net, p, q, r)
        #internal function for episim - can also be useful from REPL
        #perform a single step of the epidemic propagation process
        #loop over the nodes, updating the infection states
        #state is a vector and net is a contact network with commensurate dimensions
        # p, q and r are the epidemic parameters
        for v in vertices(net)
            if state[v] == 1
                for n in all_neighbors(net, v)
                    if state[n] == 3 && rand(Float64) .< p
                        state[v] = 2   #susceptible becomes exposed
                    elseif state[n] == 14 && rand(Float64) .< p
                        state[v] = 2   #susceptible becomes exposed
                    end
                end
            elseif state[v] == 2 && rand(Float64) .< q
                state[v] = 13 #exposed becomes infectious (state 13 is newly infectious)
            elseif state[v] == 3 && rand(Float64) .< r
                state[v] = 14 #infectious becomes removed (state 14 is newly removed)
            end
        end
        #update the new infected and removed nodes
        state[state.>10] = state[state.>10] .- 10
    end #of epistep!

    function epirealign!(st, ex, fe, rm, i, j, predays, pop)
    #internal function for episim - the following performs realignment of the
    #single predicted trajectory so that the changepoint is correctly occuring at
    #time predays, arrays st, ex, fe, and rm are updated in the process
        if i>predays #our simulation was a bit slow
                shft=i-predays
                for k in 1:predays #so slide it back to re-align
                    st[k,j] = st[k+shft,j]
                    ex[k,j] = ex[k+shft,j]
                    fe[k,j] = fe[k+shft,j]
                    rm[k,j] = rm[k+shft,j]
                end
        elseif i<predays #our simulation was quicker than expected
                shft=predays-i
                for k in (predays-shft):-1:1 #so push it forwad
                    st[k+shft,j] = st[k,j]
                    ex[k+shft,j] = ex[k,j]
                    fe[k+shft,j] = fe[k,j]
                    rm[k+shft,j] = rm[k,j]
                end
                for k in 1:shft #and add some zeros
                    st[k,j] = pop
                    ex[k,j] = 0
                    fe[k,j] = 0
                    rm[k,j] = 0
                end
        end
    end #of epirealign!

    function episim(net1, net2, epiparam, feth::Int64=100, predays::Int64=40, postdays::Int64=80, nsims::Int64=50)
        #single changepoint model
        #feth is the threhold to switch between growth and control - replace net1 with net2

        #epidemic parameters
        pop=epiparam["pop"]         #population - assumed to be a square number
        p0=epiparam["p0"]           #rate p for first phase
        p2=epiparam["p2"]           #rate p for the second phase
        q=epiparam["q"]             #rate q (latency) fixed throughout
        r0=epiparam["r0"]           #rate r for first phase
        r2=epiparam["r2"]           #rate r for second phase
        nseeds=epiparam["nseeds"]   #number of infected nodes to start with

        #initialise and get set to track the various tallies
        ndays=predays+postdays
        st=Array{UInt64,2}(undef,ndays,nsims)
        ex=Array{UInt64,2}(undef,ndays,nsims)
        fe=Array{UInt64,2}(undef,ndays,nsims)
        rm=Array{UInt64,2}(undef,ndays,nsims)

        iter = ProgressBar(1:nsims) #everyone loves a good progress bar

        #main iteratarion loop nsims simulations
        for j in iter
            #active contact network is net1 to start
            net=net1
            #initilise the rate parameters
            r=r0
            p=p0

            #reinitialise the state vector
            state=Array{Int8,2}(undef,1,pop) #this is a bit wasteful, there must be a better categorical way to do this...
            state[1:pop].=1;
            state[rand((1:pop),(1,nseeds))].=2; #seeds should be exposed cases

            i=1
            notdoneyet=true
            #loop for a single epidemic time course (total duration predays+postdays)
            while i<=(predays+postdays)

                #updating the infection states
                epistep!(state, net, p, q, r)

                #count the respective totals
                st[i,j] = count(x->x==1, state)
                ex[i,j] = count(x->x==2, state)
                fe[i,j] = count(x->x==3, state)
                rm[i,j] = count(x->x==4, state)

                ##################################################################################
                #switch network after 300 infections - this is the change-point cludge
                ##################################################################################
                inftot = pop-(st[i,j]+ex[i,j])
                ##################################################################################
                # our criteria for changepoint is that the total number of identified infections
                # (i.e. infected + removed, which we compute ad inftot) exceeds feth
                # once achieved, the transmission model switches to net2, 2r and p2.
                # but, we also need to shif the time history so that in the epidemic curves
                # st, ex, fe, and rm these transition all occur at the same time - day predays
                ##################################################################################
                if notdoneyet && inftot>feth   #notdoneyet - we only want to do this ONCE
                    net=net2 #change network structure
                    r=r2     #change the removal rate
                    p=p2     #and the infection rate
                    #realign everything - this is clunky, should be a nicer way to achieve the same thing...
                    epirealign!(st, ex, fe, rm, i, j, predays, pop)
                    i=predays           #counter set to done predays.
                    notdoneyet=false    #we don't want to come back here again
                end
                ##################################################################################
                ##################################################################################
                ##################################################################################

                i += 1  #done one more day

            end
        end

        return st, ex, fe, rm

    end #episim

    function plotquantiles(y,col,qnt::Float64=0.45)   #no nk overloading the method should work...
        #plot an infection curve with qnt*2 confidence intervals - no label
        nt,ny=size(y)
        low=Array{Float64,1}(undef,nt)
        mid=Array{Float64,1}(undef,nt)
        hig=Array{Float64,1}(undef,nt)
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
        end
        plot!(1:nt,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=false) #updates current plot
    end #of plotquantiles(y,col,qnt::Float64)

    function plotquantiles(y,col,labl,qnt::Float64=0.45)
        #plot an infection curve with qnt*2 confidence intervals - and label
        #really suprised that julia doesn't have something to do this already...
        nt,ny=size(y)
        low=Array{Float64,1}(undef,nt)
        mid=Array{Float64,1}(undef,nt)
        hig=Array{Float64,1}(undef,nt)
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
        end
        plot!(1:nt,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=labl) #updates current plot
    end #of plotquantiles(y,col,labl,qnt)

    function plotqnt(x,y,col,labl,qnt::Float64=0.45)
        #variant on the above - with x-labels given and aligned
        nt,ny=size(y)
        low=Array{Float64,1}(undef,nt)
        mid=Array{Float64,1}(undef,nt)
        hig=Array{Float64,1}(undef,nt)
        for i in 1:nt
            low[i],mid[i],hig[i] = quantile(y[i,:],[0.5-qnt, 0.5, 0.5+qnt])
        end
        plot!(x,mid,grid=false,ribbon=(mid-low,hig-mid),fillalpha=.25,lw=3, seriescolor=col, label=labl)
    end

    function getdata(country::String="Australia", state::String="Western Australia", file::String="time_series_covid19_confirmed_global.csv")
    #load data from the named CSV file for the specified country and state
    #works for the non-US data files
        x=[]
        for row in CSV.File(file) #read row
            if row[2]==country
                if length(state)==0 && row[1]===missing #if country is spefiied and state is missing/absent
                    x=row
                elseif row[1]==state #otherwise, if country and state both match
                    x=row
                end
            end
        end

        if length(x)==0   #no such country/state found
            println("Error: "*country*", "*state*" not found.")
            return
        end

        #otherwise, load the data
        z=Array{Int64,1}(undef,length(x)-4)
        for (index, value) in enumerate(x)
            if index>4
                z[index-4]=value
            end
        end
        return z
    end #of getdata


    function getdataus(locality::String="Santa Fe", state::String="New Mexico", file::String="time_series_covid19_confirmed_US.csv")
    #load data from the named CSV file for the specified country and state
    #US-only
        x=[]
        for row in CSV.File(file)
            if row[7]==state
                if length(locality)==0 && row[6]===missing #if country is spefiied and state is missing/absent
                    x=row
                elseif row[6]==locality #otherwise, if country and state both match
                    x=row
                end
            end
        end

        if length(x)==0     #no such country/state found
            println("Error: "*locality*", "*state*" not found.")
            return
        end

        #otherwise, load the data
        z=Array{Int64,1}(undef,length(x)-11)
        for (index, value) in enumerate(x)
            if index>11
                z[index-11]=value
            end
        end

        return z
    end #of getdataus


    function epipred(th,St,Et,It,Rt)
        #find the timestep of the threshold th
        (lenp,nsims)=size(St)
        tstep=Array{Int,1}(undef,nsims)
        nbad=0
        for i in 1:nsims
            j=1
            while It[j,i]+Rt[j,i]<th && j<lenp
                j += 1
            end
            if j>=lenp
                j = -1
                nbad += 1
            elseif j==1
                j = -1
                nbad += 1
            end
            tstep[i] = j
        end
        ndays=lenp-maximum(tstep)
        Sp=Array{UInt64,2}(undef,ndays,nsims-nbad)
        Ep=Array{UInt64,2}(undef,ndays,nsims-nbad)
        Ip=Array{UInt64,2}(undef,ndays,nsims-nbad)
        Rp=Array{UInt64,2}(undef,ndays,nsims-nbad)
        i=1
        ii=1
        while ii<=nsims
            if tstep[ii] > 0
                Sp[:,i] = St[tstep[ii]:(tstep[ii]+ndays-1) , ii]
                Ep[:,i] = Et[tstep[ii]:(tstep[ii]+ndays-1) , ii]
                Ip[:,i] = It[tstep[ii]:(tstep[ii]+ndays-1) , ii]
                Rp[:,i] = Rt[tstep[ii]:(tstep[ii]+ndays-1) , ii]
                i += 1
            end
            ii += 1
        end
        return (Sp,Ep,Ip,Rp)
    end


###############################################################################
###############################################################################
###############################################################################
# This code is a massive kludge - essentially, what follows is when you have a
# model with a single changepoint a design a hack to introduce a second
###############################################################################
###############################################################################
    function episim3(net1, net2, net3, epiparam, feth::Int64=100, predays::Int64=40, middays::Int64=30, postdays::Int64=80, nsims::Int64=50)
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
            alldone=false
            donecounter=0
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

                #########################################
                #switch again after middays infections
                if ~notdoneyet && ~alldone
                    donecounter += 1
                    if donecounter>middays
                        net=net3
                        alldone=true
                    end
                end

                i += 1
            end
        end

        return st, ex, fe, rm

    end #of episim3
###############################################################################
###############################################################################
###############################################################################

end # of module
