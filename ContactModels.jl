#module ContactModels

#	export fullmixing, limitpop, limitmix, isolation, socialdist, nomassmix, covidsafe, rewire!
	using LightGraphs, ProgressBars

	function fullmixing(pop)
	    #model a population of pop fully mixed nodes
	    #a BA scale-free network model, mean degree 4
	    if pop<5
	    	net=complete_graph(pop)
	    	return net
	    end
	    net=barabasi_albert(pop, 3, 2)
	    #need to add six extra edges
	    verts=vertices(net)
	    for i in 1:6
	        edg=Edge(rand(verts),rand(verts))
	        while edg in edges(net)
	            edg=Edge(rand(verts),rand(verts))
	        end
	        add_edge!(net,edg)
	    end
	    return net
	end

	function limitpop(net,gth)
	    #prune contact patterns to limit fatherings on net to 
	    #no more than gth people. edges are redistributed to preserve mean degree
	    nlost=0
	    nedges=ne(net)
	    #prune edges from hubs
	    for vert in vertices(net)
	        neigh=neighbors(net,vert)
	        while length(neigh)>gth
	            rem_edge!(net,vert,rand(neigh))
	            nlost +=1
	            neigh=neighbors(net,vert) 
	        end
	    end
	    #redistribute them elsewhere
	    vert=vertices(net)
	    for i in 1:nlost #will complain if there are not enough edges to add back in
	        verts=vert#[degree(net).<(gth-1)] 
	        v1=rand(verts)
	        v2=v1
	        while v2==v1
	            v2=rand(verts)
	        end
	        edg=Edge(v1,v2)
	        add_edge!(net,edg)
	    end
	    #maybe a few short, add those too
#	    while ne(net)<nedges
#	        edg=Edge(rand(verts),rand(verts))
#	        while edg in edges(net)
#	            edg=Edge(rand(verts),rand(verts))
#	        end
#	        add_edge!(net,edg)
#nomass	    end
	    return net
	end

	function limitmix(pop,gth)
	    #model a population of pop fully mixed nodes
	    #a BA scale-free network model, mean degree 4
	    #No more than gth people
	    return limitpop(fullmixing(pop),gth)
	end

	function isolation(pop)
	    #model a population of pop (must be a square number) nodes
	    #on a 2D grid, mean degree 4
	    spop=Int(floor(sqrt(pop)))
	    if abs(spop^2-pop)>0
	        spop1=Int(floor(sqrt(pop)))
	        spop2=Int(ceil(sqrt(pop)))
	        println("population not square ",pop-spop1*spop2," node(s) isolated")
	        net=LightGraphs.grid([spop1, spop2], periodic=true)
	        add_vertices!(net,pop-spop1*spop2)
	    else
	        net=LightGraphs.grid([spop, spop], periodic=true)
	    end
	    return net
	end

	function socialdist(pop,cmpl)
	    #model a population of pop (must be a square number) nodes
	    #on a 2D grid, with social distancing compliance cmpl, mean degree 4
	    q=1-cmpl^(1/4)  #as elsewhere, this whole thing is assume a mean degree of 4
	    net=isolation(pop)
	    rewire!(net,q)
	    return net
	end

	function nomassmix(pop)
	    #model a population of pop fully mixed nodes - no superspreaders
	    #a ER random graph model, mean degree 4
	    if pop<5
	    	return complete_graph(pop)
	    end
	    return erdos_renyi(pop, 2*pop)
	end

	function covidsafe(net,cvd)
	    #"covidsafe" at fraction cvd effectiveness
	    iso=cvd^2 #covidsafe implies removing iso fraction of all edges of graph
	    for edg in edges(net)
	        if rand(Float64) .< iso
	            rem_edge!(net,edg)
	        end
	    end
	    return net
	    #unlike the preceeding constructions, this one (and only this one) will reduce the overall degree
	end

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

	function gc_dist(pos1,pos2)
		π180=π/180
		ϕ1=pos1[2]*π180
		ϕ2=pos2[2]*π180
		θ1=pos1[1]*π180
		θ2=pos2[1]*π180
		d=( sin((ϕ2-ϕ1)/2) )^2 + cos(ϕ1)*cos(ϕ2)*( sin((θ2-θ1)/2) )^2
		d=asin(sqrt(d))
		d=2*6378*d
		d=d+1
		return d
	end

	function similaritytown(df,i,j,posn)
	    # similaritytown(df,i,j)
	    # Define a "similarity" between two locations. This code is an ad-hoc rule that is an 
	    # attempt to model the amount of movement between discrete regions. Dataframe df is 
	    # assumed to be in the "Orlando format" as illustrated above.
	    bonus1=1000  #effectively, the larger the bonus, the more we include non-indigeneos travel.
	    bonus2=10   #bonus1 is for towns in the same region, bonus2 for the same LGA
	    if i<j
	        score=similaritytown(df,j,i,posn)
	        return score
	    elseif i<=majorlocalities && j<=majorlocalities
	        x1= convert(Matrix{Int}, df[2*i-1:2*i,6:end])
	        x2= convert(Matrix{Int}, df[2*j-1:2*j,6:end])
	        languagesimilarity=x1[:]'*x2[:]        #or not
	        regionbonus= bonus1/gc_dist(posn[i,:],posn[j,:])
	    elseif i>majorlocalities && j>majorlocalities
	        x1= convert(Vector{Int}, df[majorlocalities+i,6:end])
	        x2= convert(Vector{Int}, df[majorlocalities+j,6:end])
	        languagesimilarity=x1[:]'*x2[:]        #or not
	        regionbonus= bonus1/gc_dist(posn[i,:],posn[j,:]) #same region
	        regionbonus= regionbonus + (df[i+majorlocalities,5]==df[j+majorlocalities,5])*bonus2 #and same locality
	    else
	        x1= convert(Vector{Int}, df[majorlocalities+i,6:end])
	        x1= reshape(x1,1,length(x1))
	        x2= convert(Matrix{Int}, df[2*j-1:2*j,6:end])
	        languagesimilarity= (sum(x2,dims=1) * x1[:])[1]      #or not
	        regionbonus= bonus1/gc_dist(posn[i,:],posn[j,:])#same region
	    end
	    score=Int(floor(sqrt(languagesimilarity+regionbonus))) # sqrt? why not?
	    return score
	end

	function similaritytown(df,i,j)
	    # similaritytown(df,i,j)
	    # Define a "similarity" between two locations. This code is an ad-hoc rule that is an 
	    # attempt to model the amount of movement between discrete regions. Dataframe df is 
	    # assumed to be in the "Orlando format" as illustrated above.
	    bonus1=100  #effectively, the larger the bonus, the more we include non-indigeneos travel.
	    bonus2=10   #bonus1 is for towns in the same region, bonus2 for the same LGA
	    if i<j
	        score=similaritytown(df,j,i)
	        return score
	    elseif i<=majorlocalities && j<=majorlocalities
	        x1= convert(Matrix{Int}, df[2*i-1:2*i,6:end])
	        x2= convert(Matrix{Int}, df[2*j-1:2*j,6:end])
	        languagesimilarity=x1[:]'*x2[:]        #or not
	        regionbonus= (df[2*i,3]==df[2*j,3])*bonus1
	    elseif i>majorlocalities && j>majorlocalities
	        x1= convert(Vector{Int}, df[majorlocalities+i,6:end])
	        x2= convert(Vector{Int}, df[majorlocalities+j,6:end])
	        languagesimilarity=x1[:]'*x2[:]        #or not
	        regionbonus= (df[i+majorlocalities,3]==df[j+majorlocalities,3])*bonus1 #same region
	        regionbonus= regionbonus + (df[i+majorlocalities,5]==df[j+majorlocalities,5])*bonus2 #and same locality
	    else
	        x1= convert(Vector{Int}, df[majorlocalities+i,6:end])
	        x1= reshape(x1,1,length(x1))
	        x2= convert(Matrix{Int}, df[2*j-1:2*j,6:end])
	        languagesimilarity= (sum(x2,dims=1) * x1[:])[1]      #or not
	        regionbonus= (df[i+majorlocalities,3]==df[2*j,3])*bonus1 #same region
	    end
	    score=Int(floor(sqrt(languagesimilarity+regionbonus))) # sqrt? why not?
	    return score
	end


	function vulnerability(It)
	    (ndays, nsims, nlocs) = size(It)
	    vln=Array{Float64,2}(undef,nlocs,nsims)
	    for i in 1:nlocs
	        for j in 1:nsims
	            idx=findall(x -> x>0, It[:,j,i])
	            if isempty(idx)
	                vln[i,j]=Inf
	            else
	                vln[i,j]=minimum(idx)
	            end
	        end
	    end
	    vln=median(vln,dims=2)
	    vln=vln./vln[3] #median vulnerability comapred to Gero
	    return vln
	end


	function buildstate(statedata,netbuilder, p, smalltown::Bool=true, richclub::Bool=false, posn::Any=false)
	    # buildstate(statedata, netbuilder, p)
	    ##############################################################################
	    #statedata is a csv table loaded from the structure /format above (AKA the Orlando Format)
	    #netbuilder is the intra-locale network construction rule (a function)
	    #smalltown==true will bias toward travel to/from remote communities, otherwise travel will be predominantly to the largest places
	    #richclub==true preferences "popular" (hihgly connected) nodes as travellers
	    ##############################################################################
	    # 
	    #extract the data from CSV table
	    #majorlocalities=58
	    brk=majorlocalities*2
	    locale=[statedata[1:2:brk,1]; statedata[brk+1:end,1]]
	    townpop1=wapopint=parse.(Int, replace.(statedata[1:2:brk,2], r","=> ""))
	    townpop2=wapopint=parse.(Int, replace.(statedata[brk+1:end,2], r","=> ""))
	    popl=[townpop1; townpop2]
	    npl=length(popl)
	    net=SimpleGraph() #empty graph
	    transit=Array{Int64,2}(undef,npl,npl) #number of transits beween locale[i] and locale[j]
	    #################################
	    # Population connectivity within regions
	    #
	    for (i,town) in enumerate(locale)
	        println("Adding ",town," (population: ",popl[i],")")
	        #add the intralocale links
	        if popl[i]>0
	            netadd=netbuilder(popl[i]) #use the preassigned method to add the new component
	            net=blockdiag(net,netadd)
	        end
	        #compute the amount of transit
	        for j in 1:(i-1) #number of transit between here and every other previous part
	        	if typeof(posn)==Bool
	        		transit[i,j] = Int(floor(similaritytown(statedata,i,j)))
	            else
	            	transit[i,j] = Int(floor(similaritytown(statedata,i,j,posn)))
	            end
	            transit[j,i]=transit[i,j]
	        end
	        transit[i,i]=0
	    end
	    #################################
	    # Population connectivity between regions
	    #    
	    rr=[0; cumsum(popl)]
	#    tpopl=sum(sqrt.(popl))/length(popl)
	    tpopl=sum(popl)/sum(sqrt.(popl))
	    nedges_add = 0
	    println("Connecting towns")
	    iter = ProgressBar(1:npl) #everyone loves a good progress bar

        #main iteratarion loop nsims simulations
	    for i in iter
	        #do one of the following two lines
	        if smalltown
	            addlink=Int(floor(minimum([popl[i],tpopl*p*sqrt(popl[i])]))) #biased to small communities -  testing the effect of the hypothesis of more movement in these communities
	        else
	            addlink=Int(floor(p*popl[i])) #unbiased
	        end
	        edg1=rand(collect(rr[i]+1:rr[i+1]),addlink)
	        destpdf=transit[i,:] .* maximum([popl[i]./popl[:];1])
	        destpdf[i]=0
	        edg2 = sample(1:npl, Weights(destpdf) , addlink)
	        for k in 1:addlink
	            rrs=rr[edg2[k]]+1:rr[edg2[k]+1]
	            #do one of the following two lines
	            if richclub
	                edg2k = sample(rrs, Weights(degree(net[rrs]))) #biased by target degree
	            else
	                edg2k = rand(collect(rrs))   #random sample of the relevant community
	            end
	            edg2[k] = edg2k
	        end
	        #need to do it twice, so that the node degrees don't grow   
	        for k in 1:addlink
	                    add_edge!(net,edg1[k],edg2[k]) #this is grossly inefficient
	                    nedges_add += 1
	        end
	    end
	    #all done
	    println(" $nedges_add edges added")
	    return net, transit, locale, popl, nedges_add
	end
#end #of module

