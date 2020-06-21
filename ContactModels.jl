#module ContactModels

#	export fullmixing, limitpop, limitmix, isolation, socialdist, nomassmix, covidsafe, rewire!
	using LightGraphs

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

#end #of module

