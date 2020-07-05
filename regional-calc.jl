using LightGraphs, Plots, GraphPlot, Compose, Random, StatsBase, CSV, JLD2, FileIO, Statistics


include("ContactModels.jl")
include("EpiSim.jl")

wapop=CSV.read("all-wa.csv")
waloc=CSV.read("all-wa-pos.csv")


#########################################
#the following parameter is hardcoded here, even though it is used within ContactModels.jl just to ensure it is not forgotten....
global majorlocalities=58    #global scope intended. 
global brk=majorlocalities*2
#########################################


#### remove empty towns
(nr,nc)=size(wapop)
townpop1=parse.(Int, replace.(wapop[1:2:brk,2], r","=> ""))
townpop2=parse.(Int, replace.(wapop[brk+1:end,2], r","=> ""))
popl=[townpop1; townpop2]
rowidx=Int64[]
for i in majorlocalities+1:length(popl)
    if popl[i]!=0
        append!(rowidx, i)
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

########
#'reasonable' parameters ################################
epiparam=Dict()
epiparam["p0"]=0.2 #a guess - tuned to match observed data 
epiparam["q"]=1/7 #"up to" two weeks
epiparam["r0"]=1/14 #about two weeks for mild, 3-6 for severe
epiparam["nseeds"]=5 #probably too many, consider dropping. 
epiparam["pop"]=sum(popl) # =size(net)[1]
#########################################################

contact=[x -> covidsafe(nomassmix(x),0.5),
 x -> covidsafe(nomassmix(x),0.3),
 x -> nomassmix(x),
 x -> limitmix(x,50),
 x -> limitmix(x,500),
 x -> fullmixing(x)
]

labels=["covidsafe50", "covidsafe30", "nomassmix", "limitmix50", "limitmix500", "fullmixing" ]

function b2s(var::Bool)
	if var
		return "T"
	else
		return "F"
	end
end



######## build contact graph 
pt=0.01 #probability of travel/proportion of population connecting multiple regions
ndays=210
nsims=250
locpop=cumsum([1; popl],dims=1)

vresult=Array{Float64,2}(undef,length(contact)*4,length(popl))

i=1
for (ci,cfunc) in enumerate(contact)
	for smallc in [true false]
		for prefpa in [true false]
			global i
			filename="regional_"*labels[ci]*"_"*b2s(smallc)*b2s(prefpa)*"_p"*string(Int(0.01*100); base=10,pad=2)*".jld"
			net, transit, locale, popl, nedges_added = 
			                 buildstate(wapop, cfunc, pt, smallc, prefpa, posn)
			St,Et,It,Rt = 
					EpiSim.episimcom(net, locpop, epiparam, ndays, nsims) 
			vv=vulnerability(It)
			vresult[i,:]=vv'
		#	save(filename,"St",St,"Et",Et,"It",It,"Rt",Rt) 
			@save "regional_progress_3.jld"
			save("regional_output_3.jld","vuln",vresult)
			println("Done $i of $(length(contact)*4)")
			i += 1
		end
	end
end

exit()





