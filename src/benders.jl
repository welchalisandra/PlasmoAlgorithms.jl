using JuMP
using Gurobi
using Plasmo
using DataFrames

function fix(var,value)
  setlowerbound(var,value)
  setupperbound(var,value)
end

##Place MP and SP into PlasmoGraph
mp = Model(solver = GurobiSolver())
sp = Model(solver = GurobiSolver())

@variable(mp,y>=0)
@variable(mp,θ>=0)
@objective(mp,Min,2y+θ)

@variable(sp,x[1:2]>=0)
@variable(sp,y>=0)
@constraint(sp,2x[1]-x[2]+3y>=4)
@constraint(sp,x[1]+2x[2]+y>=3)
@objective(sp,Min,2x[1]+3x[2])

## Plasmo Graph
g = PlasmoGraph()
g.solver = GurobiSolver()
n1 = add_node(g)
setmodel(n1,mp)
n2 = add_node(g)
setmodel(n2,sp)

## Linking constraints between MP and SP
@linkconstraint(g, n1[:y] == n2[:y])

function benderssolve(graph::PlasmoGraph;max_iterations=2)
# bendersolve(g, max_iterations=10)

  ########## 1. Initialize ########
  tic()
  starttime = time()
  # Results outputs
  #df = DataFrame(Iter=[],Time=[],α=[],step=[],UB=[],LB=[],Hk=[],Zk=[],Gap=[])
  res = Dict()

  #TODO changed in new version of Plasmo?
  mp = getmodel(graph.nodes[1])
  sp = getmodel(graph.nodes[2])

  solve(mp)

  mflat = create_flat_graph_model(graph)
  mflat.solver = graph.solver
  #print(mflat.linconstrDuals[1])
  print(mflat)
  solve(mflat,relaxation=true)
  links = getlinkconstraints(graph)
  nmult = length(links)

  λ = mflat.linconstrDuals[end-nmult+1:end]
  status = solve(sp)
  if status != :Optimal
    @constraint(mp,fc, 0>=λ*(getupperbound(y)-y))
    println(fc)
  else
    θk = getobjectivevalue(sp)
    @constraint(mp,test,θ >= θk + λ[1]*(getupperbound(y)-links[1].terms.vars[1]))
    println("theta",θk)
    println(test)
  end
  #solve(mp)
  #push!(df,[iter,round(time()-starttime),α,step,UB,LB,Hk,Zk,gap])
  return res
end

benderssolve(g)
