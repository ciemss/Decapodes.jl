### Drawing of Decapodes
import Catlab.Graphics.Graphviz
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs
using Catlab.Graphs
using Catlab.Graphs.BasicGraphs

#= reg_to_sub = Dict('0'=>'₀', '1'=>"₁", '2'=>'₂', '3'=>'₃', '4'=>'₄',
    '5'=>'₅', '6'=>'₆','7'=>'₇', '8'=>'₈', '9'=>'₉', 'r'=>'•', 'l'=>'L', ) =#

reg_to_sub = Dict("Form0"=>'₀', "Form1"=>"₁", "Form2"=>'₂', "Form3"=>'₃', 
                  "DualForm0"=>'₀', "DualForm1"=>"₁", "DualForm2"=>'₂', "DualForm3"=>'₃', 
                  "infer"=>'•', "Literal"=>'L', "Constant"=>'C', "Parameter"=>'P', )

toSub(digit::String) = get(reg_to_sub, digit, digit)

spacename(d, v) = begin
    t = d[v, :type]
    subscript = toSub(String(t))
    dom = startswith(String(t), "Dual") ? "Ω̃" : "Ω"
    # TODO: To get rid of omega for non-form types
    # dom = endswith(String(t), r"[0-9]") ? dom : ""
    return "$dom$subscript"
end
varname(d, v) = "$(d[v, :name]):$(spacename(d, v))"

# TODO: Change orientation to print 
"""    Graphics.to_graphviz(F::AbstractDecapode; directed = true, kw...)

Visualize the given Decapode through Graphviz. Ensure that you have called `using Catlab.Graphics` before-hand, and have a way of visualizing SVG files in your current environment.
"""
Graphics.to_graphviz(F::AbstractDecapode; directed = true, kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(F; directed, kw...))

decapode_edge_label(s::Symbol) = String(s)
decapode_edge_label(s::Vector{Symbol}) = join(String.(s), "⋅")
decapode_edge_label(s::String) = s
decapode_edge_label(s::Vector{String}) = join(s, "⋅")


function Catlab.Graphics.to_graphviz_property_graph(d::AbstractNamedDecapode, directed = true; kw...)
    pg = PropertyGraph{Any}(;kw...)
    vids = map(parts(d, :Var)) do v
      add_vertex!(pg, label=varname(d,v))
    end

    # Add entry and exit vertices and wires
    if(directed)
      for state in infer_states(d)
        tempv = add_vertex!(pg)
        add_edge!(pg, tempv, vids[state])
      end
    end

    map(parts(d, :Op1)) do op
      s, t = d[op, :src], d[op, :tgt]
      add_edge!(pg, vids[s],vids[t], label=decapode_edge_label(d[op,:op1]))
    end

    map(parts(d, :Op2)) do op
      s, t = d[op, :proj1], d[op, :proj2]
      r = d[op, :res]
      v = add_vertex!(pg, label="$(spacename(d, s))×$(spacename(d,t))", shape="rectangle")

      # If in directed mode, sources point into the projection field
      # Else, everything points out
      if(directed)
        add_edge!(pg, vids[s], v, label="π₁", style="dashed")
        add_edge!(pg, vids[t], v, label="π₂", style="dashed")
      else
        add_edge!(pg, v, vids[s], label="π₁", style="dashed")
        add_edge!(pg, v, vids[t], label="π₂", style="dashed")
      end
      add_edge!(pg, v, vids[r], label=decapode_edge_label(d[op, :op2]))
    end

    return pg
end

function Catlab.Graphics.to_graphviz_property_graph(d::SummationDecapode; directed = true, prog = "dot", node_attrs=Dict(), edge_attrs=Dict(), graph_attrs=Dict(), node_labels = true, kw...)
    
    default_graph_attrs = Dict(:rankdir => "TB")
    default_edge_attrs = Dict()
    default_node_attrs = Dict(:shape => "oval")

    G = to_graphviz_property_graph(Catlab.Graphs.Graph(0); prog, node_labels,
      node_attrs = merge!(default_node_attrs, node_attrs),  
      edge_attrs = merge!(default_edge_attrs, edge_attrs),  
      graph_attrs = merge!(default_graph_attrs, graph_attrs))

    vids = map(parts(d, :Var)) do v
      add_vertex!(G, label=varname(d,v))
    end

    # Add entry and exit vertices and wires
    if(directed)
      tempin = add_vertex!(G, shape = "none", label = "")
      for state in infer_states(d)
        add_edge!(G, tempin, vids[state])
      end

      tempout = add_vertex!(G, shape = "none", label = "")
      for tvar in d[:incl]
        add_edge!(G, vids[tvar], tempout)
      end
    end

    map(parts(d, :Op1)) do op
      s, t = d[op, :src], d[op, :tgt]
      add_edge!(G, vids[s],vids[t], label=decapode_edge_label(d[op,:op1]))
    end

    map(parts(d, :Op2)) do op
      s, t = d[op, :proj1], d[op, :proj2]
      r = d[op, :res]
      v = add_vertex!(G, label="$(spacename(d, s))×$(spacename(d,t))", shape="rectangle")

      # If in directed mode, sources point into the projection field
      # Else, everything points out
      if(directed)
        add_edge!(G, vids[s], v, label="π₁", style="dashed")
        add_edge!(G, vids[t], v, label="π₂", style="dashed")
      else
        add_edge!(G, v, vids[s], label="π₁", style="dashed")
        add_edge!(G, v, vids[t], label="π₂", style="dashed")
      end
      add_edge!(G, v, vids[r], label=decapode_edge_label(d[op, :op2]))
    end

    findvid(G, d, v) = incident(G.graph, [Dict{Symbol, Any}(:label=>varname(d, v))], :vprops)
    white_nodes = map(parts(d, :Σ)) do s
        v = add_vertex!(G, label="Σ$s", shape="circle")
        u = d[s, :sum]
        matches = first(findvid(G, d, u))
        length(matches) == 1 || error("did not find a unique vertex match for Σ$s")
        uG = first(matches)
        add_edge!(G, v, uG, label="+")
        return v
    end
    for e in parts(d, :Summand)
        # If directed, point summands into the sum
        # Else, everything points outward
        if(directed)
          e = add_edge!(G, d[e, :summand], white_nodes[d[e, :summation]], style="dashed")
        else
          e = add_edge!(G, white_nodes[d[e, :summation]], d[e, :summand], style="dashed")
        end
    end
    return G
end

savevizsvg(g, fname::String) = open(fname, "w") do fp
  run_graphviz(fp, to_graphviz(to_graphviz_property_graph(nsdp)), prog="neato", format="svg")
end
