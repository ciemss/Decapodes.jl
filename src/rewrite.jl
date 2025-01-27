using Catlab
using AlgebraicRewriting
using AlgebraicRewriting: rewrite

#""" function get_valid_op1s(deca_source::SummationDecapode, varID)
#Searches SummationDecapode, deca_source, at the request varID
#and returns all op1s which are allowed to be averaged. Returns
#an array of indices of valid op1 sources.
#
#Namely this is meant to exclude ∂ₜ from being included in an average.
#"""
function get_valid_op1s(deca_source::SummationDecapode, varID)
    # skip_ops = Set([:∂ₜ])
    indices = incident(deca_source, varID, :tgt)
    return filter!(x -> deca_source[x, :op1] != :∂ₜ, indices)
end
  
#""" function is_tgt_of_many_ops(d::SummationDecapode, var)
#Return true if there are two or more distinct operations leading
#into Var var (not counting ∂ₜ).
#"""
function is_tgt_of_many_ops(d::SummationDecapode, var)
  op1Count = length(get_valid_op1s(d, var))
  op2Count = length(incident(d, var, :res))
  sumCount = length(incident(d, var, :sum))

  op1Count + op2Count + sumCount >= 2
end
  
#""" function find_tgts_of_many_ops(d::SummationDecapode)
#Searches SummationDecapode, d, for all Vars which have two or
#more distinct operations leading into the same variable.
#"""
function find_tgts_of_many_ops(d::SummationDecapode)
  filter(var -> is_tgt_of_many_ops(d, var), parts(d, :Var))
end
  
#""" function get_preprocess_indices(deca_source::SummationDecapode)
#Searches SummationDecapode, deca_source, for all Vars which are
#valid for average rewriting preprocessing. Namely this just includes
#all op2 and summation operations. Returns two arrays, first is 
#array of valid Op2 ids, second is array of valid Σ ids.
#"""
function get_preprocess_indices(deca_source::SummationDecapode)
    targetOp2 = []
    targetSum = []

    targetVars = find_tgts_of_many_ops(deca_source)

    for var in targetVars
        append!(targetOp2, incident(deca_source, var, :res))
        append!(targetSum, incident(deca_source, var, :sum))
    end

    return targetOp2, targetSum
end
  
#""" function preprocess_average_rewrite(deca_source::SummationDecapode)
#Preprocesses SummationDecapode, deca_source, for easier average 
#rewriting later on. Specifically, all op2 and summation results are
#stored in variables called "Temp" and results are then passed off to 
#their original result along an op1 called "temp". This "temp" operation
#is equivalent to an identity function.
#"""
function preprocess_average_rewrite(deca_source::SummationDecapode)
    targetOp2, targetSum = get_preprocess_indices(deca_source)

    # If we don't need to preprocess then don't
    if(length(targetOp2) == 0 && length(targetSum) == 0)
        return deca_source
    end

    LHS = []
    RHS = []

    SuperMatch = []
    SuperVarMap = Vector{Int}()
    SuperOp2Map = Vector{Int}()
    SuperSigmaMap = Vector{Int}()
    SuperSummandMap = Vector{Int}()

    serial = 0
    # Process all of the target rewrites for op2
    for opID in targetOp2

        vars = [deca_source[opID, :proj1], deca_source[opID, :proj2], deca_source[opID, :res]]
        types = deca_source[vars, :type]
        names = deca_source[vars, :name]
        op2name = deca_source[opID, :op2]

        Match = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = 3
            type = types
            name = names

            Op2 = 1
            proj1 = [1]
            proj2 = [2]
            res = [3]
            op2 = op2name
        end

        I = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = 3
            type = types
            name = names
        end

        Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = 4
            type = vcat(types, types[end])
            name = vcat(names, Symbol("Temp_", serial))

            Op1 = 1
            src = [4]
            tgt = [3]
            op1 = [:temp]

            Op2 = 1
            proj1 = [1]
            proj2 = [2]
            res = [4]
            op2 = op2name
        end

        serial += 1

        L = ACSetTransformation(I, Match, Var = 1:3)
        R = ACSetTransformation(I, Sub, Var = 1:3)

        push!(LHS, L)
        push!(RHS, R)
        push!(SuperMatch, Match)

        append!(SuperVarMap, vars)
        push!(SuperOp2Map, opID)
    end
    
    # Process all of the target rewrites for sums
    for sumID in targetSum
        summandIDs = incident(deca_source, sumID, :summation)
        vars = vcat(deca_source[summandIDs, :summand], deca_source[sumID, :sum])
        types = deca_source[vars, :type]
        names = deca_source[vars, :name]
        
        rewrite_size = length(vars)
        nary = rewrite_size - 1

        Match = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = rewrite_size
            type = types
            name = names

            Σ = 1
            sum = [rewrite_size]

            Summand = nary
            summand = 1:nary
            summation = fill(1, nary)
        end

        I = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = rewrite_size
            type = types
            name = names
        end

        Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = rewrite_size + 1
            type = vcat(types, types[end])
            name = vcat(names, Symbol("Temp_", serial))

            Op1 = 1
            src = [rewrite_size + 1]
            tgt = [rewrite_size]
            op1 = [:temp]

            Σ = 1
            sum = [rewrite_size + 1]

            Summand = nary
            summand = 1:nary
            summation = fill(1, nary)
        end

        serial += 1

        L = ACSetTransformation(I, Match, Var = 1:rewrite_size)
        R = ACSetTransformation(I, Sub, Var = 1:rewrite_size)

        push!(LHS, L)
        push!(RHS, R)
        push!(SuperMatch, Match)

        append!(SuperVarMap, vars)
        push!(SuperSigmaMap, sumID)
        append!(SuperSummandMap, summandIDs)
    end

    # Combine all rules in parallel and apply
    rule = Rule(oplus(LHS), oplus(RHS))
    m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op2 = SuperOp2Map, Σ = SuperSigmaMap, Summand = SuperSummandMap)

    rewrite_match(rule, m)
end
  
#""" function process_average_rewrite(deca_source::SummationDecapode)
#Rewrites SummationDecapode, deca_source, by including averages
#of redundent operations. While this function only searches for op1s
#to match on, because of preprocessing, this indirectly includes op2 
#and summations in the final result.
#"""
function process_average_rewrite(deca_source::SummationDecapode)
    targetVars = find_tgts_of_many_ops(deca_source)

    # If no rewrites available, then don't rewrite
    if(length(targetVars) == 0)
        return deca_source
    end

    LHS = []
    RHS = []

    SuperMatch = []
    SuperVarMap = Vector{Int}()
    SuperOp1Map = Vector{Int}()

    varSerial = 0
    sumSerial = 0

    # Create rewrite rules for multiple op1s leading into target
    for varID in targetVars
        targetOp1 = get_valid_op1s(deca_source, varID)
        vars = vcat(deca_source[targetOp1, :src], varID)

        num_nodes_match = length(vars)
        nary_of_rewrite = num_nodes_match - 1

        result_index = num_nodes_match
        sum_index = 2 * result_index

        variable_types = deca_source[vars, :type]
        variable_var =  deca_source[vars, :name]
        variable_op1 = deca_source[targetOp1, :op1]

        Match = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = num_nodes_match
            type = variable_types
            name = variable_var

            Op1 = nary_of_rewrite
            src = 1:nary_of_rewrite
            tgt = fill(result_index, nary_of_rewrite)
            op1 = variable_op1
        end

        I = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = num_nodes_match
            type = variable_types
            name = variable_var
        end 

        Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
            Var = 2 * num_nodes_match
            type = vcat(variable_types, variable_types)
            name = vcat(variable_var, map(x -> Symbol("••", varSerial + x), 1:nary_of_rewrite), [Symbol("••sum", sumSerial)])
            Op1 = nary_of_rewrite + 1
            src = vcat(1:nary_of_rewrite, sum_index)
            tgt = vcat(num_nodes_match+1:sum_index-1, [result_index])
            op1 = vcat(variable_op1, Symbol(:avg, nary_of_rewrite))
            Σ = 1
            sum = [sum_index]
            Summand = nary_of_rewrite
            summand = num_nodes_match+1:sum_index-1
            summation = fill(1, nary_of_rewrite)
        end

        varSerial += nary_of_rewrite
        sumSerial += 1

        L = ACSetTransformation(I, Match, Var = 1:num_nodes_match);
        R = ACSetTransformation(I, Sub, Var = 1:num_nodes_match);
        
        push!(LHS, L)
        push!(RHS, R)
        push!(SuperMatch, Match)

        append!(SuperVarMap, vars)
        append!(SuperOp1Map, targetOp1)
    end

    # Combine all rewrite rules and apply at once
    rule = Rule(oplus(LHS), oplus(RHS))

    m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op1 = SuperOp1Map)
    rewrite_match(rule, m)
end
  
"""    function average_rewrite(deca_source::SummationDecapode)

Compute each quantitity in the given Decapode by the average of all computation paths leading to that node.
"""
function average_rewrite(deca_source::SummationDecapode)
    return process_average_rewrite(preprocess_average_rewrite(deca_source))
end

#""" function find_variable_mapping(deca_source, deca_tgt)
#Returns array of match on variables between from a decapode
#source to a target, also returns if mapping is valid
#WARNING: This assumes that variable names are unique.
#If variable names are not unique or do not exist, 
#corrsponding mapping value is set to 0.
#"""
function find_variable_mapping(deca_source, deca_tgt)

    mapping = []
    valid_matching = true
  
    for varID in parts(deca_source, :Var)
      varName = deca_source[varID, :name]
      matches = incident(deca_tgt, varName, :name)
      if(length(matches) >= 2)
        # println("Name for variable at index $varID named $varName is not unique. Setting to 0.")
        valid_matching = false
        push!(mapping, 0)
      elseif(length(matches) == 0)
        # println("Variable at index $varID named $varName does not exist in target. Setting to 0.")
        valid_matching = false
        push!(mapping, 0)
      else
        push!(mapping, only(matches))
      end
    end
  
    return mapping, valid_matching
end
