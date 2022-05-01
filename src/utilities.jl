###DEPRECATED###WILL BE REMOVE SOON###

"""
    find_metabolite_by_formula(query::String, database)

Searching function for every metabolite that matches the chemical formula. Eases process of finding the desired metabolites faster in larger data arrays.
The `query` string should preferably be formatted according to the Hill notation. The `database` needs to be of one of the model types.
Please note that the pattern matching is exact and thus is prone to typing errors or alternate formulation (e.g. "NaCl" instead of "ClNa").
"""
function find_metabolite_by_formula(database, query::String)
    res_ls = String[]
    for met in metabolites(database)
        refstr = ""
        for (elem, elnum) in metabolite_formula(database, met)
            if elnum != 1
                refstr = string(refstr, elem, elnum)
            else
                refstr = string(refstr, elem)
            end
        end
        if refstr == query
            push!(res_ls, met)
        end
    end
    if length(res_ls) == 0
        println("No metabolite with matching formula found.")
    else
        return res_ls
    end
end

function print_formula(stoichiometry)
    leftstr = string()
    rightstr = string()
    for (k, v) in stoichiometry
        if v < 0
            leftstr = string(leftstr, string(abs(v)), k, " + ")
        elseif v > 0
            rightstr = string(rightstr, string(abs(v)), k, " + ")
        end
    end
    println(string(leftstr[1:end-3], " => ", rightstr[1:end-3]))
end
