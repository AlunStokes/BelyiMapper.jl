# Copyright (c) 2022 Alun Cennyth Stokes
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function generate_dessin(
    B::Vector{<:Integer},
    W::Vector{<:Integer},
    F::Vector{<:Integer},
    b_distinguished_val::Complex{<:Real}=0 + 0im,
    w_distinguished_val::Complex{<:Real}=1 + 0im,
    f_distinguished_val::Complex{<:Real}=Inf * im,
    dessin_name::AbstractString=get_name(B, W, F),
    base_folder::AbstractString="dessin-data",
    verbose::Bool=true,
    output_files::Bool=true,
    use_λ::Bool=true;
    index::Union{Integer,Nothing}=nothing
)

    if !isnothing(index)
        dessin_name *= "ind-$index"
    end

    @assert sum(B) == sum(W) == sum(F)
    deg = sum(B)

    dessin_type = "valency"

    b_solns = Vector{Tuple{<:Integer,Complex}}()
    w_solns = Vector{Tuple{<:Integer,Complex}}()
    f_solns = Vector{Tuple{<:Integer,Complex}}()

    push!(b_solns, (1, b_distinguished_val))
    push!(w_solns, (1, w_distinguished_val))
    push!(f_solns, (1, f_distinguished_val))

    use_F_ = false
    F_ = Nothing
    if isinf(f_distinguished_val)
        use_F_ = true
        F_ = F[2:end]
    end

    B_poly, b, z = get_real_poly(B, "b")
    W_poly, w, _ = get_real_poly(W, "w")
    F_poly = Nothing
    f = Nothing
    if use_F_
        F_poly, f, _ = get_real_poly(F_, "f", 1)
    else
        F_poly, f, _ = get_real_poly(F, "f")
    end

    B_poly = evaluate(B_poly, [b[1]] => [b_distinguished_val])
    W_poly = evaluate(W_poly, [w[1]] => [w_distinguished_val])
    if !use_F_
        F_poly = evaluate(F_poly, [f[1]] => [f_distinguished_val])
    end

    @var λ

    P = λ * B_poly - λ * W_poly - F_poly

    if use_λ
        println("B(z) = $B_poly")
        println("W(z) = $W_poly")
        println("F(z) = $F_poly\n")
        println("0 = $P")
        println("\n")
    else
        println("B(z) = $(replace(string(B_poly), "λ"=>"c"))")
        println("W(z) = $(replace(string(W_poly), "λ"=>"c"))")
        println("F(z) = $(replace(string(F_poly), "λ"=>"c"))")
        println("0 = $(replace(string(P), "λ"=>"c"))")
        println("\n")
    end

    eqns = generate_equations(P, z)
    #eqns = [horner(eqn) for eqn in eqns]

    Sys = System(eqns)
    Sys = optimize(Sys)
    #if verbose
    #    display(Sys)
    #end
    #vars = vcat(b, f, w, [λ])
    vars = variables(Sys)

    if verbose
        if use_λ
            println(join([string(v) for v in vars], " "))
        else
            @var c
            println(join([string(v) for v in vcat(vars[1:end-1], [c])], " "))
        end
        println("\n")
    end

    if verbose
        print("[")
        for (i, eqn) in enumerate(eqns)
            if !isequal(i, size(eqns, 1))
                if use_λ
                    print("$eqn,\n")
                else
                    print("$(replace(string(eqn), "λ"=>"c")),\n")
                end
            else
                if use_λ
                    print("$eqn]")
                else
                    print("$(replace(string(eqn), "λ"=>"c"))]")
                end
            end
        end
    end

    println("")

    result = solve(Sys; show_progress=verbose)

    if verbose
        println()
        println(result)
    end

    if nresults(result) == 0
        return
    end

    cert_result = Nothing
    try
        cert_result = certify(Sys, result)
    catch
        println("Failed to certify result.")
        return
    end

    if verbose
        display(cert_result)
    end

    certs = [cert for cert in certificates(cert_result) if is_certified(cert)]
    paths = [path for path in path_results(result) if is_success(path) && is_nonsingular(path)]
    solns = [solution_candidate(cert) for cert in certs]

    if verbose
        println(solns[1])
    end

    if output_files

        if output_files
            hcjl_output_folder = Nothing
            if isequal(dessin_type, "valency")
                hcjl_output_folder = "./$base_folder/$dessin_type/deg_$deg/$dessin_name/hcjl/"
            else
                hcjl_output_folder = "./$base_folder/$dessin_type/$dessin_name/hcjl/"
            end

            mkpath(hcjl_output_folder)
            println("\nSaving in $hcjl_output_folder\n")
        end

        for (n, (soln, path)) in enumerate(zip(solns, paths))
            #display(path)
            #println("\n----------------------\n")
            b_solns_res = copy(b_solns)
            w_solns_res = copy(w_solns)
            f_solns_res = copy(f_solns)
            for (var, val) in zip(variables(Sys), soln)
                var = string(var)
                if var[1:1] == "b"
                    num = parse(Int64, var[2:end])
                    push!(b_solns_res, (num, val))
                elseif var[1:1] == "w"
                    num = parse(Int64, var[2:end])
                    push!(w_solns_res, (num, val))
                elseif var[1:1] == "f"
                    num = parse(Int64, var[2:end])
                    push!(f_solns_res, (num, val))
                else
                end
            end

            b_output = ColouredOutputObject(B, b_solns_res)
            w_output = ColouredOutputObject(W, w_solns_res)
            f_output = ColouredOutputObject(F, f_solns_res)

            λ_val = soln[end]

            insert!(soln, 1, b_distinguished_val)
            insert!(soln, size(b)[1] + 2, f_distinguished_val)
            insert!(soln, size(b)[1] + size(f)[1] + 3, w_distinguished_val)

            output = OutputObject(n, vars, soln, b_output, w_output, f_output, λ_val, B_poly, W_poly, F_poly)

            output_solution(output, eqns, path, hcjl_output_folder, string(n), use_λ)

        end
    end

    return solns

end
