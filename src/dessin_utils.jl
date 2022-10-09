# Copyright (c) 2022 Alun Cennyth Stokes
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

struct ColouredOutputObject
    valencies::Vector{Int64}
    vals::Vector{Tuple{Int64,Complex}}
end

struct OutputObject
    solution_num::Int
    vars::Vector{Variable}
    vals::Vector{Complex}
    b_output::ColouredOutputObject
    w_output::ColouredOutputObject
    f_output::ColouredOutputObject
    λ_val::Complex
    B_poly::Expression
    W_poly::Expression
    F_poly::Expression
end

function get_monic_poly(degree::Integer, power::Integer, var::AbstractString, var_offset::Integer)
    @var z
    v = [Variable("$(var)$i") for i = (var_offset+1):(var_offset+degree)]

    poly = z^degree
    for (j, (i)) in enumerate(1:1:(degree))
        poly += v[i] * z^(degree - j)
    end
    poly = poly^power
    return poly
end

#Groups m repeated v-valencies into a single, monic, degree m poly., raised to exponent v
function get_real_poly(valencies::Vector{<:Integer}, var::AbstractString, var_offset::Integer=0)
    @var z
    v = [Variable("$(var)$i") for i = 1:(size(valencies)[1])]

    valencies = copy(valencies)

    poly = 1
    var_counter = 1 + var_offset

    valency_counts = Dict([(i, count(x -> x == i, valencies)) for i in unique(valencies)])

    for i in unique(valencies)
        deg = i
        count = valency_counts[i]
        poly *= get_monic_poly(count, deg, var, var_counter - 1)
        var_counter += count
    end

    return poly, v, z
end

function generate_equations(P::Expression, var::Variable)
    _, p_coefs = exponents_coefficients(P, [var])

    return p_coefs
end

function format_equation_string(eqn::Union{AbstractString,Expression}, use_λ::Bool=true)
    eqn_string = string(eqn)
    eqn_string = replace(eqn_string, "*im" => "j")
    eqn_string = replace(eqn_string, "^" => "**")
    if !use_λ
        eqn_string = replace(eqn_string, "λ" => "c")
    end
    return eqn_string
end


function leading_zeros(s::AbstractString)
    i = 1
    while i <= length(s)
        if !isequal(s[i:i], "0")
            break
        end
        i += 1
    end
    return i - 1
end


function output_solution(output::OutputObject, eqns::Vector{Expression}, path::PathResult, folder_path::AbstractString, file::AbstractString, use_λ::Bool=true)

    b_valencies = copy(output.b_output.valencies)
    w_valencies = copy(output.w_output.valencies)
    f_valencies = copy(output.f_output.valencies)

    mkpath(folder_path * "solns/")


    residual_string = string(residual(path))

    if occursin("e", residual_string)
        exponent_string = split(split(residual_string, ".")[end], "e")[end][2:end]
        b10_precision = parse(Int64, exponent_string) - 1
    else
        exponent_string = split(residual_string, ".")[end]
        b10_precision = leading_zeros(exponent_string)
    end

    value_precision_string = "{:1.$(b10_precision)e}"
    value_precision_fe = FormatExpr(value_precision_string)


    info = Dict(
        "solution_num" => output.solution_num,
        "unknowns_are_roots" => false,
        "includes_min_polys" => false,
        "uses_lambd" => use_λ,
        "b10_precision" => b10_precision,
    )

    valencies = Dict(
        "B" => b_valencies,
        "W" => w_valencies,
        "F" => f_valencies
    )

    equations = [format_equation_string(eqn, use_λ) for eqn in eqns]

    solutions = Dict{String,Any}()
    for (n, v) in sort(output.b_output.vals, by=first)
        solutions["b$n"] = [format(value_precision_fe, real(v)), format(value_precision_fe, imag(v))]
    end
    for (n, v) in sort(output.w_output.vals, by=first)
        solutions["w$n"] = [format(value_precision_fe, real(v)), format(value_precision_fe, imag(v))]
    end
    for (n, v) in sort(output.f_output.vals, by=first)
        if isinf(v)
            solutions["f$n"] = ["Inf", "Inf"]
        else
            solutions["f$n"] = [format(value_precision_fe, real(v)), format(value_precision_fe, imag(v))]
        end

    end
    if use_λ
        solutions["λ"] = [format(value_precision_fe, real(output.λ_val)), format(value_precision_fe, imag(output.λ_val))]
    else
        solutions["c"] = [format(value_precision_fe, real(output.λ_val)), format(value_precision_fe, imag(output.λ_val))]
    end

    min_polys = Dict{String,Any}()
    for (n, v) in sort(output.b_output.vals, by=first)
        min_polys["b$n"] = 0
    end
    for (n, v) in sort(output.w_output.vals, by=first)
        min_polys["w$n"] = 0
    end
    for (n, v) in sort(output.f_output.vals, by=first)
        min_polys["f$n"] = 0
    end
    if use_λ
        min_polys["λ"] = 0
    else
        min_polys["c"] = 0
    end

    belyi_polys = Dict(
        "B_poly" => format_equation_string(output.B_poly, use_λ),
        "W_poly" => format_equation_string(output.W_poly, use_λ),
        "F_poly" => format_equation_string(output.F_poly, use_λ)
    )

    toml_dict = Dict(
        "info" => info,
        "valencies" => valencies,
        "equations" => equations,
        "solutions" => solutions,
        "min_polys" => min_polys,
        "belyi_polys" => belyi_polys
    )

    open(folder_path * "solns/$(file).toml", "w") do io
        TOML.print(io, toml_dict)
    end
end

function get_name(B::Vector{<:Integer}, W::Vector{<:Integer}, F::Vector{<:Integer})
    return "B-$(join(B, "-"))_W-$(join(W, "-"))_F-$(join(F, "-"))"
end