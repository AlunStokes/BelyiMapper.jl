# Copyright (c) 2022 Alun Cennyth Stokes
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Test, BelyiMapper

base_folder = "/tmp"
verbose = true
output_files = false
use_λ = false

B = [5, 4, 4]
W = [2, 2, 2, 2, 2, 2, 1]
F = [3, 3, 3, 2, 2]

b_dist_val = 0 + 0im
w_dist_val = -1 + 0im
f_dist_val = Inf * im

dessin_name = BelyiMapper.get_name(B, W, F)

println(dessin_name)


println("B = $B | $(sum(B))")
println("W = $W | $(sum(W))")
println("F = $F | $(sum(F)) \n")

solutions = BelyiMapper.generate_dessin(
    B,
    W,
    F,
    b_dist_val,
    w_dist_val,
    f_dist_val,
    dessin_name,
    base_folder,
    verbose,
    output_files,
    use_λ)

display(solutions)

@test size(solutions)[1] == 6

println("Test passed.")