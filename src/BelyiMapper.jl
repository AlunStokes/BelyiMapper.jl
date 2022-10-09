# Copyright (c) 2022 Alun Cennyth Stokes
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module BelyiMapper

export generate_dessin,
    get_name

using Revise
using HomotopyContinuation

using Formatting
import JSON
import TOML

include("./dessin_utils.jl")
include("./generate_dessin.jl")

end