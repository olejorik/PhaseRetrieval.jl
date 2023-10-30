
abstract type Domain{N} end

"""
    PDplan(plan, diversity)

Construct `P`upil `D`iversity plan, which, if multiplied by the array `a` of proper dimensions, computes ``fft( diversity *a)``
"""
struct PDplan
    plan
    diversity::Array{<:Number}
end

PDplan(diversity::Array{<:Number}) = PDplan(plan_fft(diversity), diversity)

*(pd::PDplan, a::Array{<:Number}) = pd.plan * (pd.diversity .* a)
# export Base.:*
