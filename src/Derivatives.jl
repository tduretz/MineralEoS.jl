using Enzyme

function compute_density_derivative(EoS, P, dP, T, params, opts)
    Enzyme.autodiff(Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params; options=opts), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
end

function compute_density_derivative(EoS, P, dP, T, params)
    Enzyme.autodiff(Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
end