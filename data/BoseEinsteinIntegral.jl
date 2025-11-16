using QuadGK
using Printf

# ----------------------------------------------------------
# General integral:
#   I3(y0,y) = ∫_{y0}^{y} x^3/(e^x-1) dx
#
# Note: integral between y0 and y = full(0,y) - full(0,y0)
# ----------------------------------------------------------

# ζ(4) = π⁴ / 90
const ZETA4 = pi^4/90

# ----------------------------------------------------------
# 1) Quadrature evaluation between y0 and y
# ----------------------------------------------------------
function I3_quad(y0, y)
    f(x) = x^3 / (exp(x) - 1)
    val, _ = quadgk(f, y0, y, rtol=1e-12, atol=1e-14)
    return val
end

# ----------------------------------------------------------
# 2) Series evaluation between y0 and y
#
# Using:
#   ∫₀^y x^3/(e^x-1) dx = 6ζ(4) − Σ e^{-k y} (...)  
#
# So the integral between y0 and y is:
#   I3_series(y0,y) = I3_series_full(0,y) − I3_series_full(0,y0)
# ----------------------------------------------------------

function I3_series_full(y; tol=1e-12, kmax=10_000)
    S = 6 * ZETA4
    for k in 1:kmax
        ek = exp(-k*y)
        term = ek * (y^3/k + 3y^2/k^2 + 6y/k^3 + 6/k^4)
        S -= term
        if abs(term) < tol
            return S
        end
    end
    @warn "Series did not converge before kmax"
    return S
end

function I3_series(y0, y; tol=1e-12)
    return I3_series_full(y; tol=tol) - I3_series_full(y0; tol=tol)
end

# ----------------------------------------------------------
# Comparison + wall clock timing
# ----------------------------------------------------------

ys  = 0.1:0.1:1.1
y0  = 0.05        # lower bound (example)

println("Comparison of I3(y0,y) between y0 = $y0 and y ∈ [0.1, 1.1]")
println("-----------------------------------------------------------------------------------------")
@printf("%6s   %18s   %18s   %12s   %12s\n",
        "y", "quadgk", "series", "t_quad (s)", "t_series (s)")
println("-----------------------------------------------------------------------------------------")

for y in ys
    tq = @elapsed i_quad   = I3_quad(y0, y)
    ts = @elapsed i_series = I3_series(y0, y)

    @printf("%6.2f   %18.12e   %18.12e   %12.4e   %12.4e\n",
            y, i_quad, i_series, tq, ts)
end

println("-----------------------------------------------------------------------------------------")
