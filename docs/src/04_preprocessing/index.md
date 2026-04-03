# [Preprocessing](@id preprocessing)

```@contents
Pages = ["index.md"]
Depth = 3
```

All preprocessing is configured through `FitOptions` and applied automatically
by `kinbiont_fit`, or explicitly by calling `preprocess(data, opts)` directly.

```julia
using Kinbiont

# Every field has a sensible default — only set what you need
opts = FitOptions(
    smooth               = true,
    smooth_method        = :rolling_avg,
    blank_subtraction    = true,
    blank_from_labels    = true,
    cut_stationary_phase = true,
)

# Apply preprocessing without fitting
processed = preprocess(data, opts)
# processed.curves  → smoothed, blank-subtracted, truncated matrix
# processed.times   → time points (may change if :gaussian + custom grid)
# processed.labels  → labels (may change if average_replicates=true)
```

---

## Smoothing

Kinbiont supports five smoothing methods. All are controlled by `smooth_method`:

| Method | Key option | Best for |
|---|---|---|
| `:lowess` (default) | `lowess_frac` | General purpose; adapts to local density |
| `:rolling_avg` | `smooth_pt_avg` | Fast; uniform noise; fixed window |
| `:gaussian` | `gaussian_h_mult`, `gaussian_time_grid` | When you need a custom output time grid |
| `:boxcar` | `boxcar_window` | Symmetric window; preserves all time points |
| `:none` | — | Disable smoothing explicitly |

```julia
using Kinbiont, Plots, Random
Random.seed!(1)

# Simulate a noisy growth curve
sim  = ODE_sim("aHPM", [0.05], 0.0, 48.0, 1.0, [0.4, 0.1, 1.2, 2.0])
times = Float64.(sim.t)
od    = Float64.(reduce(hcat, sim.u)[1, :]) .+ 0.01 .* randn(length(times))
data  = GrowthData(reshape(od, 1, length(od)), times, ["well_A1"])

methods = [
    (:lowess,      FitOptions(smooth=true, smooth_method=:lowess,      lowess_frac=0.1)),
    (:rolling_avg, FitOptions(smooth=true, smooth_method=:rolling_avg, smooth_pt_avg=7)),
    (:gaussian,    FitOptions(smooth=true, smooth_method=:gaussian,    gaussian_h_mult=2.0)),
    (:boxcar,      FitOptions(smooth=true, smooth_method=:boxcar,      boxcar_window=7)),
]

p = scatter(times, od; label="raw", color=:grey, markersize=2, alpha=0.5,
            xlabel="Time (h)", ylabel="OD", title="Smoothing methods")
for (name, opts) in methods
    proc = preprocess(data, opts)
    plot!(p, proc.times, proc.curves[1, :]; label=string(name), linewidth=2)
end
display(p)
```

!!! note ":gaussian with a custom time grid"
    ```julia
    grid = collect(0.0:2.0:48.0)   # evenly spaced output times
    opts = FitOptions(smooth=true, smooth_method=:gaussian,
                      gaussian_time_grid=grid)
    proc = preprocess(data, opts)
    length(proc.times) == length(grid)   # true
    ```

---

## Blank subtraction

Subtract the blank OD from all curves. Two modes:

**Mode 1 — constant value** (no annotation file needed):

```julia
opts = FitOptions(blank_subtraction=true, blank_value=0.01)
```

**Mode 2 — automatic from labels** (requires wells labelled `"b"` in the annotation):

```julia
using Kinbiont, CSV

# Load data with annotation labels
data_raw = GrowthData("data_examples/ecoli_sample.csv")
ann      = CSV.File("data_examples/ecoli_annotation.csv"; header=false)
data     = GrowthData(data_raw.curves, data_raw.times, String.(ann[:Column2]))
# data.labels[1] is now "b" (Curve11728 is the blank well)

opts = FitOptions(blank_subtraction=true, blank_from_labels=true)
proc = preprocess(data, opts)
```

!!! warning "Order matters"
    Blank resolution from labels happens **before** replicate averaging, because
    averaging drops `"b"`-labelled wells from the matrix. Always load annotation
    labels and call `blank_from_labels=true` before setting `average_replicates=true`.

---

## Replicate averaging

When multiple wells measure the same biological condition, average them into one curve:

```julia
using Kinbiont, CSV

# Annotation format (no header): well_name,label
# Curve11729,condition_A
# Curve11730,condition_A   ← these two will be averaged
# Curve11728,b             ← blank: excluded from output
# Curve11737,X             ← discard: excluded from output

data_raw = GrowthData("data_examples/ecoli_sample.csv")
ann      = CSV.File("data_examples/ecoli_annotation.csv"; header=false)
data     = GrowthData(data_raw.curves, data_raw.times, String.(ann[:Column2]))

opts = FitOptions(
    average_replicates = true,
    blank_subtraction  = true,
    blank_from_labels  = true,
)
proc = preprocess(data, opts)
# proc.labels now contains only unique non-b/non-X labels
```

---

## Negative value correction

After blank subtraction, some values may go negative. Three correction methods:

```julia
# :remove — drop time points where OD < 0 (default)
opts = FitOptions(correct_negatives=true, negative_method=:remove)

# :thr_correction — floor all values at negative_threshold
opts = FitOptions(correct_negatives=true, negative_method=:thr_correction,
                  negative_threshold=0.001)

# :blank_correction — replace negatives with the blank mean
opts = FitOptions(correct_negatives=true, negative_method=:blank_correction)
```

---

## Stationary phase cutting

Truncate each curve at the onset of stationary phase before fitting. Kinbiont
detects the cutoff by finding where the specific growth rate (SGR) falls below
`stationary_percentile_thr × max(SGR)`, then snaps to the nearest OD peak.

```julia
using Kinbiont, Plots

# Simulate a triphasic curve with a clear stationary phase
sim  = ODE_sim("triple_piecewise_adjusted_logistic", [0.1], 0.0, 500.0, 15.0,
               [0.06, 1.0, 200, 0.5, 0.001, 450, -0.0002])
times = Float64.(sim.t)
od    = Float64.(reduce(hcat, sim.u)[1, :])
data  = GrowthData(reshape(od, 1, length(od)), times, ["well_A1"])

opts = FitOptions(
    cut_stationary_phase            = true,
    stationary_percentile_thr       = 0.05,   # SGR threshold (5% of max)
    stationary_pt_smooth_derivative = 10,      # smoothing window for SGR
    stationary_win_size             = 5,       # look-ahead for OD peak
    stationary_thr_od               = 0.02,   # ignore time points below this OD
)
proc = preprocess(data, opts)

p = scatter(times, od;       label="full curve",    markersize=2, alpha=0.5,
            xlabel="Time", ylabel="OD")
scatter!(p, proc.times, proc.curves[1, :]; label="after cutoff", markersize=3)
display(p)
println("Kept $(length(proc.times)) / $(length(times)) time points")
```

---

## Scattering correction

Correct OD values using an instrument-specific calibration curve (maps raw OD
to true OD measured by an independent method). See `data_examples/cal_curve_example.csv`
for the expected format.

```julia
opts = FitOptions(
    scattering_correction = true,
    calibration_file      = "data_examples/cal_curve_example.csv",
    scattering_method     = :interpolation,   # or :exp_fit
)
proc = preprocess(data, opts)
```

!!! warning "Instrument-specific"
    The calibration file must be generated on the same microplate reader used to
    collect the data. Do not reuse calibration files across instruments.

---

## Real-data example

Apply a full preprocessing pipeline to the *E. coli* Keio knockout dataset
(Reding-Roman et al., *Nature Scientific Data* 2026,
doi:10.1038/s41597-026-07075-9):

```julia
using Kinbiont, Plots, CSV

data_raw = GrowthData("data_examples/ecoli_sample.csv")
ann      = CSV.File("data_examples/ecoli_annotation.csv"; header=false)
data     = GrowthData(data_raw.curves, data_raw.times, String.(ann[:Column2]))

opts = FitOptions(
    smooth               = true,
    smooth_method        = :rolling_avg,
    smooth_pt_avg        = 7,
    blank_subtraction    = true,
    blank_from_labels    = true,
    correct_negatives    = true,
    negative_method      = :thr_correction,
    negative_threshold   = 0.001,
    cut_stationary_phase = true,
)

proc = preprocess(data, opts)

p = plot(layout=(1,2), size=(800,320),
         title=["Raw" "Preprocessed"],
         xlabel="Time (h)", ylabel="OD")
for i in axes(data.curves, 1)
    plot!(p[1], data.times, data.curves[i, :]; label=false, alpha=0.6)
    plot!(p[2], proc.times, proc.curves[i, :]; label=false, alpha=0.6)
end
display(p)
```
