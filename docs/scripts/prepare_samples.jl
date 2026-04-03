# docs/scripts/prepare_samples.jl
# Run from repo root: julia --project=. docs/scripts/prepare_samples.jl

using CSV, DataFrames, XLSX

CURVES_CSV   = "../curves-data/all_curves.csv"
BATCH_FIT    = "../chemical-media-analysis/results/batch_fit_results.csv"
MEDIUM_XLSX  = "../chemical-media-dataset/xlsx_raw/BW25113_Medium composition.xlsx"
OUT_DIR      = "data_examples"

# ── 1. E. coli sample (10 wells, all time points) ────────────────────────────
println("Extracting ecoli_sample.csv …")
raw = CSV.read(CURVES_CSV, DataFrame; header=true)
# Column 1 is unnamed row-index; column 2 is "Time（h）"; rest are curves
rename!(raw, names(raw)[1] => "index", names(raw)[2] => "Time_h")
wells = ["Curve11728","Curve11729","Curve11730","Curve11731","Curve11732",
         "Curve11733","Curve11734","Curve11735","Curve11736","Curve11737"]
sample = raw[:, vcat(["Time_h"], wells)]
CSV.write(joinpath(OUT_DIR, "ecoli_sample.csv"), sample)
println("  → $(nrow(sample)) rows × $(ncol(sample)) cols")

# ── 2. Annotation for ecoli_sample ───────────────────────────────────────────
println("Writing ecoli_annotation.csv …")
ann = DataFrame(
    well  = wells,
    label = ["b", "condition_A", "condition_A", "condition_B", "condition_B",
             "condition_C", "condition_C", "condition_D", "condition_D", "X"],
)
CSV.write(joinpath(OUT_DIR, "ecoli_annotation.csv"), ann; header=false)
println("  → $(nrow(ann)) rows")

# ── 3. Batch fit sample (200 converged rows) ──────────────────────────────────
println("Extracting batch_fit_sample.csv …")
fit_df = CSV.read(BATCH_FIT, DataFrame)
conv   = filter(r -> r.converged, fit_df)
samp   = conv[1:min(200, nrow(conv)), :]
CSV.write(joinpath(OUT_DIR, "batch_fit_sample.csv"), samp)
println("  → $(nrow(samp)) rows")

# ── 4. Medium composition (matched to batch_fit_sample, 5 compounds) ──────────
println("Extracting medium_composition_sample.csv …")
xf   = XLSX.readxlsx(MEDIUM_XLSX)
sh   = xf[XLSX.sheetnames(xf)[1]]
data = sh[:]
raw_hdrs = [string(h) for h in data[1, :]]

# Sanitise headers: replace newlines, spaces, parentheses with underscores
clean_hdrs = [replace(h, r"[\n\s\(\)]+" => "_") for h in raw_hdrs]

# Build DataFrame
med_df = DataFrame(
    [Symbol(clean_hdrs[j]) => data[2:end, j] for j in 1:length(clean_hdrs)]
)

# The first column is the curve Label
label_col = Symbol(clean_hdrs[1])
rename!(med_df, label_col => :Label)
med_df.Label = string.(med_df.Label)

println("  Medium cols: ", names(med_df)[1:min(8, ncol(med_df))])

# Find compound columns: try common names after sanitisation
candidates = [:Glucose_mM_, :K2HPO4_mM_, :KH2PO4_mM_, :Na2HPO4_mM_, :_NH4_2SO4_mM_]
compound_cols = filter(c -> c in Symbol.(names(med_df)), candidates)

if isempty(compound_cols)
    # Fallback: just take columns 4-8 (skipping Label, Assay_ID, Condition_ID)
    all_cols = Symbol.(names(med_df))
    compound_cols = all_cols[4:min(8, length(all_cols))]
    @warn "Expected compound columns not found; using fallback: $compound_cols"
end

select!(med_df, vcat([:Label], compound_cols))

# Rename to clean names regardless of fallback
target_names = [:Glucose_mM, :K2HPO4_mM, :KH2PO4_mM, :Na2HPO4_mM, :NH4_2SO4_mM]
for (i, old) in enumerate(compound_cols)
    i > length(target_names) && break
    rename!(med_df, old => target_names[i])
end

# Keep only labels present in batch_fit_sample
keep = Set(samp.label)
med_sample = filter(r -> r.Label in keep, med_df)
CSV.write(joinpath(OUT_DIR, "medium_composition_sample.csv"), med_sample)
println("  → $(nrow(med_sample)) rows × $(ncol(med_sample)) cols")

println("\nDone. Files written to $(OUT_DIR)/")
