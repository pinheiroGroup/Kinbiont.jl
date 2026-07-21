@testset "Preprocessing" begin
    # Small synthetic data: 5 curves × 10 timepoints
    Random.seed!(42)
    times  = collect(0.0:1.0:9.0)
    curves = abs.(rand(5, 10)) .+ 0.1   # all positive
    data   = GrowthData(curves, times, ["c$i" for i in 1:5])

    @testset "No-op (all steps disabled)" begin
        processed = preprocess(data, FitOptions())
        @test processed.curves ≈ data.curves
        @test processed.times  == data.times
        @test processed.clusters === nothing
    end

    @testset "Blank subtraction" begin
        opts = FitOptions(blank_subtraction=true, blank_value=0.05)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves .- 0.05
    end

    @testset "Auto blank detection from labels" begin
        blank_od  = 0.08
        # Build a dataset with two blank wells (label "b") and three signal wells
        n_tp      = length(times)
        blank_row = fill(blank_od, n_tp)
        all_curves = vcat(curves, blank_row', blank_row')
        all_labels = ["c1","c2","c3","c4","c5","b","b"]
        data_with_blanks = GrowthData(all_curves, times, all_labels)

        opts = FitOptions(blank_subtraction=true, blank_from_labels=true)
        processed = preprocess(data_with_blanks, opts)
        # every curve should be shifted by the mean blank (= blank_od)
        @test processed.curves ≈ data_with_blanks.curves .- blank_od
    end

    @testset "Auto blank: warns and uses 0.0 when no 'b' labels" begin
        opts = FitOptions(blank_subtraction=true, blank_from_labels=true)
        @test_warn r"no wells labelled" preprocess(data, opts)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves   # subtracting 0.0
    end
  
    @testset "Replicate averaging collapses duplicate labels" begin
        # 4 curves: A appears twice, B appears twice
        rep_curves = vcat(curves[1:2, :], curves[3:4, :])
        rep_labels = ["A", "A", "B", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test size(processed.curves, 1) == 2          # collapsed to 2 unique labels
        @test processed.labels == ["A", "B"]
        @test processed.curves[1, :] ≈ mean(rep_curves[1:2, :], dims=1)[:]
        @test processed.curves[2, :] ≈ mean(rep_curves[3:4, :], dims=1)[:]
    end

    @testset "Replicate averaging drops 'b' and 'X' wells" begin
        rep_curves = vcat(curves[1:2, :], curves[3:4, :], curves[5:5, :])
        rep_labels = ["A", "A", "b", "X", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test processed.labels == ["A", "B"]
        @test size(processed.curves, 1) == 2
    end

    @testset "Replicate averaging is a no-op when disabled" begin
        opts = FitOptions(average_replicates=false)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves
        @test processed.labels == data.labels
    end

    @testset "Clustering assigns cluster ids, centroids and WCSS" begin
        n_k  = 2
        opts = FitOptions(cluster=true, n_clusters=n_k, cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters  isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
        # centroid matrix: one row per cluster, one column per timepoint
        @test processed.centroids isa Matrix{Float64}
        @test size(processed.centroids) == (n_k, length(data.times))
        # centroids are in z-normalised space: each non-empty row has ~zero mean
        for k in 1:n_k
            if any(==(k), processed.clusters)
                @test abs(mean(processed.centroids[k, :])) < 0.1
            end
        end
        # WCSS is a non-negative scalar
        @test processed.wcss isa Float64
        @test processed.wcss >= 0.0
    end

    @testset "WCSS equals total SSE to assigned cluster centroid (all methods)" begin
        # z-score each row the way _cluster does (no smoothing → curves unchanged);
        # StatsBase.zscore uses the corrected (n-1) std, matching (v-mean)/std here.
        _zrow(v) = (v .- mean(v)) ./ std(v)
        Z = reduce(vcat, [reshape(_zrow(data.curves[i, :]), 1, :) for i in 1:size(data.curves, 1)])
        _sse_to_centroid(lbls) = begin
            total = 0.0
            for lbl in unique(lbls)
                idx = findall(==(lbl), lbls)
                sub = Z[idx, :]
                c   = vec(mean(sub, dims=1))
                for i in axes(sub, 1)
                    total += sum((sub[i, :] .- c) .^ 2)
                end
            end
            total
        end
        # For k-medoids in particular this guards against the old bug where WCSS was
        # a sum of *unsquared* distances to medoids (result.totalcost).
        for method in (:kmeans, :kmedoids, :hclust)
            opts = FitOptions(cluster=true, n_clusters=2, cluster_method=method,
                              cluster_trend_test=false)
            p = preprocess(data, opts)
            @test isapprox(p.wcss, _sse_to_centroid(p.clusters); rtol=1e-8)
        end
    end

    @testset "WCSS counts set-aside flat curves (trend/sentinel path)" begin
        _zrow(v) = (v .- mean(v)) ./ std(v)
        tt = collect(0.0:5.0)
        cc = [0.10  0.20  0.35  0.55  0.80  1.10;    # growing
              0.12  0.25  0.40  0.60  0.85  1.15;    # growing
              0.100 0.101 0.099 0.100 0.102 0.101;   # flat
              0.050 0.051 0.049 0.050 0.052 0.051]   # flat
        dd = GrowthData(cc, tt, ["g1", "g2", "f1", "f2"])
        opts = FitOptions(cluster=true, n_clusters=2, cluster_trend_test=true,
                          cluster_trend_p_thr=0.05)
        p = preprocess(dd, opts)
        Z = reduce(vcat, [reshape(_zrow(cc[i, :]), 1, :) for i in 1:size(cc, 1)])
        expected = 0.0
        for lbl in unique(p.clusters)
            idx = findall(==(lbl), p.clusters)
            sub = Z[idx, :]
            c   = vec(mean(sub, dims=1))
            for i in axes(sub, 1)
                expected += sum((sub[i, :] .- c) .^ 2)
            end
        end
        # WCSS is over ALL 4 curves (both dynamic and the reserved flat sentinel),
        # not just the dynamic subset that the algorithm clustered.
        @test isapprox(p.wcss, expected; rtol=1e-8)
    end

    @testset "Replicate averaging collapses duplicate labels" begin
        # 4 curves: A appears twice, B appears twice
        rep_curves = vcat(curves[1:2, :], curves[3:4, :])
        rep_labels = ["A", "A", "B", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test size(processed.curves, 1) == 2          # collapsed to 2 unique labels
        @test processed.labels == ["A", "B"]
        @test processed.curves[1, :] ≈ mean(rep_curves[1:2, :], dims=1)[:]
        @test processed.curves[2, :] ≈ mean(rep_curves[3:4, :], dims=1)[:]
    end

    @testset "Replicate averaging drops 'b' and 'X' wells" begin
        rep_curves = vcat(curves[1:2, :], curves[3:4, :], curves[5:5, :])
        rep_labels = ["A", "A", "b", "X", "B"]
        rep_data   = GrowthData(rep_curves, times, rep_labels)

        opts = FitOptions(average_replicates=true)
        processed = preprocess(rep_data, opts)

        @test processed.labels == ["A", "B"]
        @test size(processed.curves, 1) == 2
    end

    @testset "Replicate averaging is a no-op when disabled" begin
        opts = FitOptions(average_replicates=false)
        processed = preprocess(data, opts)
        @test processed.curves ≈ data.curves
        @test processed.labels == data.labels
    end

    @testset "Exponential prototype cluster labels are valid positive integers" begin
        opts = FitOptions(cluster=true, n_clusters=3,
                          cluster_trend_test=false, cluster_exp_prototype=true)
        processed = preprocess(data, opts)
        # exp prototype may add one extra cluster beyond n_clusters if that slot
        # was already used by k-means, so upper bound is n_clusters + 1
        @test all(1 .<= processed.clusters .<= opts.n_clusters + 1)
    end

    @testset "Exponential prototype reassigns a clearly exponential curve" begin
        # Build a plate where one curve is strongly exponential
        t = data.times
        exp_curve = 0.01 .* exp.(0.8 .* t)
        mixed = vcat(data.curves, reshape(exp_curve, 1, :))
        mixed_data = GrowthData(mixed, t, ["c$i" for i in 1:6])
        opts = FitOptions(cluster=true, n_clusters=3,
                          cluster_trend_test=false, cluster_exp_prototype=true)
        processed = preprocess(mixed_data, opts)
        # The exponential prototype always occupies the highest label index;
        # the clearly exponential curve (last row) should be assigned to it.
        @test processed.clusters[end] == maximum(processed.clusters)
    end

    @testset "WCSS decreases as n_clusters increases" begin
        # More clusters → lower total SSE (elbow-plot property)
        wcss_vals = map(2:4) do k
            opts = FitOptions(cluster=true, n_clusters=k, cluster_trend_test=false)
            preprocess(data, opts).wcss
        end
        @test issorted(wcss_vals, rev=true)
    end

    @testset "Clustering labels always within 1..n_clusters (cluster_trend_test=true)" begin
        n_k = 3
        opts = FitOptions(cluster=true, n_clusters=n_k, cluster_trend_test=true)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test all(1 .<= processed.clusters .<= n_k)
        @test size(processed.centroids) == (n_k, length(data.times))
    end

    @testset "cluster_trend_test: lag+growth+stationary not mislabeled as flat" begin
        n_lag  = 15
        n_grow = 3
        n_stat = 15
        od_low  = 0.05
        od_high = 1.5
        grow_vals = [od_low + (od_high - od_low) * i / (n_grow + 1) for i in 1:n_grow]
        growing_curve = vcat(fill(od_low, n_lag), grow_vals, fill(od_high, n_stat))
        flat_curve    = fill(od_low, n_lag + n_grow + n_stat)
        n_tp      = n_lag + n_grow + n_stat
        tp_times  = collect(0.0:(n_tp - 1))
        tp_curves = Matrix(hcat(growing_curve, flat_curve)')  # 2 × n_tp
        tp_data   = GrowthData(tp_curves, tp_times, ["growing", "flat"])

        n_k  = 2
        opts = FitOptions(cluster=true, n_clusters=n_k, cluster_trend_test=true)
        processed = preprocess(tp_data, opts)

        # The flat curve must receive the flat label (n_clusters)
        @test processed.clusters[2] == n_k
        # The growing curve must NOT receive the flat label
        @test processed.clusters[1] != n_k
    end

    @testset "cluster_trend_test does not reserve sentinel when no flat curves are found" begin
        trend_times = collect(0.0:9.0)
        trend_curves = [0.10 0.20 0.32 0.48 0.67 0.90 1.12 1.35 1.55 1.75;
                        0.12 0.24 0.38 0.56 0.78 1.00 1.25 1.48 1.70 1.92;
                        1.80 1.62 1.44 1.26 1.08 0.90 0.72 0.54 0.36 0.18;
                        1.95 1.75 1.55 1.35 1.15 0.95 0.75 0.55 0.35 0.15]
        trend_data = GrowthData(trend_curves, trend_times, ["up1", "up2", "down1", "down2"])

        opts = FitOptions(cluster=true, n_clusters=2, cluster_trend_test=true)
        processed = preprocess(trend_data, opts)

        @test length(unique(processed.clusters)) == 2
        @test all(1 .<= processed.clusters .<= 2)
        @test all(count(==(k), processed.clusters) > 0 for k in 1:2)
    end

    @testset "cluster_trend_test reserves sentinel for flat curves" begin
        trend_times = collect(0.0:5.0)
        trend_curves = [0.100 0.101 0.099 0.100 0.102 0.101;
                        0.110 0.109 0.112 0.110 0.111 0.109;
                        0.100 0.160 0.280 0.450 0.650 0.820;
                        0.120 0.190 0.330 0.520 0.740 0.930;
                        0.140 0.220 0.390 0.610 0.850 1.050;
                        0.160 0.250 0.450 0.700 0.970 1.200]
        trend_data = GrowthData(trend_curves, trend_times,
                                ["F1", "F2", "G1", "G2", "G3", "G4"])

        opts = FitOptions(cluster=true, n_clusters=2,
                          cluster_trend_test=true, cluster_trend_p_thr=0.05)
        processed = preprocess(trend_data, opts)

        @test all(1 .<= processed.clusters .<= 2)
        @test processed.clusters[1:2] == [2, 2]
        @test all(processed.clusters[3:6] .!= 2)
    end

    @testset "cluster_trend_test leaves k=1 as all-data baseline" begin
        trend_times = collect(0.0:5.0)
        trend_curves = [0.100 0.101 0.099 0.100 0.102 0.101;
                        0.110 0.109 0.112 0.110 0.111 0.109;
                        0.100 0.160 0.280 0.450 0.650 0.820;
                        0.120 0.190 0.330 0.520 0.740 0.930]
        trend_data = GrowthData(trend_curves, trend_times, ["F1", "F2", "G1", "G2"])

        baseline = preprocess(trend_data, FitOptions(cluster=true, n_clusters=1,
                                                     cluster_trend_test=false))
        processed = preprocess(trend_data, FitOptions(cluster=true, n_clusters=1,
                                                      cluster_trend_test=true))

        @test processed.clusters == fill(1, size(trend_curves, 1))
        @test processed.wcss ≈ baseline.wcss
        @test processed.wcss > 0
    end

    @testset "cluster_trend_p_thr controls post-hoc flat reassignment" begin
        weak_trend = [0.1 + 0.01 * t + 0.05 * sin(1.7 * t) for t in times]
        flat_curve = fill(0.1, length(times))
        trend_curves = Matrix(hcat(weak_trend, flat_curve)')
        start_labels = [1, 1]

        lower_thr = Kinbiont._apply_trend_labels(trend_curves, times, start_labels, 2, 0.01)
        default_thr = Kinbiont._apply_trend_labels(trend_curves, times, start_labels, 2, 0.05)

        @test lower_thr[1] == 2
        @test default_thr[1] == 1
        @test lower_thr[2] == 2
        @test default_thr[2] == 2
    end

    @testset "Constant pre-screening keeps labels within 1..n_clusters" begin
        # Mix flat and growing curves so pre-screening has something to detect
        flat_curves = hcat(fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3), fill(0.1, 3))
        mixed = vcat(data.curves, flat_curves)
        mixed_data = GrowthData(mixed, data.times, ["c$i" for i in 1:8])
        n_k  = 3
        opts = FitOptions(cluster=true, n_clusters=n_k,
                          cluster_prescreen_constant=true, cluster_trend_test=false)
        processed = preprocess(mixed_data, opts)
        @test all(1 .<= processed.clusters .<= n_k)
        @test size(processed.centroids) == (n_k, length(data.times))
        # flat curves should be assigned to the last label
        @test all(processed.clusters[6:8] .== n_k)
    end

    @testset "Constant pre-screening leaves k=1 as all-data baseline" begin
        flat_curves = hcat(fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3),
                           fill(0.1, 3), fill(0.1, 3), fill(0.1, 3), fill(0.1, 3))
        mixed = vcat(data.curves, flat_curves)
        mixed_data = GrowthData(mixed, data.times, ["c$i" for i in 1:8])

        baseline = preprocess(mixed_data, FitOptions(cluster=true, n_clusters=1,
                                                     cluster_trend_test=false))
        processed = preprocess(mixed_data, FitOptions(cluster=true, n_clusters=1,
                                                      cluster_prescreen_constant=true,
                                                      cluster_trend_test=false))

        @test processed.clusters == fill(1, size(mixed, 1))
        @test processed.wcss ≈ baseline.wcss
        @test processed.wcss > 0
    end

    @testset "Constant pre-screening does not reserve sentinel when no constant curves are found" begin
        prescreen_times = collect(0.0:9.0)
        prescreen_curves = [0.10 0.20 0.32 0.48 0.67 0.90 1.12 1.35 1.55 1.75;
                            0.12 0.24 0.38 0.56 0.78 1.00 1.25 1.48 1.70 1.92;
                            1.80 1.62 1.44 1.26 1.08 0.90 0.72 0.54 0.36 0.18;
                            1.95 1.75 1.55 1.35 1.15 0.95 0.75 0.55 0.35 0.15]
        prescreen_data = GrowthData(prescreen_curves, prescreen_times,
                                    ["up1", "up2", "down1", "down2"])

        opts = FitOptions(cluster=true, n_clusters=2,
                          cluster_prescreen_constant=true, cluster_trend_test=false)
        processed = preprocess(prescreen_data, opts)

        @test length(unique(processed.clusters)) == 2
        @test all(1 .<= processed.clusters .<= 2)
        @test all(count(==(k), processed.clusters) > 0 for k in 1:2)
    end

    # -----------------------------------------------------------------------
    # cluster_method variants
    # -----------------------------------------------------------------------

    @testset "cluster_method=:kmedoids produces valid labels" begin
        opts = FitOptions(cluster=true, n_clusters=3, cluster_method=:kmedoids,
                          cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
        @test all(1 .<= processed.clusters .<= 3)
    end

    @testset "cluster_method=:hclust produces valid labels" begin
        opts = FitOptions(cluster=true, n_clusters=3, cluster_method=:hclust,
                          cluster_hclust_linkage=:ward, cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
        @test all(1 .<= processed.clusters .<= 3)
    end

    @testset "cluster_method=:dbscan produces non-negative labels" begin
        opts = FitOptions(cluster=true, cluster_method=:dbscan,
                          cluster_dbscan_eps=2.0, cluster_dbscan_minpts=2,
                          cluster_trend_test=false)
        processed = preprocess(data, opts)
        @test processed.clusters isa Vector{Int}
        @test length(processed.clusters) == size(data.curves, 1)
        @test all(processed.clusters .>= 0)   # 0 = noise
    end

    @testset "cluster_method=:kmedoids with prescreen assigns sentinel" begin
        flat_curves = repeat([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1], 2)
        mixed_data  = GrowthData(vcat(data.curves, flat_curves), data.times,
                                  ["c$i" for i in 1:7])
        opts = FitOptions(cluster=true, n_clusters=3, cluster_method=:kmedoids,
                          cluster_prescreen_constant=true, cluster_trend_test=false)
        processed = preprocess(mixed_data, opts)
        @test all(processed.clusters[6:7] .== 3)
        @test all(processed.clusters[1:5] .!= 3)
    end

    @testset "cluster_method=:hclust with trend_test reassigns flat curves" begin
        times_test = collect(0.0:1.0:9.0)
        flat   = fill(0.3, 1, 10)
        growth = [range(0.0, 1.0; length=10)'; range(0.0, 1.5; length=10)']
        gd     = GrowthData(vcat(flat, growth), times_test, ["flat","g1","g2"])
        opts   = FitOptions(cluster=true, n_clusters=3, cluster_method=:hclust,
                            cluster_trend_test=true)
        processed = preprocess(gd, opts)
        @test processed.clusters[1] == maximum(processed.clusters)
    end

    @testset "Unknown cluster_method throws" begin
        opts = FitOptions(cluster=true, n_clusters=2, cluster_method=:bogus,
                          cluster_trend_test=false)
        @test_throws ErrorException preprocess(data, opts)
    end

    @testset "Pure function — original data not mutated" begin
        original_copy = copy(data.curves)
        _ = preprocess(data, FitOptions(blank_subtraction=true, blank_value=0.1))
        @test data.curves == original_copy
    end

    @testset "Smoothing (rolling_avg) does not error" begin
        opts = FitOptions(smooth=true, smooth_method=:rolling_avg, smooth_pt_avg=3)
        processed = preprocess(data, opts)
        # rolling_avg may shorten the time dimension; number of curves is preserved
        @test size(processed.curves, 1) == size(data.curves, 1)
        @test size(processed.curves, 2) <= size(data.curves, 2)
    end

    @testset "Smoothing (boxcar) preserves shape and time grid" begin
        opts = FitOptions(smooth=true, smooth_method=:boxcar, boxcar_window=3)
        processed = preprocess(data, opts)
        # boxcar is length-preserving: shape and times identical to input
        @test size(processed.curves) == size(data.curves)
        @test processed.times == data.times
        # smoothing must change at least some values
        @test processed.curves != data.curves
    end

    @testset "Smoothing (gaussian) keeps original times when no grid given" begin
        opts = FitOptions(smooth=true, smooth_method=:gaussian, gaussian_h_mult=2.0)
        processed = preprocess(data, opts)
        @test size(processed.curves) == size(data.curves)
        @test processed.times == data.times
    end

    @testset "Gaussian smoothing with custom time grid changes times and columns" begin
        new_times = collect(0.0:0.5:9.0)   # finer grid → 19 points vs original 10
        opts = FitOptions(
            smooth=true,
            smooth_method=:gaussian,
            gaussian_time_grid=new_times,
        )
        processed = preprocess(data, opts)
        @test processed.times == new_times
        @test size(processed.curves, 1) == size(data.curves, 1)   # same n_curves
        @test size(processed.curves, 2) == length(new_times)
    end

    @testset "Public clustering preprocessing helpers" begin
        raw_times = [collect(0.0:1.0:4.0), collect(1.0:1.0:5.0)]
        raw_curves = [[0.0, 1.0, 2.0, 3.0, 4.0], [10.0, 11.0, 12.0, 13.0, 14.0]]
        grid = common_time_grid(raw_times; n_grid=5, q_low=0.0, q_high=1.0)
        interp = interpolate_curves_to_grid(raw_curves, raw_times, grid)

        @test grid == collect(range(0.0, 5.0; length=5))
        @test size(interp) == (2, 5)
        @test interp[1, 1] == raw_curves[1][1]
        @test interp[2, end] == raw_curves[2][end]

        corrected = apply_blank_timeseries(copy(data.curves), fill(0.1, length(times));
                                           method=:pointbypoint)
        @test corrected ≈ data.curves .- 0.1

        dirty = copy(data.curves)
        dirty[1, 1] = NaN
        filled = fill_nonfinite_colmean(dirty)
        @test isfinite(filled[1, 1])

        blank_rows = vcat(fill(0.01, 1, length(times)),
                          [0.2 + 0.1 * t for t in times]')
        blanks = detect_blank_indices(blank_rows, times;
                                      flat_range_thr=0.005,
                                      od_percentile=0.5)
        @test blanks == [1]

        q = cluster_quality_indices(data.curves, [1, 1, 2, 2, 2])
        @test haskey(q, "silhouette_mean")
        @test haskey(q, "davies_bouldin")
    end

    @testset "apply_blank_timeseries handles NaN blanks" begin
        # Non-finite blank entries should not propagate into corrected curves
        # for :pointbypoint (matches :shift/:clip filtering behaviour).
        curves_local = [0.1 0.2 0.3; 0.4 0.5 0.6]
        blank_with_nan = [0.05, NaN, 0.05]
        corrected = apply_blank_timeseries(curves_local, blank_with_nan;
                                           method=:pointbypoint)
        @test all(isfinite, corrected)
        @test corrected[1, 1] ≈ 0.05
        @test corrected[1, 2] == 0.2   # NaN blank → no correction
        @test corrected[2, 3] ≈ 0.55
    end

    @testset "prepare_clustering_data from CSV" begin
        mktempdir() do dir
            csv_file = joinpath(dir, "input.csv")
            open(csv_file, "w") do io
                println(io, "time,well_blank,well_a,well_b")
                for t in 0.0:1.0:9.0
                    println(io, "$t,0.05,$(0.05 + 0.05 * t),$(0.1 + 0.05 * t)")
                end
            end

            gd = prepare_clustering_data(; csv_path=csv_file,
                                          auto_detect_blanks=true,
                                          subtract_blank=false,
                                          blank_od_percentile=0.5)
            @test gd isa GrowthData
            @test size(gd.curves, 1) == 2   # blank well auto-detected and removed
            @test "well_blank" ∉ gd.labels
            @test length(gd.times) == 10
        end
    end

    @testset "prepare_clustering_data drops fully-NaN curves" begin
        # Regression: previously a fully-NaN curve was fabricated as [0.0, 0.0]
        # and silently passed through. It should be dropped instead.
        mktempdir() do dir
            csv_file = joinpath(dir, "input.csv")
            open(csv_file, "w") do io
                println(io, "time,well_a,well_dead,well_b")
                for (i, t) in enumerate(0.0:1.0:4.0)
                    dead = "NaN"
                    println(io, "$t,$(0.1 * i),$dead,$(0.2 * i)")
                end
            end

            gd = prepare_clustering_data(; csv_path=csv_file,
                                          auto_detect_blanks=false)
            @test "well_dead" ∉ gd.labels
            @test size(gd.curves, 1) == 2
        end
    end
end

# 2.3 — negative_value_correction (old-API utility)
@testset "negative_value_correction" begin
    times = collect(0.0:1.0:9.0)
    od    = [-0.05, 0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0, 1.2, 1.3]
    data_mat = Matrix(transpose(hcat(times, od)))   # 2×10 old-API format

    corrected = negative_value_correction(data_mat, Float64[];
                                          method="thr_correction", thr_negative=0.01)
    @test all(corrected[2, :] .>= 0.0)
    @test size(corrected, 1) == 2
end

# 2.4 — specific_gr_evaluation on exponential data
@testset "specific_gr_evaluation" begin
    gr_true = 0.5
    t  = collect(0.0:0.5:10.0)
    od = 0.05 .* exp.(gr_true .* t)
    data_mat = Matrix(transpose(hcat(t, od)))   # 2×N old-API format

    sgr = specific_gr_evaluation(data_mat, 7)

    @test length(sgr) > 0
    @test all(isfinite.(sgr))
    # mean SGR should be close to the true growth rate
    @test abs(mean(sgr) - gr_true) < 0.15
end

# 2.6 — stationary phase cutting (exercised through kinbiont_fit)
@testset "Stationary phase cutting via kinbiont_fit" begin
    Random.seed!(42)
    # Dense logistic curve that clearly reaches a plateau by t=20
    sim   = ODE_sim("logistic", [0.05], 0.0, 30.0, 0.3, [0.5, 1.5])
    times = collect(0.0:0.3:30.0)
    curve = sim[1, :]
    data  = GrowthData(reshape(curve, 1, :), times, ["A"])

    spec      = ModelSpec([MODEL_REGISTRY["NL_logistic"]], [[1.5, 0.3, 0.0]])
    opts_cut  = FitOptions(cut_stationary_phase=true, loss="RE")
    opts_full = FitOptions(cut_stationary_phase=false, loss="RE")

    res_cut  = kinbiont_fit(data, spec, opts_cut)
    res_full = kinbiont_fit(data, spec, opts_full)

    # With cutting, the fitted curve should cover fewer time points
    @test length(res_cut[1].times) < length(res_full[1].times)
end
