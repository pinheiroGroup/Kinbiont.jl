
using MLJ
using SymbolicRegression
using Zygote 
using DelimitedFiles
jm_res_test = readdlm("E:/Lavoro/Monod_AA_res/ODE/exp_4/ODE_exp_4_parameters_aHPM.csv", ',')
annotation_test = CSV.File("E:/Lavoro/Monod_AA_res/Monod_AA_detection/exp_4/annotation.csv")
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:V1],annotation_test[:V3])
jmaki_results =jm_res_test
function downstream_symbolic_regression(jmaki_results,
    feature_matrix,
    row_to_learn;
    average_replicate = true,
    do_plots == true, 
    binary_operators=[+, *, /],
    unary_operators=nothing,
    batching=false,
    batch_size=1000,
    output_folder="/res/",
    output_file="output.txt",
    niterations=100,
    enable_autodiff =true,
    should_simplify=true,
    should_optimize_constants=false)

    mkpath(output_folder)

    names_of_the_wells_res = jmaki_results[3, 2:end]
    names_of_the_wells_annotation = feature_matrix[1:end, 1]
    wells_to_use = intersect(names_of_the_wells_res, names_of_the_wells_annotation)


    index_res = Int
    index_annotation = Int




    for i in wells_to_use
       if i == wells_to_use[1]
        index_res = findfirst(names_of_the_wells_res .== i )
        index_annotation = findfirst(names_of_the_wells_annotation .== i )
    
       else

        index_res = hcat(index_res,findfirst(names_of_the_wells_res .== i ))
        index_annotation = hcat(index_annotation,findfirst(names_of_the_wells_annotation .== i ))

       end
      

    end

    # order results and annotation by well in the same order






        ML_model = SRRegressor(
            niterations=niterations,
            binary_operators=binary_operators,
           # unary_operators=unary_operators,
            should_simplify=should_simplify,
            should_optimize_constants=should_optimize_constants,
            batching=batching,
            batch_size=batch_size,
            output_file=output_file
        )



        output = convert.(Float64, jmaki_results[row_to_learn, index_res[1,:] .+1])

        predictors = convert.(Float64, feature_matrix[index_annotation[1,:],2:end])




        mach = machine(ML_model, predictors, output)
        MLJ.fit!(mach)
        # report the fit 
        res_gr = report(mach)
        res_1 = res_gr.equation_strings
        res_2 = res_gr.equations
        res_3 = res_gr.equations[res_gr.best_idx]
        res_4 = res_gr.losses
        res_5 = res_gr.complexities
        res_6 = res_gr.scores
        index_min_score = findall(res_6 .==minimum(res_gr.scores ))


         if size(feature_importances)[2]==2 & do_plots == true



           display( scatter(predictors,output))
            for k in 1:length( res_gr.equations)

             a1= MLJ.predict(mach, (data=predictors, idx=k))
                scatter!(predictors,a1)
            end

         elseif  size(feature_importances)[2]!=2 & do_plots == true
            
            println("Warning: the plots are possible only for 2D data (1 feature, 1 output)")
            
         end   


    return res_1, res_2, res_3, res_4, res_5, res_6,res_gr.equations[index_min_score]

end





if average_replicate == true

    list_replicate = unique(feature_matrix[:,2:end])
    new_data = Any 

    for replicate_temp in list_replicate

       names_replicate_temp =findall(
                    x -> x == replicate_temp,
                    feature_matrix[:,2:end],
                )
       for i in names_replicate_temp
        #to do
         

       end
          
      if 
        new_data = hcat(replicate_temp, replicate_mean)

      else
        new_data = hcat(new_data, replicate_mean)

      end

    end
end






function downstream_decision_tree(
)


    predictors = hcat(convert.(Float64, size_at_start[2:end]), convert.(Float64, ribo_at_start[2:end]))

    # create output variable 

    output = convert.(Float64, final_size[2:end])

    model = build_tree(output, predictors)

    DecisionTree.print_tree(model, 2)

end


