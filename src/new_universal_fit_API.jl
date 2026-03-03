using Kinbiont

"""

"""
function Kinbiont_fit(
    Kinbiont_data, # data array counld be a matrix A of any,  a path in string "path_1", an array of matrix and string array ["path_1",Matrix_1,...]
    Kinbiont_model, #  Kinbiont fit model data struct 
    Kinbiont_options_preprocessing,
    Kinbiont_options_fit
    ; # 
    global_label_exp = missing,# array of the experiment used to save data, Now optional if not fixed should be assigned to ["data_1","data_2",...]
    integrator=Tsit5(), 

    )
    
    # input integrity 
    validate_kinbiont_inputs(Kinbiont_data,Kinbiont_model)

 
    # Step 3 preprocesse
    Kinbiont_data = Kinbiont_data_preprocessing(Kinbiont_data,Kinbiont_options_preprocessing)
      # Step 3.1 multiple scattering correction


      # Step 3.2 average replicates
      # Step 3.3 blank computation 
      # Step 3.4 blank subtraction
      # Step 3.5 negative value removal

    # Step 4 fit


end



function validate_kinbiont_inputs(data::Kinbiont_data, models::Kinbiont_models_to_fit)

    # -------------------------
    # 1. Check Times length == number of columns in Data_matrix
    # -------------------------
    n_cols = size(data.Data_matrix, 2)
    if length(data.Times) != n_cols
        error("Mismatch: length(Times) = $(length(data.Times)) but Data_matrix has $n_cols columns.")
    end

    # -------------------------
    # 2. Check Labels length == number of rows in Data_matrix
    # -------------------------
    n_rows = size(data.Data_matrix, 1)
    if length(data.Labels) != n_rows
        error("Mismatch: length(Labels) = $(length(data.Labels)) but Data_matrix has $n_rows rows.")
    end

    # -------------------------
    # 3. Check model_list, CI, UB, LB have same length
    # -------------------------
    L = length(models.model_list)
    fields_to_check = [
        (:models_CI, models.models_CI),
        (:models_UB, models.models_UB),
        (:models_LB, models.models_LB),
    ]

    for (name, vec) in fields_to_check
        if vec !== missing && length(vec) != L
            error("Mismatch: length($name) = $(length(vec)) but model_list has length $L.")
        end
    end

    # -------------------------
    # 4. Check each model has same number of parameters in CI, UB, LB
    # -------------------------
    for i in 1:L
        ci = models.models_CI[i]

        # UB and LB may be missing or a vector containing missing
        ub = (models.models_UB === missing || models.models_UB[i] == [missing]) ? nothing : models.models_UB[i]
        lb = (models.models_LB === missing || models.models_LB[i] == [missing]) ? nothing : models.models_LB[i]

        p = length(ci)
        pd = length(lb)

        # Check UB
        if ub !== nothing && length(ub) != p
            error("Model $i: UB length $(length(ub)) does not match CI length $p.")
        end

        # Check LB
        if lb !== nothing && length(lb) != p
            error("Model $i: LB length $(length(lb)) does not match CI length $p.")
        end


          if ub !== nothing && length(ub) != pd
            error("Model $i: UB length $(length(ub)) does not match LB length $pd.")
        end
    end

    
end

function generic_model_selector(model_function)

    if typeof(model_function) == String
    
        model_string = Kinbiont_models[model_function].name
        model_function = Kinbiont_models[model_string].func
    
    
    else
    
        model_string = "custom"
    
    
    end  
    return model_function, model_string
    end
    
    