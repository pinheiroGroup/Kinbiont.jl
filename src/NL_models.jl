# from "Statistical evaluation of mathematical models for microbial growth"
# non linear model 
# exponential



# Logistic

function NL_model_logistic(p,times)

     u = p[1] ./ (1 .+ exp.( .- p[2] .* ( times .- p[3] )))
   
    return u

end 



function NL_model_Gompertz(p,times)

     u =  p[1] .* exp.( - exp.( - p[2] .* ( times .- p[3] )))

    return u

end 


function NL_model_Bertalanffy(p,times)

    u =  p[1] .+ (p[2] .- p[1]) .* (1 .- exp.( .- p[3] .* times) ).^(1 ./p[4])

    return u

end 
#Richards

function NL_model_Richards(p,times)

     u = p[1] ./ (1 .+ p[2] .* exp.( .- p[3] .* ( times .- p[4] ))).^(1 ./ p[5])

    return u

end 
#Morgan


function NL_model_Morgan(p,times)

    u= (p[1] .* p[2].^ p[3] .+ p[4] .* times .^ p[3]) ./( p[2] .^p[3] .+ times .^ p[3])


    return u

end


#Weibull


function NL_model_Weibull(p,times)

    u = p[1] .- (p[1] .-p[2]) .* exp.( - (p[3] .* times ) .^p[4]) 


    return u

end



