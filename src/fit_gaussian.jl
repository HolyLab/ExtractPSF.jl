#Fit gaussians according to the equation: y = e^-(x^2/2std^2)

#1D gaussian centered at zero, peak height of 1
eval_gaussian(x, std) = exp(-((x^2)/(2*std^2)))

function gaussian_obj(std, shift, scale, bias, data::OffsetArray{T,1}) where {T}
    if std > 0 && scale > 0
        p = zero(std)
        for i in first(axes(data))
            #squared error of fitted gaussian and data, evaluated at each index of the data array
            p += (data[i]-(bias+scale*eval_gaussian(i+shift, std)))^2
        end
        return p
    else #infeasible std, scale
        m = min(std, scale)
        return sum(data.^2)*(1+m^2)
    end
end

#std, shift, and scale packed into params array
gaussdecode(m) = (m[1], m[2], m[3])

gaussian_obj(params, bias, data) = gaussian_obj(gaussdecode(params)..., bias, data)

#the fitted shift parameter is in the same coordinate system as the data array
#Therefore if you want the returned shift to reflect the distance from the center of the data vector,
#data[0] should correspond with the center of the array
function fit_gaussian_1d(data::OffsetArray{T,1}, bias; guess_std=1.0, guess_scale=1.0) where {T}
    gobj = x->gaussian_obj(x, bias, data)
    firstguess = [guess_std; 0; guess_scale]
    result = optimize(gobj, firstguess, LBFGS(); autodiff=:forward)
    m = Optim.minimizer(result)
    fval = Optim.minimum(result)
    _std, shift, scale = gaussdecode(m)
    return _std, shift, scale, fval
end
