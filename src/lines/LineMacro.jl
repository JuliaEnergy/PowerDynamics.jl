using MacroTools

macro Line(typedef, prep, func_body)
    return create(typedef, prep, func_body)
end

macro Line(typedef, func_body)
    return create(typedef, nothing, func_body)
end

function create_line(name, prep, parameters, func_body)
    struct_exp = create_struct(name, parameters)

    ce_call = :(construct_edge(par::$(name)))
    extracted_parameters = map(sym -> :( $sym = par.$sym ), parameters)
    cl_body = quote end
    append!(cl_body.args, extracted_parameters)
    if (prep !== nothing)
        append!(cl_body.args, prep.args)
    end
    cl_function = Expr(:function, ce_call, cl_body)

    rhscall = :(rhs!(e,v_s,v_d,p,t))
    rhsbody = quote end
    rhsbody.args[1] = func_body.args[1]
    append!(rhsbody.args, [:(source_voltage = v_s[1] + v_s[2]*im)])
    append!(rhsbody.args, [:(destination_voltage = v_d[1] + v_d[2]*im)])
    append!(rhsbody.args, func_body.args)

    es_real = [:(e[1] = real(current_vector[1]))]
    es_imag = [:(e[2] = imag(current_vector[1]))]
    ed_real = [:(e[3] = real(current_vector[2]))]
    ed_imag = [:(e[4] = imag(current_vector[2]))]

    append!(rhsbody.args, [es_real; es_imag; ed_real; ed_imag])

    rhs_function_exp = Expr(:function, rhscall, rhsbody)
    edge_exp = :(return StaticEdge(f! = rhs!, dim = 4))
    append!(cl_function.args[2].args, [rhs_function_exp, edge_exp])

    ret = quote
        $(struct_exp)
        $(cl_function)
    end
    return ret
end

function create_struct(name, parameters)
    struct_def = Expr(
        :struct, false,
        :($name <: AbstractLine),
        Expr(:block, parameters..., # set all the parmeters as fields in the struct
            Expr(:(=), # define the constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...)
            )
        )
    )
end

function create_showdefinition(exp, name)
    mainexstr = "$(copy(exp)|>rmlines|> MacroTools.striplines)"
    return :(showdefinition(io::IO, ::Type{$name}) = println(io, $mainexstr))
end

function create(typedef, prep, func_body)
    @capture(typedef, name_(parameters__))
    mainex = create_line(name, prep, parameters, func_body)
    showex = create_showdefinition(mainex, name)
    append!(mainex.args, [showex])
    return esc(mainex)
end
