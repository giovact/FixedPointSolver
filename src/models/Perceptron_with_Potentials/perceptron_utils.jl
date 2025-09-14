generalization_error(Model::TM) where TM <: Perceptron = acos(Model.O[:R]) / π


function setup_Perceptron(modeltype::String,ruletype::String;
                            β::FT = 0.0,
                            α::FT = 0.0,
                            γ::FT = 0.0,
                            y::Int = 100000,
                            θ::FT=0.0,
                            R::FT = 1e-2,
                            Q::FT = 1e-2,
                            P::FT = 1e-2,
                            Q0::FT = 1e-2,
                            δ::FT = 1e-2,
                            Xint::TI = NoInt()) where TI <: IntegrationMethod

    rule = set_rule(ruletype)
    if modeltype == "BinaryPerceptronRS"
        return BinaryPerceptronRS(β,α,rule,[R;Q])
    elseif modeltype == "BinaryPerceptron1RSB"
        @assert θ!=0.0
        @assert Q0> Q 
        return BinaryPerceptron1RSB(β,α,θ,rule,[R;Q;Q0])
    elseif modeltype == "CoupledBinaryPerceptrons"
        @assert y<=10
        @assert P >Q
        return CoupledBinaryPerceptrons(β,α,γ,y,rule,[R;Q;P])
    elseif modeltype == "CoupledBinaryPerceptronsNoTrace"
        @assert y<=10
        @assert P >Q
        return CoupledBinaryPerceptronsNoTrace(β,α,γ,y,rule,[R;Q;P])
    elseif modeltype == "BinaryPerceptronYInfNoDelta"
        return BinaryPerceptronYInfNoDelta(β,α,γ,rule,[R;Q])
    elseif modeltype == "BinaryPerceptronYInf"
        return BinaryPerceptronYInf(β,α,γ,rule,[R;Q;δ],Xint)
    elseif modeltype == "BinaryPerceptron1RSBDyntheta1"
        return BinaryPerceptron1RSBDyntheta1(β,α,rule,[0.0])
    end
end

_is_Teacher(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:R] >= 0.98
_is_reallyTeacher(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:R] >= 0.99
_is_NOT_Teacher(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:R] < 0.98
_is_NOT_TeacherReally(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:R] < 0.999

function _is_Teacher_NOT_Dominant(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod}
    Model_teacher = create_model(deepcopy(Model), Model.params)
    Model_teacher.O[:R] = 0.9999
    Model_teacher.O[:Q] = 0.9999
    if :P in names(Model.O,1)
        Model_teacher.O[:P] = 0.9995
    end
    
    solve_SPequations!(Model_teacher, X; tol = 1e-7, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 3000)
    if abs(Model.aux["f"] - Model_teacher.aux["f"]) < 1e-4
        return true
    end
    return Model.aux["f"] < Model_teacher.aux["f"] && ( Model_teacher.O[:R] - Model.O[:R] ) > 0.05
end

function _is_Teacher_Dominant(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod}
    Model_NOTteacher = create_model(deepcopy(Model), Model.params)
    Model_NOTteacher.O[:R] = 1e-1
    Model_NOTteacher.O[:Q] = 1e-1
    if :P in names(Model.O,1)
        Model_teacher.O[:P] = 2e-1
    end
    solve_SPequations!(Model_NOTteacher, X; tol = 1e-7, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 3000)

    return Model.aux["f"] < Model_NOTteacher.aux["f"] && ( Model.O[:R] - Model_NOTteacher.O[:R]) > 0.04
end

_is_Entropy_Positive(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.aux["entropy"] >0
_is_Entropy_Negative(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.aux["entropy"] <0
_is_Free_Energy_Negative(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.aux["f"] <0 
_isRSstable(Model::Perc,X::TI,Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.aux["dAT"] <0 
_is_complexity_zeroem4(Model::Perc, X::TI, Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.aux["complexity"] < 1e-4 

_students_correlated_notwithTeacher(Model::Perc, X::TI, Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:R] < 0.98 && Model.O[:Q] > 0.99

_above_TDyn(Model::Perc, X::TI, Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:Q0] - Model.aux["Q1"] < 0.01
_Q0notfrozen(Model::Perc, X::TI, Oold::NamedVec,tol::FT) where {Perc<:Perceptron,TI <: IntegrationMethod} = Model.O[:Q0]  < 0.98


# this should be deprecated
dict_functions = Dict("is_Teacher"=>_is_Teacher,
                    "is_reallyTeacher"=>_is_reallyTeacher,
                    "is_NOT_Teacher"=>_is_NOT_Teacher,
                    "is_Teacher_NOT_Dominant"=>_is_Teacher_NOT_Dominant,
                    "is_Teacher_Dominant"=>_is_Teacher_Dominant,
                    "is_NOT_TeacherReally"=>_is_NOT_TeacherReally,
                    "is_Entropy_Positive"=>_is_Entropy_Positive,
                    "is_Entropy_Negative"=>_is_Entropy_Negative,
                    "is_Free_Energy_Negative"=>_is_Free_Energy_Negative,
                    "isRSstable"=>_isRSstable,
                    "is_complexity_zeroem4"=>_is_complexity_zeroem4,
                    "students_correlated_notwithTeacher"=>_students_correlated_notwithTeacher,
                    "above_TDyn"=>_above_TDyn,
                    "Q0notfrozen"=>_Q0notfrozen
)

