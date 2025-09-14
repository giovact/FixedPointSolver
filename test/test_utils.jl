function test_model(Model::TM,X::TI) where {TM<: FPModel, TI <: IntegrationMethod}
    conv, iter, vareps = solve_SPequations!(Model, X; tol = 1e-5, niter = 10000,showinfo = false,printprogress = false)
	return conv, iter
end

perceptron_types=[BinaryPerceptronRS;
                    BinaryPerceptron1RSB;
                    CoupledBinaryPerceptrons;
                    CoupledBinaryPerceptronsNoTrace
                    ]


models_retrieval = [RBMBinaryReLU;
                    RBMBinaryReLUT0;
                    RBMBinaryBinaryLL;
                    HopfieldLL;
                    HopfieldHL;
                    HopfieldHLT0;
                    BAMLL;
                    BAMHLT0Mixtures]


parametric_models = vcat(perceptron_types, models_retrieval)
