import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(main_model_part, custom_settings):
    return StructuralMechanicsExtract(main_model_part, custom_settings)

#class StructuralMechanicsExtract(StructuralMechanicsAnalysis):
class StructuralMechanicsExtract(MechanicalSolver):
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        km.Logger.PrintInfo("::[CustomScipyBaseSolver]:: ",
                            "Construction finished")

    def GetDefaultParameters(cls):
        this_defaults = km.Parameters("""{"scheme_type"        : "dynamic"}""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def SolveSolutionStep(self):
        self.M = self._MatrixComputation(mat="mass")
        self.K = self._MatrixComputation(mat="stiff")

    def _create_solution_scheme(self):
        """Create the scheme for the scipy solver.
        The scheme determines the mass and stiffness matrices
        """
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "dynamic":
            solution_scheme = ksm.EigensolverDynamicScheme()
        else:  # here e.g. a stability scheme could be added
            err_msg = ("The requested scheme type \"" + scheme_type +
                       "\" is not available!\n")
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)
        return solution_scheme

    def _create_linear_solver(self):
        ''' Linear solver will not be used. But eventually the solution strategy calls the solver's clear function.
        To avoid crashing linear solver is provided here'''
        return km.LinearSolver()

    def _create_mechanical_solution_strategy(self):
        if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            warn_msg = ("In case an eigenvalue problem is computed an" +
                        "elimantion builder shall be used to ensure boundary" +
                        " conditions are applied correctly!")
            km.Logger.PrintWarning("CustomScipyBaseSolver", warn_msg)

#    def _create_mechanical_solution_strategy(self):
#        if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
#            warn_msg = ("In case an eigenvalue problem is computed an" +
#                        "elimantion builder shall be used to ensure boundary" +
#                        " conditions are applied correctly!")
#            km.Logger.PrintWarning("CustomScipyBaseSolver", warn_msg)

        eigen_scheme = self.get_solution_scheme()  # The scheme defines the matrices
        computing_model_part = self.GetComputingModelPart()
        builderSolver = self.get_builder_and_solver()
        linear_solver = self.get_linear_solver()
        return km.ResidualBasedLinearStrategy(computing_model_part,
                                              eigen_scheme, linear_solver,
                                              builderSolver, False, False,
                                              False, False)

    def _MatrixComputation(self, mat="stiff"):  # add option for constrianed or unconstrained
        space = km.UblasSparseSpace()
        if mat == "mass":
            self.GetComputingModelPart().ProcessInfo.SetValue(ksm.BUILD_LEVEL,
                                                              1)
        elif mat == "stiff":
            self.GetComputingModelPart().ProcessInfo.SetValue(ksm.BUILD_LEVEL,
                                                              2)
        scheme = self.get_mechanical_solution_strategy().GetScheme()
        aux = self.get_mechanical_solution_strategy().GetSystemMatrix()
        space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = space.CreateEmptyVectorPointer()
        space.ResizeVector(b, space.Size1(aux))
        space.SetToZeroVector(b)
        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector(xD, space.Size1(aux))
        space.SetToZeroVector(xD)

        # Build matrix
        builderSolver = self.get_builder_and_solver()
        builderSolver.Build(scheme, self.GetComputingModelPart(), aux, b)

        # Apply Constraints
        builderSolver.ApplyConstraints(scheme, self.GetComputingModelPart(),
                                       aux, b)
        # Apply Boundary Conditions
        builderSolver.ApplyDirichletConditions(scheme,
                                               self.GetComputingModelPart(),
                                               aux, xD, b)
        # Return matrix converted to scipy
        return km.scipy_conversion_tools.to_csr(aux)


with open("ProjectParameters.json", 'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())
model = km.Model()
simulation = StructuralMechanicsExtract(model, parameters)
simulation.Initialize()
simulation.SolveSolutionStep()
simulation.Finalize()
simulation.Run()
print(simulation.M)
