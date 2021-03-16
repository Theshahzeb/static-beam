#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 10:16:44 2021

@author: shahzeb
"""
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Import scipy modules
import KratosMultiphysics.scipy_conversion_tools

with open("ProjectParameters.json", 'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters(parameter_file.read())
model = KratosMultiphysics.Model()
# simulation = CustomScipyBaseSolver(model, parameters)
simulation = StructuralMechanicsAnalysis(model, parameters)        
simulation.Run()
print(simulation)
print("done")
print("-"*50)