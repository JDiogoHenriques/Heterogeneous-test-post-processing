# Heterogeneous-test-post-processing
Python code for doing the post-processing of numerical heterogeneous tests for solid elements.


## Features
- Extracts the results directly from ABAQUS "*.odb files".
- Plots principal stress and strain diagrams.
- Generates MatchID synthetic image generation input files from ABAQUS.
- Generates synthetic images and applies DIC to obtain virtual experiment results.
- Calculates the results of two heterogeneous criteria: (1) IT1 (2) Rotation angle and plots the results.
- Plots the yield surface on the normalised stress space with the material points.

## Usage
1. Update the test name and the user-specific options on the "PythonScript_main.py" script.
2. Run the "PythonScript_main.py" script.
