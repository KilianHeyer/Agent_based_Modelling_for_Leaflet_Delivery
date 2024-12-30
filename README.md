# Masterthesis_Kilian_Heyer
Agent based Modelling for Leaflet Delivery 
This repository contains the code, data, and results for simulating the distribution of advertising material (e.g., brochures). 
The simulation estimates the delivery time using a Nearest Neighbor Algorithm and other models, allowing for comparison of various delivery areas. 
The project was created as part of the Master's thesis by Kilian Heyer.

/
|-- /code/                     # Contains GAMA models and related code
|   |-- Concorde_export.gaml    
|   |-- Following_TSP_Solver_address.gaml 
|   |-- Model_1.gaml             
|   |-- Model_2_for_all_parking_spots.gaml
|   |-- Model_3_for_random_parking.gaml
|
|-- /data/                     # Input data for simulations
|   |-- Data_for_area_502010/  # Data for area 502010
|   |-- Data_for_area_502021/  # Data for area 502021
|   |-- Data_for_area_502023/  # Data for area 502023
|   |-- Data_for_area_502601/  # Data for area 502601
|   |-- fake_timetable.csv     # delivery timetable with fake times
|
|-- /results/                  # Simulation results
|   |-- 502010_model_2_results.xlsx
|   |-- 502011_model_3_results.xlsx
|   |-- 502023_model_3_results.xlsx
|   |-- 502601_model_3_results.xlsx
|
|-- /tsp_files_and_concorde/   # Files related to TSP solver
|   |-- Data_for_Concorde_area_502010_with_walking_graph.tsp
|   |-- solution_for_area_502010_on_walking_graph
|   |-- solution_for_area_502010_on_walking_graph_adjustet to GAMA
|
|-- README.md                  # Documentation (this file)



1. Code

    Location: /code/
    Files:
        Concorde_export.gaml: Code for exporting TSP-related data to the Concorde solver. This model was created during and for the Masther Thesis of Kilian Heyer it is a working file and not designed to be                                 used outside the context of the research for the thesis and lacks proper documentation.
        Following_TSP_Solver_address.gaml: Follows the TSP solver's results for route optimization. This model was created during and for the Masther Thesis of Kilian Heyer it is a working file and not                                                   designed to be used outside the context of the research for thesis and lacks proper documentation.
        Model_1.gaml, Model_2_for_all_parking_spots.gaml, Model_3_for_random_parking.gaml: Main GAMA simulation models with varying approaches for parking spot and delivery optimization.

2. Data

    Location: /data/
    Contents:
        Subdirectories and zip files for different areas (e.g., Data_for_area_502010).
        A sample fake_timetable.csv for delivery times. The fake_timetable is a dummy data containing fake data. Because of a data protection agreement with Feibra GmbH the origianl Data must not be                   published. A exact reproduktion of the thesis results is not possible with this time table but proportionate result should be generated.


3. Results

    Location: /results/
    Contents:
        Excel files with simulation results which were used for all statistics.

4. TSP Files

    Location: /tsp_files_and_concorde/
    Contents:
        Data for the Concorde TSP solver (e.g., Data_for_Concorde_area_502010_with_walking_graph.tsp).
        Solutions provided by the TSP solver (e.g., solution_for_area_502010_on_walking_graph).
        Solutions adapted to be useable by the Model Following_TSP_Solver_address.gaml (e.g., solution_for_area_502010_on_walking_graph_adjustet to GAMA)
   
1. Prerequisites

    GAMA Platform: To run the model GAMA version 1.8.2 is needed. It may funktion with higher versionst but this was not testet. For further Information and download go to: https://gama-platform.org/download
    Concorde TSP Solver: Required for generating optimal solutions for the TSP Problem. Ensure it's installed and configured if needed. 
    For further Information and download go to: https://www.math.uwaterloo.ca/tsp/concorde/index.html

2. Running the Models

    Open GAMA: Open the .gaml files in the /code/ directory in GAMA.
    Load Data: Adjust file paths in the model to point to the appropriate data in the /data/ folder. 
               To run the models sucessfully they must be given: the adresses, streets, start_point,road for deliverer as shp and the delivery_time as csv. 
               Also A csv must be stated for the results to be saved.
               An error will appear if the file to save the data does not exist.  
    Run the Simulation: Execute the model and analyze results in GAMA or from the /results/ folder.

3. Using the TSP Solver

    Export data for the TSP solver using Concorde_export.gaml.
    Solve the TSP problem with the Concorde solver.
    Use the results in the simulation by changeing "," to "." and deleting the first(0) and the highest number. Then copy them into Following_TSP_Solver_address.gaml as optimal_path in line 82.

If you have any questions or need support, please reach out to Kilian Heyer at kilian.heyer@gmail.com.
