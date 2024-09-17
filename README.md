## Overview
The content within this repository comprises the necessary code to replicate the figures outlined in the article titled "Estimating the epidemic threshold under individual vaccination behaviour and adaptive social connections: A game-theoretic complex network model" by Viney Kumar, Chris T. Bauch, and Samit Bhattacharyya.

## Time Series Simulation Code
This section encompasses various .py files to generate the time series data corresponding to each result presented in the aforementioned article.

**Graph_visualization.py**: This script facilitates the assignment of all possible states to the graph, including Susceptible (S), Infection (I), Recovered (R), and Vaccinated (V).

**nbd_index.py**: This file generates the index numbers of all neighbors for each node within the network.

**fermi_delta_payoff.py**: Utilizing the model framework outlined, this script computes the payoff for vaccinators and non-vaccinators and the payoff differences.

**av_vaccine_risk.py**: Defined to update the vaccination state and perceived vaccine risk for all nodes within the virtual network.

**Removing_old_connection.py**: This script removes old connections within the virtual network.

**Adding_new_connection.py**: Designed to add new virtual connections to the network.

**Time_series_code.py**: This script generates the time series data for all states outlined within our model framework. This code will provide us .csv file of all the states.

## Epidemic Threshold Simulation
Within this section, you will find a collection of .py files dedicated to computing the epidemic threshold associated with each state outlined in the article.

**epithres_LI_fermi_nbd.py**, **epithres_MI_fermi_nbd.py**, and **epithres_HI_fermi_nbd.py**: These files generate the necessary probabilities for each state, as elaborated in the "theoretical analysis of epidemic threshold" section, catering to LI, MI, and HI states, respectively.

**epidemic_LI_adaptiness.py**, **epidemic_MI_adaptiness.py**, and **epidemic_HI_adaptiness.py**: Responsible for calculating the final epidemic threshold value for LI, MI, and HI states, respectively.



## Figures Plot
We have simulated Matlab code to generate plots for all the .csv files obtained from executing the **Time_series_code.py** script.

## Contact information
**Viney Kumar:** [vk981@snu.edu.in](mailto:vk981@snu.edu.in)
**Chris T. Bauch:** [cbauch@uwaterloo.ca](mailto:cbauch@uwaterloo.ca)
**Samit Bhattacharyya:** [samit.b@snu.edu.in](mailto:samit.b@snu.edu.in)
