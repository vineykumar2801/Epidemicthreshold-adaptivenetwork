## Overview
The content within this repository comprises the necessary code to replicate the figures outlined in the article titled "Estimating epidemic threshold under individual vaccination behavior and adaptive social connections: A game-theoretic complex network model" by Viney Kumar and Samit Bhattacharyya.

## Time Series Simulation Code
This section encompasses various .py files to generate the time series data corresponding to each result presented in the aforementioned article.

**Graph_visualization.py**: This script facilitates the assignment of all possible states to the graph, including Susceptible (S), Infection (I), Recovered (R), and Vaccinated (V).

**nbd_index.py**: This file generates the index numbers of all neighbors for each node within the network.

**fermi_delta_payoff.py**: Utilizing the model framework outlined, this script computes the payoff for vaccinators and non-vaccinators and the payoff differences.

**av_vaccine_risk.py**: Defined to update the vaccination state and perceived vaccine risk for all nodes within the virtual network.

**Removing_old_connection.py**: This script removes old connections within the virtual network.

**Adding_new_connection.py**: Designed to add new virtual connections to the network.

**Time_series_code.py**: This script generates the time series data for all states outlined within our model framework.
