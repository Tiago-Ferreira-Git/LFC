
# 🌐 Load Frequency Control

[comment]: [![N|Solid](https://cldup.com/dTxpPi9lDf.thumb.png)](https://nodesource.com/products/nsolid)

Welcome to the **Load Frequency Control** repository! This project aims to provide a comprehensive toolkit for simulating and implementing load frequency control (LFC) strategies. Designed for researchers and students, this repository offers tools to model and analyze the dynamic behavior of power systems under varying load and Renewable Energy sources (RES) conditions, ensuring stable and efficient operation.


## ⚡ Key Features

- Dynamic Simulation: Simulate the behavior of power systems under different operational scenarios, including generation, transmission, and distribution.
- Load Frequency Control (LFC): Implement and test various LFC strategies to maintain system frequency within desired limits, ensuring stability and reliability.
- Real data: The nominal load and RESs profiles were taken from the ENTSO-E database to simulate a real evolution of these parameters, applicable to any MATPOWER testcase.
- Visualization Tools: Generate clear and informative plots to visualize system performance, frequency deviations, and control actions.
 
## 🧱 Dependencies

To run the code as is, this repository uses a control algorithm that is available in the toolbox: https://decenter2021.github.io/download/ . Moreover, the load and RESs profiles were taken from the ENTSO-E database that is available here: https://transparency.entsoe.eu/ .

## 🌳 Repository Tree
```bash
📦area partition
 📦LFC
 ┣ 📂data
 ┃ ┣ 📂Load and Renewables Generation Data
 ┃ ┃ ┣ 📂Load # ENTSO-E load data
 ┃ ┃ ┃ ┣ 📂2022 # 2022 data
 ┃ ┃ ┃ ┣ 📂2024 # 2024 data
 ┃ ┃ ┃ ┗ 📜processing.m # Processes ENTSO-E load data
 ┃ ┃ ┗ 📂RES
 ┃ ┃ ┃ ┣ 📂2022 # 2022 ENTSO-E RESs data
 ┃ ┃ ┃ ┣ 📂IEEE118 # Data available for the IEEE118 tescase 
 ┃ ┃ ┃ ┃ ┣ 📂Forecast
 ┃ ┃ ┃ ┃ ┣ 📂Measured
 ┃ ┃ ┃ ┃ ┣ 📜IEEE118.m # Processes this data
 ┃ ┃ ┃ ┣ 📜processing.m # Processes ENTSO-E RESs data
 ┃ ┣ 📂results # Folder to keep the results
 ┃ ┣ 📜IEEE14.mat # IEEE 14 partition example 
 ┃ ┣ 📜IEEE14_res.mat # IEEE 14 partition example with RESs
 ┃ ┣ 📜IEEE118.mat # IEEE 14 partition example 
 ┃ ┣ 📜load_profile.mat # Result of load processing.m
 ┃ ┣ 📜permute_matrix.m
 ┃ ┣ 📜ren_profile.mat # Result of load processing.m
 ┃ ┣ 📜res_profile.mat # Result of res processing.m
 ┃ ┣ 📜res_profile_118.mat # Result of IEEE18.m
 ┣ 📂fig # Folder to keep the figures
 ┣ 📜.gitignore
 ┣ 📜area.m
 ┣ 📜area_partitioning.m
 ┣ 📜compute_tielines.m
 ┣ 📜discrete_dynamics.m
 ┣ 📜get_disturbance_profile.m
 ┣ 📜get_g.m
 ┣ 📜get_global_ss.m
 ┣ 📜initial_conditions.m
 ┣ 📜nonlinear_model.m
 ┣ 📜permute_matrix.m
 ┣ 📜plot_metrics.m
 ┣ 📜README.md
 ┣ 📜testing_nonlinear.m
 ┣ 📜t_settling.m
 ┗ 📜update_dynamics.m
```
## 🧪 Results

<p align="center">
<img src="https://raw.githubusercontent.com/andresz1/size-limit-action/master/assets/pr.png"
  alt="Size Limit comment in pull request about bundle size changes"
  width="686" height="289">
</p>

## 🚧 Future work


## 🗨️ Contacts

## License
This project is licensed under the MIT License. See the LICENSE file for more details.


