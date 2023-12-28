
# Welcome
``````
.
├── functions.py
├── input_data
│   └── neff_Voltage_sweep_TE.csv
├── optical_frequency_response.py
├── quality_factor_sweep.py
├── README.md
├── requirements.txt
└── user_parameters
    └── user_inputs.py

``````
## Why this Repository?
This archive contains simple Python scripts used to analytically calculate the properties of a PN or PIN ring modulator.

## Quick Calculation
1. Install the requirements.txt
2. After downloading the repository, navigate to the `user_parameters/user_inputs.py` directory and edit the `user_input.py` file.
3. You can refer to the `function.py` file to explore the available functions used to calculate the properties of the ring modulator.
4. Execute the `optical_frequency_response.py` example to get the frequency response of the add-drop ring modulator.

### References
1. Lukas Chrostowski, and Michael Hochberg. **Silicon Photonics Design**, 1st. Cambridge University Press, 2015, 217-258.
2. W. Bogaerts et al., "**Silicon microring resonators**," *Laser & Photon. Rev.*, vol. 6, pp. 47-73, 2012.
