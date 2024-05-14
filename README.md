# Axions

1. Setup axionCAMB from https://github.com/dgrin1/axionCAMB and use the instructions in the AxionCAMB_documentation folder.

2. Setup axionHMcode from https://github.com/SophieMLV/axionHMcode. Modify the file path in the example_file.py to use axionCAMB.

3. example_file.py (in this repository), pk_generate.py, plot.py are examples of varying the paramenters in input_file.py, generating non-linear power spectra with axionHMcode, and plotting the power spectra.

4. To get axionCAMB also output the Hubble parameters, replace the axion_background.f90 file in axionCAMB with the one on this page to print out Hubble parameters. It is better to start from a new copy of axionCAMB so that it does not get confused with the one used in axionHMcode. The generated files are aosc_param.dat and littleh_test.dat. Both files are read in littleh_test.py. aosc_param.dat allows you to "manually" calculate the Hubble parameter, littleh_test.py directly outputs the Hubble parameter.

5. new_finite_der.py gives an example of running axionCAMB in a python script. It also sets up the calculation of finite derivatives.

6. Calculate finite derivatives with new_diffder_cal.py, and plot them in Jupyter Notebooks that contain "compder". 

7. Alternatively, Pipeline_clean.ipynb calculates the derivatives from polynomial fitting. Best to try Pipeline.ipynb and Fisher_Analysis.ipynb in https://github.com/etrott12/Axion first on a computer with sufficient memory first.

