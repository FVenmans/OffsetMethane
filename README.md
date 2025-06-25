# OffsetMethane
This is the code and the supplementary material complementing the article "Temporary carbon removals to offset methane emissions."
The Matlab file "IRFmain.m" makes all the figures as well as Table 1. It also calibrates FAIR2 on historic emissions. It uses the functions IRF_FAIR2 and FAIR2. 
The function "IRF_FAIR2.m" makes an impact response function for a pulse of any greenhouse gas.  
The function "FAIR2.m" takes in any emission path (for different greenhouse gases and exogenous forcing) and creates several outputs such as temperature (atmosphere and ocean at different depths), forcing, gas concentrations (in different 'reservoirs' or boxes for carbon). 
The excel file "Value of an offset with variable carbon storage" allows you to introduce the tonnes of carbon stored in each year on sheet 2 column M to obtain the value of an impermanent carbon removal in cell N4. Lenght of the project can also be set on sheet 1 B3 and failure rate on sheet 1 B5.  
