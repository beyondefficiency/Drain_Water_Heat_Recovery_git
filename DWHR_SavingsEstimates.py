# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:13:21 2017

This script predicts the energy savings of DWHR devices. It works based on the following principles:
    
1. It uses glob and a for loop to iterate through and repeat all calculations on a list of draw profiles.This means that the analysis can be performed on 
several different draw profiles all at once by storing the draw profiles in the specified folder. The draw profiles are specified in the Path_DrawProfiles 
variable (Line 56 at this time)

2. It calculates the effectiveness using the methods documented at https://title24stakeholders.com/measures/cycle-2019/drain-water-heat-recovery/
A brief summary: Curve coefficients define a performance map predicting a a correction factor stating the performance of a generic DWHR device as
a function of the drain-side and cold-side water flow rates. This correction factor states the effectiveness of the device as a percentage (In decimal
form) of the performance at 2.51 gal/min. For instance, a correction factor of 1.1 for a DWHR device with rated effectiveness of 45.6% (at 2.51 gal/min)
under a certain set of conditions means that the effectiveness under those conditions is 1.1 * 45.6% = 50.16%

3. It uses these methods to calculate the effectiveness for each draw in each draw profile (During the for loop) assuming it's installed in the equal
flow configuration. This is simply using the performance map to find the correction factor, multiplying by the rated effectiveness to find the effectiveness
in that draw, and calculating Q = Effectiveness * m * C_p * dT

4. It performs the same calculations assuming that the DWHR device is installed in the Unequal - Water Heater configuration. The main difference is that the 
flow rate through the cold side is now less than the flow rate through the drain side, and must be calculated. The function Calculate_Fraction_Cold_ThroughDWHR
performs an energy balance to identify the fraction of water passing through the fixture that came from the cold side, based on the outlet temperature from the 
cold side of the DWHR device. That fraction is then multiplied by the total flow rate to identify the cold-side flow rate. The same process as in step 3 can then 
be used

5. It calculates the performance of a DWHR device installed in the Unequal - Shower configuration. This is similar to the Unequal - Shower configuration in that
the cold-side flow rate will be less than the drain-side flow rate, and in that the cold side flow rate is a function of the cold-side outlet water temperature.
The difference is that, this time, the mixture of hot/cold water at the fixture is determined by the temperature of the water. This means that the outlet temperature
impacts the flow rate, which impacts the temperature, which impacts the flow rate, which... It's an iterative solution. This section sets up some starting assumptions,
then uses an while loop to implement the iterative solution. The while loop uses the flow rate through the cold side to search for conversion - When the assumed flow
rate before the calculations and the flow rate after the calculations differ by less than 0.01 gal/min, the solution is accepted. It's imperative that this code 
converge quickly, to reduce computation time. To assist this, the correction factor method was replaced with another method from Ramin Manouchehri because it was
found to converge in fewer iterations. That work is documented at: https://uwspace.uwaterloo.ca/handle/10012/10035

6. It then saves the file with the new calculations to a file with '_Analyzed' added to the end of the file name.

Desired future changes:
    -It would be nice to program a version of this script as a function. This could be combined with the draw profile generation script (Which would also be programmed
        as a function). Doing so would allow the user to set up a list of draw profiles they want to analyse, then iterate through the list to both create the profile
        and run it through these calculations to predict the performance.

@author: pgrant
"""

#%%-------------------IMPORT STATEMENTS---------------

import pandas as pd
import numpy as np
import itertools
import math
import glob

#%%------------------INPUTS------------------------
    
Path_DrawProfiles = r'C:\Users\Peter Grant\DWHR\Profiles' #The folder where all of the draw profiles to be analyzed are located

Coefficients_Generic_Vertical_Unequal = np.fromfile(r'C:\Users\Peter Grant\Dropbox (Beyond Efficiency)\Beyond Efficiency Team Folder\Frontier Energy-TRC 2022 Title 24 MF CASE Support\DWHR\Analysis\Coefficients\Generic_Vertical_Unequal.csv') #The path to the coefficients for the unequal flow generic model of DWHR devices
Coefficients_Generic_Vertical_Equal = np.fromfile(r'C:\Users\Peter Grant\Dropbox (Beyond Efficiency)\Beyond Efficiency Team Folder\Frontier Energy-TRC 2022 Title 24 MF CASE Support\DWHR\Analysis\Coefficients\Generic_Vertical_Equal.csv') #The path to the coefficients for the equal flow generic model of DWHR devices

Effectiveness_Rated = 0.42 #The rated effectiveness of the DWHR device tested under CSA B55.1 or IAPMO IGC 346-2017

#%%------------------CONSTANTS------------------------

SpecificHeat_Water = 0.998 #Btu/(lb_m-F) @ 80 deg F, http://www.engineeringtoolbox.com/water-properties-d_1508.html
Density_Water = 8.3176 #lb-m/gal @ 80 deg F, http://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
    
#Hot water temperature constants are taken from pg B-3 of the 2016 CBECC ACM reference manual
#These constants can be changed if wanting to try different arrangements (E.g. A different water heater set temperature)
Temperature_Shower = 105 #deg F
Temperature_WaterHeater = 115 #deg F

#The temperature of water entering the DWHR device will be lower than the shower temperature
Temperature_Drain_Inlet = 100.4 #deg F, the temperature of water entering the drain side of the device
    
#%%-------------------DEFINE FUNCTIONS----------------

def FourthOrder(List, Condition): #Returns the value of a fourth order polynomial given the curve-fit parameters, and the value of the independent variable
    return List[0] * Condition**4 + List[1] * Condition**3 + List[2] * Condition**2 + List[3] * Condition + List[4]
    
#This functions performs an energy balance at the fixture, using the specified flow rates, temperatures, and installation configuration to predict the
#fraction of total flow that is cold water, and therefore passes through the cold side of the DWHR device
def Calculate_Fraction_Cold_ThroughDWHR(Flow_Shower, Temperature_Mains, Temperature_Shower, Temperature_WaterHeater, Configuration):
    if Configuration == 'Equal': #If the DWHR device is installed in an equal flow configuration
        Flow_Drain = Flow_Shower #Then the drain side flow rate is equal to the total flow rate
        Flow_Cold = Flow_Shower #Then the cold side flow rate is equal to the total flow rate
        Fraction_Cold_ThroughDWHR = 1 #All of the water going down the drain passes through the cold side of the DWHR device
    elif Configuration == 'Unequal_WaterHeater' or Configuration == 'Unequal_Fixture': #If the DWHR device is installed in either the Unequal-WaterHeater or Unequal-Shower configuration
        Flow_Drain = Flow_Shower #Then the flow rate through the drain is equal to the shower flow rate
        Flow_Cold = Flow_Shower * (Temperature_Shower - Temperature_WaterHeater) / (Temperature_Mains - Temperature_WaterHeater) #Calculates the flow rate through the cold side of the DWHR device by performing an energy balance
        Fraction_Cold_ThroughDWHR = Flow_Cold / Flow_Drain #Calculates the fraction of water passing through the DWHR device that has passed through the cold side of the device
    return Fraction_Cold_ThroughDWHR

#This function evaluates a 2-dimensional polynomial with the inputs x, y and the coefficients m
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

#This makes it slightly easier to calculate the natural log of something in the code? Doesn't seem worth having as a function, may delete
def log(x):
    return math.log(x)

#This function finds the text between two different entries in a string. The inputs are he full string, the text preceeding the portion you want, and the text
#folowwing the portion you want. For instance, calling the function on this comment with First = 'For instance, ' and Last = ' on this comment' would return
#'calling the function'
def Find_Between(String, First, Last):
    try:
        Start = String.index(First) + len(First)
        End = String.index(Last, Start)
        return String[Start:End]
    except ValueError:
        return ""


#%%------------------CALCULATIONS--------------------

Draw_Profiles = glob.glob(Path_DrawProfiles +  '/*.csv') #Use glob to create a list of all of the .csv files in the Path_DrawProfiles folder

for i in Draw_Profiles: #Repeat once for each file in Draw_Profiles. Note that i will be the actual filename, not an index number
    Draw_Profile = pd.read_csv(i) #Open the .csv file 
    #Draw_Profile = Draw_Profile[Draw_Profile['Fixture'] == 'SHWR'] #Filter the draw profile to only include shower draws. This line will likely be removed before use

    Draw_Profile['Mixed Water Volume (gal)'] = Draw_Profile['Flow Rate (gpm)'] * Draw_Profile['Duration (min)'] #Calculate the total volume of each draw

    #Equal Flow calculations

    Draw_Profile['Effectiveness Generic Equal (-)'] = FourthOrder(Coefficients_Generic_Vertical_Equal, Draw_Profile['Flow Rate (gpm)'])
    Draw_Profile['Savings Generic Equal (Btu)'] = Draw_Profile['Effectiveness Generic Equal (-)'] * Draw_Profile['Mixed Water Volume (gal)'] * Density_Water * SpecificHeat_Water * (Temperature_Drain_Inlet - Draw_Profile['Mains Temperature (deg F)'])

    Savings_Generic_Equal = Draw_Profile['Savings Generic Equal (Btu)'].sum()/100000

    #Unequal-WaterHeater calculations
    

    Draw_Profile['Potable Flow Rate Unequal-WaterHeater (gal/min)'] = Draw_Profile['Flow Rate (gpm)'] * Calculate_Fraction_Cold_ThroughDWHR(Draw_Profile['Flow Rate (gpm)'], Draw_Profile['Mains Temperature (deg F)'], Temperature_Shower, Temperature_WaterHeater, 'Unequal_WaterHeater') #Use the Calculate_Fraction_Cold_ThroughDWHR function to identify the flow rate of water through the cold side of the DWHR device
    Draw_Profile['Potable Flow Unequal-WaterHeater (gal)'] = Draw_Profile['Mixed Water Volume (gal)'] * Calculate_Fraction_Cold_ThroughDWHR(Draw_Profile['Flow Rate (gpm)'], Draw_Profile['Mains Temperature (deg F)'], Temperature_Shower, Temperature_WaterHeater, 'Unequal_WaterHeater') #Perform the same calculation for the volume instead of the flow rate

    Draw_Profile['Effectiveness Generic Unequal_WaterHeater (-)'] = polyval2d(Draw_Profile['Flow Rate (gpm)'], Draw_Profile['Potable Flow Rate Unequal-WaterHeater (gal/min)'], Coefficients_Generic_Vertical_Unequal) * Effectiveness_Rated #Call the polyval2d function using the specified flow rates and the loaded coefficients to identify the effectiveness correction factor for this scenario. Then multiply by the rated effectiveness to find the actual effectiveness in this draw
    Draw_Profile['Savings Generic Unequal_WaterHeater (Btu)'] = Draw_Profile['Effectiveness Generic Unequal_WaterHeater (-)'] * Draw_Profile['Potable Flow Unequal-WaterHeater (gal)'] * Density_Water * SpecificHeat_Water * (Temperature_Drain_Inlet - Draw_Profile['Mains Temperature (deg F)']) #Multiply the effectiveness by the total available heat to identify the amount of energy saved

    Savings_Generic_Unequal_WaterHeater = Draw_Profile['Savings Generic Unequal_WaterHeater (Btu)'].sum()/100000 #Sum the total energy savings and divide by 100000 to convert to therms

    #Unequal-Fixture calculations

    #Initilizing iterative loop calculations
    Draw_Profile['Fixture Cold Temperature Generic (deg F)'] = Draw_Profile['Mains Temperature (deg F)'] #Make the starting assumption that the temperature of water exiting the cold side of the device is equal to the mains inlet temperature
    Draw_Profile['Flow_Delta Generic (gal/min)'] = 9999 #Flow_Delta represents the difference between the assumed cold-side flow rate and the calcualted cold-side flow rate. This is needed as unequal-shower configurations change the cold-side outlet temperature, which changes the flow rate, which changes the temperature, ... . It's an iterative solution that is considered solved when this difference is very small
    Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] = Draw_Profile['Flow Rate (gpm)'] - Draw_Profile['Hot Water Flow Rate (gpm)'] #Set the cold side flow rate equal to the total flow rate - the hot flow rate

    Draw_Profile['FlowRate Effectiveness Minimum'] = 0.5 #Set a lower bound to the flow rate used to calculate the effectiveness. This helps stabilize the iterative solution, brings about faster results
    Draw_Profile['FlowRate Effectiveness Maximum'] = 7.5 #Set an upper bound to the flow rate used to calculate the effectiveness. This helps stabilize the iterative solution, bring about faster results

    j = 0 #I'm not sure that this is still useful

    #Repeat this process iteratively until the solution is achieved
    while Draw_Profile['Flow_Delta Generic (gal/min)'].max() >= 0.01: #Convergence is identified when the cold side flow rate before and after calculations differ by less than 0.01 gal/min. Continue these calculations until that is achieved
        
        Draw_Profile['FlowRate_Cold_Effectiveness Generic'] = np.where(Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] >= Draw_Profile['FlowRate Effectiveness Minimum'], Draw_Profile['Flow_Cold_ThroughDWHR_Generic'], Draw_Profile['FlowRate Effectiveness Minimum']) #In all rows where the cold side flow rate is larger than the minimum, keep the flow rate. Otherwise, replace the flow rate with the minimum
        Draw_Profile['FlowRate_Cold_Effectiveness Generic'] = np.where(Draw_Profile['FlowRate_Cold_Effectiveness Generic'] <= Draw_Profile['FlowRate Effectiveness Maximum'], Draw_Profile['FlowRate_Cold_Effectiveness Generic'], Draw_Profile['FlowRate Effectiveness Maximum'])     #In all rows where the cold side flow rate is larger than the maximum, keep the flow rate. Otherwise, replace the flow rate with the maximum

        Draw_Profile['FlowRate_Draw_Effectiveness Generic'] = np.where(Draw_Profile['Flow Rate (gpm)'] >= Draw_Profile['FlowRate Effectiveness Minimum'], Draw_Profile['Flow Rate (gpm)'], Draw_Profile['FlowRate Effectiveness Minimum']) #Performs the same process for the mixed water flow rates
        Draw_Profile['FlowRate_Draw_Effectiveness Generic'] = np.where(Draw_Profile['FlowRate_Draw_Effectiveness Generic'] <= Draw_Profile['FlowRate Effectiveness Maximum'], Draw_Profile['FlowRate_Draw_Effectiveness Generic'], Draw_Profile['FlowRate Effectiveness Maximum'])
    
        Draw_Profile['FlowRate_Minimum'] = Draw_Profile['Flow Rate (gpm)'].combine(Draw_Profile['Flow_Cold_ThroughDWHR_Generic'], min, 0) #Ensure that under no circumstances is a flow rate less than 0 gal/min used     
        
        Draw_Profile['Effectiveness_Draw_Generic Equal, Flow=Cold'] = FourthOrder(Coefficients_Generic_Vertical_Equal, Draw_Profile['Flow_Cold_ThroughDWHR_Generic']) * Effectiveness_Rated #Calculates the effectiveness of the device under equal flow conditions, with the flow rate equal to the cold side flow rate. Per Eqn 5.7 on pg 104 of Manouchehri's thesis
        
        Draw_Profile['HeatRecovered_Draw Generic Equal, Flow=Cold'] = Draw_Profile['Effectiveness_Draw_Generic Equal, Flow=Cold'] * Draw_Profile['FlowRate_Minimum'] * Density_Water * SpecificHeat_Water * (Temperature_Drain_Inlet - Draw_Profile['Mains Temperature (deg F)']) * Draw_Profile['Duration (min)'] #Calcualte the energy recovered in the draw using Q_dot = m_dot * C_p * dT
            
        Draw_Profile['HeatRecoveryRate_Equal Generic (Btu/min)'] = Draw_Profile['Effectiveness_Draw_Generic Equal, Flow=Cold'] * Draw_Profile['FlowRate_Minimum'] * Density_Water * SpecificHeat_Water * (Temperature_Drain_Inlet - Draw_Profile['Mains Temperature (deg F)']) #Calcualte the heat recovery rate of each draw in Btu/min
            
        Draw_Profile['Flow Ratio Generic'] = Draw_Profile['Flow Rate (gpm)'] / Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] #Calculate the flow ratio for this draw, per Ramin Manoucherri's paper and calculation method
            
        Draw_Profile['HeatRecoveryRate_Unequal-Fixture Generic (Btu/min)'] = Draw_Profile['HeatRecoveryRate_Equal Generic (Btu/min)'] * (0.3452 * Draw_Profile['Flow Ratio Generic'].apply(log) + 1) #Calculate the heat recovery rate per Ramin Manoucheri's paper and method
        
        Draw_Profile['HeatRecovered_Draw_Generic_Unequal-Fixture (Btu)'] = Draw_Profile['HeatRecoveryRate_Unequal-Fixture Generic (Btu/min)'] * Draw_Profile['Duration (min)'] #Calculate the energy recovered in each draw
        
        Draw_Profile['Fixture Cold Temperature Generic (deg F)'] = Draw_Profile['HeatRecoveryRate_Unequal-Fixture Generic (Btu/min)'] / (Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] * SpecificHeat_Water* Density_Water) + Draw_Profile['Mains Temperature (deg F)'] #Calculate the outlet temperature of the DWHR device in this draw
        
        Draw_Profile['Flow_Cold_ThroughDWHR_PostCalcuations Generic'] = Calculate_Fraction_Cold_ThroughDWHR(Draw_Profile['Flow Rate (gpm)'], Draw_Profile['Fixture Cold Temperature Generic (deg F)'], Temperature_Shower, Temperature_WaterHeater, 'Unequal_Fixture') * Draw_Profile['Flow Rate (gpm)'] #Calculate the cold-side flow rate through the DWHR device using the newly calculated cold-side outlet temperature
            
        Draw_Profile['Flow_Delta Generic (gal/min)'] = abs(Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] - Draw_Profile['Flow_Cold_ThroughDWHR_PostCalcuations Generic']) #Calculate the difference between the cold-side flow rate before and after calculations

        Draw_Profile['Flow_Cold_ThroughDWHR_Generic'] = Draw_Profile['Flow_Cold_ThroughDWHR_PostCalcuations Generic'] #Set the cold-side flow rate through the device equal to the new cold-side flow rate through the device, so it's treated as the starting point in the next iteration

        j += 1

    Savings_Generic_Unequal_Fixture = Draw_Profile['HeatRecovered_Draw_Generic_Unequal-Fixture (Btu)'].sum()/100000 #Sum the energy savings and convert from Btu to therms

#%%----------------SAVE DATA----------------------------

    Draw_Profile.to_csv(i + '_Analyzed.csv', index = False) #Save the performed calcualtsion to a new file with the same name followed by '_Analyzed'


