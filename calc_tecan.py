from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys

interval = 30  ## measurement interval in minutes
win_max = .75  # GMR max
win_min = .15  # GMR min

filename = sys.argv[1]


def ParseFromFile(filename):  ## puts all data into "well_dict"
       
    well_dict = {}

    f = open(filename)
    for line in f:
        if line.startswith('<>'):
            continue  # Skip header lines.
            
        row_data = line.strip().split("\t")
        row_label = row_data[0]
            
        # First entry in line is the row label.
        for i, cell in enumerate(row_data[1:]):
            cell_label = '%s%d' % (row_label, i+1)
            if len(cell_label) == 2:
                cell_label = cell_label[0]+'0'+cell_label[1]  ## adds 0 to the begenning of cell labels <10

            cell_data = float(cell.strip()) or '0.0'

            well_dict.setdefault(cell_label, []).append(cell_data)
        
    return well_dict


def SmoothData(data_dict):  ## smoothes the curves by taking a rolling mean

    smooth_dict = {}

    for key in data_dict:
        
        smooth_dict[key] = []
        
        for pos, val in enumerate(data_dict[key]):

            if pos == 0 or pos == len(data_dict[key])-1:  #  don't try to average on the first or last value in the timecourse
                continue

            ## smooth data

            smooth_dict[key].append(np.mean([data_dict[key][pos-1], data_dict[key][pos], data_dict[key][pos+1]]))
            
            ## ensure monotonocity

            if len(smooth_dict[key]) > 1:

                if smooth_dict[key][-1] < smooth_dict[key][-2]:
                
                    smooth_dict[key][-1] = smooth_dict[key][-2]

    return data_dict

def ZeroData(data):   ## zeroes data, ensures no values less than 0

    for well_key in data:

        sorted = np.sort(data[well_key])
            
        data[well_key]-= np.mean(sorted[2:5])  # subtracts the mean of the 4 lowest values from all wells, not including the first 2

        for pos, val  in enumerate(data[well_key]):  ## ensure no values are less than 0
            if val < 0:
                data[well_key][pos] = 0.0

    return data


def CalcGrowthRates(data):  ## calculates a list of pointe estimates of growth rates to be used fror GRM

    #Geometric Mean Rate (doublings per hour) - geometric mean of the growth rates within a fixed window


    GMR_dict = {}
    inter_dict = {} ## holds the rate values for each well


    x = np.array((range(len(data[list(data.keys())[0]]))))*np.array([30])  # creates a 1D array of each timepoint
    
    for well_key in data:

        f = interp1d(x, data[well_key], kind='cubic')  ## creates an interpolation function 

        inter_values = f(range(x[-1]))  ## creates a list of interpolated values for every minute of the timecourse

        in_inters = []

        for val in inter_values:


            if val > win_min and val < win_max:   ## these numbers determine the bounding of the window that is used to calculated GMR.  
                in_inters.append(val)    



        GMR_dict[well_key] = (np.log(win_max/win_min)/(len(in_inters)))*60  ## *60 to convert from minutes to hours


    return GMR_dict


def Get_rate_lag_eff(smooth_dict):
    
    eff_dict = {}
    max_rate_dict = {}
    lag_dict = {}
    for well in smooth_dict:
        
        sorted_smooth_data = sorted(smooth_dict[well])
        eff = np.mean(sorted_smooth_data[-3:]) - np.mean(sorted_smooth_data[:3])
        eff_dict[well] = eff

        all_rates = []
        log_data = np.log(np.array(smooth_dict[well]))

        rate_list = []
        for pos, val in enumerate(log_data[:-1]):

 

            if smooth_dict[well][pos] < .04:  ## the minimum detection limit for max rate.  .04 seems a good compromise between eliminating noise and getting at the true maximum
                all_rates.append([0,0])
                continue
            
            
            regress_data = stats.linregress([(pos-1)*interval, (pos)*interval, (pos+1)*interval],[log_data[pos-1], log_data[pos], log_data[pos+1]])

            rate = regress_data[0]*60  ## * 60 to make it per hour
            y_int = regress_data[1]

            
            all_rates.append([rate, y_int])
            
            rate_list.append(rate)
            
        all_rates.sort(key=lambda z: z[0])

        

        max_rates = all_rates[-4:]      ## take the mean of the 4 highest rates
            
        mean_max_rate = np.mean([x[0] for x in max_rates])
        mean_yint = np.mean([x[1] for x in max_rates])

        mean_lag_time = (-mean_yint/mean_max_rate)/2  ## this would be in minutes, but rete has already been multiplied by 60, so this is the same as dividing the lag time by 60.  IT is /2 for 30 minutes per timepoint
        
        max_rate_dict[well] = mean_max_rate
        lag_dict[well] = mean_lag_time

        x = np.array(range(len(smooth_dict[well])))*.5
        

#        plt.plot(rate_list)
#        plt.plot(x, smooth_dict[well])
#        plt.show()
        
    return max_rate_dict, lag_dict, eff_dict



def PrintOutput(GMR_dict, max_rate_dict, lag_dict, eff_dict):

    print ("Well", 'Max Rate', "Lag", "Efficiency", "GMR")

    for well_key in sorted(GMR_dict):
        print(well_key, max_rate_dict[well_key], lag_dict[well_key], eff_dict[well_key], GMR_dict[well_key])



#### START PROGRAM ######

raw_data = ParseFromFile(filename) # this returns a dictionary with well labels as keys and data as a list in each key
smooth_data = SmoothData(raw_data)
zeroed_data = ZeroData(smooth_data)
GMR_dict = CalcGrowthRates(zeroed_data)
max_rate_dict, lag_dict, eff_dict =  Get_rate_lag_eff(smooth_data)

PrintOutput(GMR_dict, max_rate_dict, lag_dict, eff_dict)




