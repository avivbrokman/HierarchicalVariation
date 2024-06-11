#%% libraries
import numpy as np
from numpy import arange, add, array
from numpy.random import uniform
from pandas import DataFrame, concat, Categorical
import plotnine
from plotnine import ggplot, aes, geom_point, geom_jitter, ylim, theme, geom_rect, xlab, ylab, labs, element_text, geom_vline, geom_segment, coord_cartesian, scale_fill_brewer, scale_fill_hue, qplot

import json
import os
from os.path import join
import glob


#%% get all partitions

#%% data

centers1 = [0.1, 0.9, 0.1, 0.9]
partition1 = [[1,2], [3,4]]
centers2 = [0.1, 0.9, 0.9, 0.9]
partition2 = [[1,2], [3], [4]]
centers3 = centers2
partition3 = [[1], [2], [3], [4]]
centers4 = [0.9] * 4
partition4 = partition3


p1s = arange(0.5, 1 + 1/40, 1/40)

centers = [centers1] * 12 + [centers2] * 1 + [centers2] * 1 + [centers2] * 2 + [centers2] * 1 + [centers2] * 3 + [centers4] * 1

partitions = [partition1] * 12 + [partition2] * 1 + [partition2] * 1 + [partition2] * 2 + [partition2] * 1 + [partition2] * 3 + [partition4] * 1

def unlist(list_of_lists):
    
    return [item for sublist in list_of_lists for item in sublist]

def make_plotting_mini_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    
    
    

    if len(partition) == 1:
        shift_ = [p1] * fecundity
    elif len(partition) == 2:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - 0.5 * pch_width] * len(el)
            elif i == 1:
                shift_ += [p1 + 0.5 * pch_width] * len(el)
    elif len(partition) == 3:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - pch_width] * len(el)
            elif i == 1:
                shift_ += [p1] * len(el)
            elif i == 2:
                shift_ += [p1 + pch_width] * len(el)
    elif len(partition) == 4:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - 1.5 * pch_width] * len(el)
            elif i == 1:
                shift_ += [p1 - 0.5 * pch_width] * len(el)
            elif i == 2:
                shift_ += [p1 + 0.5 * pch_width] * len(el)
            elif i == 3:
                shift_ += [p1 + 1.5 * pch_width] * len(el)
    
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shift_": shift_})
    
def make_plotting_dataframe(p1s, centers, partitions, pch_width):
    mini_frames = [make_plotting_mini_dataframe(el_p1, el_cent, el_part, pch_width) for el_p1, el_cent, el_part in zip(p1s, centers, partitions)]
    
    data_frame = concat(mini_frames, ignore_index = True)
    
    return data_frame

data = make_plotting_dataframe(p1s, centers, partitions, pch_width = 0.006)
#%% rectangles
possible_partitions = [[1, 2, 3, 4]], [[1, 2, 3], [4]], [[1, 2], [3, 4]], [[1, 2], [3], [4]], [[1], [2], [3], [4]]
partition_strings = ["1/2/3/4", "1/2/3, 4", "1/2, 3/4", "1/2, 3, 4", "1, 2, 3, 4"]

unique_p1 = sorted(data['p1'].unique())
rectangles = DataFrame({
    'xmin': [x - 1/(2*40) for x in unique_p1],  
    'xmax': [x + 1/(2*40) for x in unique_p1],  
    'ymin': [0] * len(unique_p1),  
    'ymax': [1] * len(unique_p1),  
    'fill': [partition_strings[possible_partitions.index(el)] for el in partitions]  # Alternating colors
})





 


#%% read in data
os.chdir('/Users/avivbrokman/Documents/Kentucky/Grad School/ms_project/branching1')

def read_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

p1s = arange(0.5, 1 + 1/40, 1/40)


#%%
def unlist(list_of_lists):
    
    return [item for sublist in list_of_lists for item in sublist]

def calculate_patch_shifts(num_patches, pch_width):
    
    return [(i - 0.5 * (num_patches - 1)) * pch_width for i in range(num_patches)]

def calculate_partition_shift(partition, pch_width):
    
    patch_shifts = calculate_patch_shifts(len(partition), pch_width)
    
    offspring_shifts = [[patch_shifts[i]] * len(el) for i, el in enumerate(partition)]
    
    return unlist(offspring_shifts)

def make_plotting_mini_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    offsets = calculate_partition_shift(partition, pch_width)
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shifted_p1": add(p1s, offsets)})

#%% rectangles
def partition2string(partition):
    string_list = ['/'.join([str(el) for el in patch]) for patch in partition]
    return ', '.join(string_list)
    

def make_rectangles(partitions, p1s):

    rectangles = DataFrame({
        'xmin': [x - 1/(2*40) for x in p1s],  
        'xmax': [x + 1/(2*40) for x in p1s],  
        'ymin': [0] * len(p1s),  
        'ymax': [1] * len(p1s),  
        'fill': [partition2string(el) for el in partitions]
    })
    
    partition_strings = ["1, 2, 3, 4", "1/2/3, 4", "1/2, 3/4", "1/2, 3, 4", "1/2/3/4"]
    rectangles['fill'] = Categorical(rectangles['fill'], categories = partition_strings, ordered = True)
     
    return rectangles




#%% binomial alpha1 = 0.6

def make_datasets(pattern):
    files = glob.glob(pattern, recursive = True)
    mini_datasets = []
    p1_partitions = []
    extinction_probabilities = []
    for file in files:
        results = read_json(file)
        p1 = results['p1']
        centers = results['centers']
        partition = results['partition']
        mini_dataset = make_plotting_mini_dataframe(p1, centers, partition, 0.006)
        
        mini_datasets.append(mini_dataset)
        p1_partitions.append((p1, partition))
        extinction_probabilities.append(results['extinction_probability'])
        
    data = concat(mini_datasets, ignore_index = True)  
    
    p1s, partitions = zip(*p1_partitions)
    rectangles = make_rectangles(partitions, p1s)
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities})
    
    
    return data, rectangles, extinction_data

def make_extinction_probability_dataset(pattern):
    files = glob.glob(pattern, recursive = True)
    extinction_probabilities = []
    partitions = []
    p1s = []
    for file in files:
        results = read_json(file)
# =============================================================================
#         print(results)
# =============================================================================
        for el in results:
            print(el)
            p1s.append(el['p1'])
            partitions.append(partition2string(el['partition']))
            extinction_probabilities.append(el['extinction_probability'])
        
    
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities, "partition": partitions})
    
    
    return extinction_data


#%% make average_dataset
def make_plotting_mini_avg_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    offsets = calculate_partition_shift(partition, pch_width)
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shifted_p1": add(p1s, offsets)})

def make_avg_datasets(pattern):
    files = glob.glob(pattern, recursive = True)
    mini_datasets = []
    p1_partitions = []
    extinction_probabilities = []
    for file in files:
        results = read_json(file)
        p1 = results['p1']
        centers = results['centers']
        partition = results['partition']
        mini_dataset = make_plotting_mini_dataframe(p1, centers, partition, 0.006)
        
        mini_datasets.append(mini_dataset)
        p1_partitions.append((p1, partition))
        extinction_probabilities.append(results['extinction_probability'])
        
    data = concat(mini_datasets, ignore_index = True)  
    
    p1s, partitions = zip(*p1_partitions)
    rectangles = make_rectangles(partitions, p1s)
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities})
    
    
    return data, rectangles, extinction_data
#%% plotting

def plot(data, rectangles):
    g = ggplot()
    g = g + ylim((0,1))
    g = g + theme(figure_size=(10, 4))
    
    g = g + geom_rect(data = rectangles, mapping = aes(xmin = 'xmin', xmax = 'xmax', ymin = 'ymin', ymax = 'ymax', fill = 'fill'), alpha = 1) #alpha = 0.5
    
    g = g + geom_point(data, aes(x = 'shifted_p1', y = 'center', shape = 'patch'), show_legend = False)
    
    g = g + xlab('$p_1$') + ylab('Center')
    g = g + labs(fill = "Distribution")
# =============================================================================
#     g = g + scale_fill_brewer(palette='Set1')
# =============================================================================
    
    g = g + scale_fill_hue()

    g = g + theme(
        axis_title = element_text(size = 14),
        axis_text = element_text(size = 12),
        legend_text = element_text(size = 12),
        legend_title = element_text(size = 14),
        plot_title = element_text(size = 16)
    )
    
    verts = rectangles.xmin.sort_values()
    verts = verts[1:]
    for el in verts:
        g = g + geom_segment(aes(x = el, xend = el, y = 0, yend = 1), color = "#B8B8B8")
    
# =============================================================================
#     my_colors = ["#e6a8a4", "#d5e6a4", "#a4e6c2", "#a4bae6", "#dda4e6"]
# =============================================================================
    my_colors = ["#FFA8A8", "#d5e6a4", "#a4e6c2", "#96E7FF", "#dda4e6"] #e6a8a4 #A1CDE6 #99FFFF
    g = g + scale_fill_manual(values = my_colors)
    
    print(g)
    return

def extinction_plot(data):
    g = ggplot(data, aes(x = "p1", y = "extinction_probability"))
    g = g + geom_point(aes(color = "partition"))
    print(g)
    return 

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.9__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.9__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_4__beta1_2__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_4__beta1_2__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)
#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_7__beta1_3__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_7__beta1_3__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_20__beta1_10__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_20__beta1_10__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_23__beta1_7__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_23__beta1_7__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.95__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.95__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_1.8__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_1.8__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_2__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_2__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_3__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_3__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_4__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_4__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)


#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_5__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_5__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_7.5__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_7.5__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)


#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_10__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_10__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%

def plot_extinction_probability(extinction_data):
    
    q = qplot("p1", "extinction_probability", extinction_data)
    print(q)
    return

#%% raw
def get_


#%% running raw
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/robust_output.json'

files = glob.glob(pattern, recursive = True)
raw_results = [read_json(el) for el in files]

data = raw_results[0]
data = array(data)
np.mean(data, axis = 1)


