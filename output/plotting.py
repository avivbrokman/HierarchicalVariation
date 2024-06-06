#%% libraries
from numpy import arange
from numpy.random import uniform
from pandas import DataFrame, concat
import plotnine
from plotnine import ggplot, aes, geom_point, geom_jitter, ylim, theme, geom_rect, xlab, ylab

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

centers = [centers1] * 12 + [centers2] * 1 + [centers3] * 1 + [centers2] * 2 + [centers3] * 1 + [centers2] * 3 + [centers4] * 1

partitions = [partition1] * 12 + [partition2] * 1 + [partition3] * 1 + [partition2] * 2 + [partition3] * 1 + [partition2] * 3 + [partition4] * 1

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
possible_partitions = [partition1, partition2, partition3]
partition_strings = ["1/2, 3/4", "1/2, 3, 4", "1, 2, 3, 4"]

unique_p1 = sorted(data['p1'].unique())
rectangles = DataFrame({
    'xmin': [x - 1/(2*40) for x in unique_p1],  # Small offset for the width
    'xmax': [x + 1/(2*40) for x in unique_p1],  # Small offset for the width
    'ymin': [0] * len(unique_p1),  # Covers the whole y range
    'ymax': [1] * len(unique_p1),  # Covers the whole y range
    #'fill': ['#D3D3D3' if i % 2 == 0 else '#B0C4DE' for i in range(len(unique_p1))]  # Alternating colors
    'fill': [partition_strings[possible_partitions.index(el)] for el in partitions]  # Alternating colors
})

#%% plotting
from plotnine import labs, element_text, geom_vline, geom_segment, coord_cartesian

g = ggplot()
g = g + geom_point(data, aes(x = 'shift_', y = 'center', shape = 'patch'), show_legend = False)
g = g + ylim((0,1))
g = g + theme(figure_size=(10, 4))
g = g + geom_rect(data = rectangles, mapping = aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax', fill='fill'), alpha=0.5)
g = g + xlab('$p_1$') + ylab('Center')
#g = g + guides(patch = False)
g = g + labs(fill = "Distribution")

g = g + theme(
    axis_title = element_text(size = 14),
    axis_text = element_text(size = 12),
    legend_text = element_text(size = 12),
    legend_title = element_text(size = 14),
    plot_title = element_text(size = 16)
)

verts = [el - 1/80 for el in p1s]
verts.pop(0)
# =============================================================================
# g = g + geom_vline(xintercept = verts, color = 'grey')
# =============================================================================
for el in verts:
    g = g + geom_segment(aes(x = el, xend = el, y = 0, yend = 1), color = "grey")
    

g


