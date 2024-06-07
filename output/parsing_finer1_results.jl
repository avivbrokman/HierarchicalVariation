###### Bimodal
# All centers are either 0.1 or 0.9 in both hierarchical and nonhierarchical scenarios. The only differences are in how many are 0.9 and which partition is optimal.
# for alpha1 = 0.6: 
    # p1 = 0.5 - 0.675: two 0.9/0.1 patches 
    # p1 = 0.7 - 0.75: one 0.9/0.1 patch, two 0.9/patches
    # p1 = 0.775 - 1, four 0.9 patches
# for alpha1 = 0.9:
    # p1 = 0.5 - 0.775: two 0.9/0.1 patches 
    # p1 = 0.7 - 0.75: one 0.9/0.1 patch, two 0.9/patches
    # p1 = 0.8 - 1, four 0.9 patches
# Overall, patches should either have a 0.1 and a 0.9 specialist or just a 0.9 specialist. As p1 ↑, the distributions shift from two double-specialist patches, to one double-specialist and two 0.9-specialists, to four 0.9-specialist patches. This transition happens later as alpha1 ↑.
# nonhierarchical comparison
    # for alpha1 - 0.6: same results
    # for alpha2 - 0.9: 
        # at p1 = 0.5 - 0.8, two patches each with 0.1 and 0.9 specialists
        # at p1 = 0.9, one patch with 0.1 and 0.9, two patches with just 0.9
        # at p1 = 1, four patches, each with 0.9
    # Overall, almost the same as hierarchical. That's because the two single distributions are honestly fairly similar to each other
####### extreme mode varying mass everywhere else
# All centers are either 0.1 or 0.9 in both hierarchical and nonhierarchical scenarios. The only differences are in how many are 0.9 and which partition is optimal.
# for beta1 = 0.05: 
    # at p1 = 0.5 - 0.8, two patches each with 0.1 and 0.9 specialists
    # at p1 = 0.9, one patch with 0.1 and 0.9, two patches with just 0.9
    # at p1 = 1, four patches, each with 0.9
# for beta1 = 0.95: 
    # at p1 = 0.5, one patch with 0.1 and 0.9, two patches with just 0.9
    # at p1 = 0.6 - 1, four patches, each with 0.9
# Overall, similar to bimodal scenario. Patches should either have a 0.1 and a 0.9 specialist or just a 0.9 specialist. As p1 ↑, the distributions shift from two double-specialist patches, to one double-specialist and two 0.9-specialists, to four 0.9-specialist patches. This transition happens earlier as beta1 ↑.
# nonhierarchical
    # for beta1 = 0.05:
         # at p1 = 0.5 - 0.9, two patches each with 0.1 and 0.9 specialists
        # at p1 = 1, four patches, each with 0.9
    # for beta1 = 0.95
        # all setting went extinct.
        # THINK ABOUT WHAT YOU WANT TO DO HERE!!!

##### extreme mode, varying slope 
# for alpha1 = 1.2:
    # all scenarios went extinct
# for alpha1 = 1.8:
    # for p1 = 0.5 - 0.7, extinction
    # for p1 = 0.8, 95% chance of extinction. four separate patches. one with 1 center at 0.9, one at 0.87, one two at 0.81. !!!NEED TO INVESTIGATE THIS ONE FURTHER
    # for p1 = 0.9, reasonable % extinction. Three patches. One patch with 0.7 (!!confirm), 0.9 centers, and two patches each with 0.9 centers.
    # for p1 = 1, two patches, each with 0.7 (!!confirm), 0.9 centers
# for alpha1 = 2:
    # for p1 = 0.5 - 0.7, extinction
    # for p1 = 0.8, 91% chance of extinction, four separate patches. one with 1 center at 0.9, one at 0.86, one at 0.85, and one at 0.72.
    # for p1 = 0.9, 1, two patches, each with 0.7 (!!confirm), 0.9 centers
# for alpha1 = 5:
    # for p1 = 0.5, 0.6, two patches, each with 0.9/0.1 centers. 
    # for p1 = 0.7, four patches, three with 0.9 centers, and one with 0.1 center
    # for p1 = 0.8, 0.9, three patches, one with 0.9/0.1 and two with 0.9 centers
    # for p1 = 1, two patches, one with 0.9/0.7/0.5 (!!confirm) and one with 0.9 center
# nonhierarchical
    # for alpha1 = 1.2:
        # alls scenarios go extinct
    # for alpha1 = 1.8
        # p1 = 0.5 - 0.7: extinction
        # p1 = 0.8-0.9: four patches of 0.9 centers
        # p1 = 1: two patches of 0.9/0.7 centers
    # alpha1 = 2
        # p1 = 0.5 - 0.6: extinction
        # p1 = 0.7 - 0.9: 4 patches of 0.9 centers
        # p1 = 1: 2 patches of 0.9/0.7 centers
    # alpha1 = 5
        # p1 = 0.5 - 0.6: 2 patches of 0.9/0.1 centers
        # p1 = 0.7 - 0.9: 1 patch of 0.9/0.1 centers, and two patches of 0.9 centers
        # p1 = 1: 1 patch 0.9/0.7/0.5. 1 patch 0.9
# Overall: It looks like this is the pattern:
    # For small alpha1, right when p1 is big enough for nonextinction, the best strategy is 4 patches with overlapping centers from 0.7 to 0.9. As p1 ↑, the optimal strategy turns into patches of 0.9/0.7 and 0.9, with increasing 0.9/0.7 patches as p1 → 1.
    # for larger alpha1, the results were a bit chaotic. p1 = 0.7 results make a pattern elusive (!!MAYBE REDO WITH A DIFFERENT SEED). !!I should definitely run this with a finer step-size for p1. 
    

