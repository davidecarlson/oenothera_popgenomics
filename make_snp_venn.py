#!/usr/bin/env python

# make venn diagram from shared and exclusive snp counts in ../stats/vcf-compare.stats

from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt

ax = venn2(subsets = (100917, 38745, 35890), set_labels = ('biennis', 'elata'), set_colors = ('r', 'b'),
alpha = 0.5)
venn2_circles(subsets = (100917, 38745, 35890), linewidth=2)
#plt.show()
plt.savefig("../figures/snp_venn.eps", bbox_inches='tight')
plt.clf()
