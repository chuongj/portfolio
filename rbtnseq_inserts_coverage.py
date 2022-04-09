import sys
import numpy as np
import matplotlib.pyplot as plt
import time

# To access coverage of RB-Tnseq libraries, calculate how many transposon inserts we observed in genes.
# Plot a histogram of inserts per gene

pooled_readfile = sys.argv[1]
gff = '/usr2/people/carlyweiss/SparScerinfo/YS2+CBS432+plasmid_clean' ## GFF dict Carly used for annotations.

def parse_pooled_outfile(pooled_readfile):
    gene_dict = {} # will be {gene: [scer_inserts, spar_inserts]}

    f = open(gff)

    for line in f:
        line = line.split("\t")
        gene = line[0]
        #print("print gene = line[0]", gene)
        species=gene[:2]
        #print("species",species)
        y_name = gene[2:-1]
        #print("y_name",y_name)
        #exit()
        
        # put the y-names from the gff into a dictionary #                                                                                                                                    

        if y_name not in gene_dict:
            gene_dict[y_name] = [0,0] 

    f.close()

    f = open(pooled_readfile)
    next(f)
    for line in f:
        line = line.split("\t")
        gene=line[4]
        #print("gene",gene)
        species=gene[:2]
        #print("species",species)
        y_name=gene[2:-1] # cut off the strand 'W' or 'C' in case reversed in species
        #print("y_name",y_name)
        

        #add inserts counts to dictionary of y-names#

        if species=="sc":
            gene_dict[y_name][0]+=1
        if species=="sp":
            gene_dict[y_name][1]+=1
            
    #print(gene_dict['YBR050'])
    return gene_dict


def get_quick_stats(gene_dict):
    total_genes_with_1_sc=0
    total_genes_with_2_sc=0
    total_genes_with_5_sc=0
    total_genes_with_10_sc=0
    total_genes_with_1_sp=0
    total_genes_with_2_sp=0
    total_genes_with_5_sp=0
    total_genes_with_10_sp=0
    total_genes_with_1_both=0
    total_genes_with_2_both=0
    total_genes_with_5_both=0
    total_genes_with_10_both=0

    for gene in gene_dict:
        if gene_dict[gene][0]>=10: #gene_dict[gene][0] is num of scer inserts 
            total_genes_with_1_sc+=1
            total_genes_with_2_sc+=1
            total_genes_with_5_sc+=1
            total_genes_with_10_sc+=1
        elif gene_dict[gene][0]>=5:
            total_genes_with_1_sc+=1
            total_genes_with_2_sc+=1
            total_genes_with_5_sc+=1
        elif gene_dict[gene][0]>=2:
            total_genes_with_1_sc+=1
            total_genes_with_2_sc+=1
        elif gene_dict[gene][0]==1:
            total_genes_with_1_sc+=1

        if gene_dict[gene][1]>=10: #gene_dict[gene][1] is num of spar inserts
            total_genes_with_1_sp+=1
            total_genes_with_2_sp+=1
            total_genes_with_5_sp+=1
            total_genes_with_10_sp+=1
        elif gene_dict[gene][1]>=5:
            total_genes_with_1_sp+=1
            total_genes_with_2_sp+=1
            total_genes_with_5_sp+=1
        elif gene_dict[gene][1]>=2:
            total_genes_with_1_sp+=1
            total_genes_with_2_sp+=1
        elif gene_dict[gene][1]==1:
            total_genes_with_1_sp+=1

        
        if gene_dict[gene][0]>=10 and gene_dict[gene][1]>=10:
            total_genes_with_1_both+=1
            total_genes_with_2_both+=1
            total_genes_with_5_both+=1
            total_genes_with_10_both+=1
        elif gene_dict[gene][0]>=5 and gene_dict[gene][1]>=5:
            total_genes_with_1_both+=1
            total_genes_with_2_both+=1
            total_genes_with_5_both+=1
        elif gene_dict[gene][0]>=2 and gene_dict[gene][1]>=2:
            total_genes_with_1_both+=1
            total_genes_with_2_both+=1
        elif gene_dict[gene][0]==1 and gene_dict[gene][1]==1:
            total_genes_with_1_both+=1


    fraction_1_sc=total_genes_with_1_sc/len(gene_dict)
    fraction_2_sc=total_genes_with_2_sc/len(gene_dict)
    fraction_5_sc=total_genes_with_5_sc/len(gene_dict)
    fraction_10_sc=total_genes_with_10_sc/len(gene_dict)

    fraction_1_sp=total_genes_with_1_sp/len(gene_dict)
    fraction_2_sp=total_genes_with_2_sp/len(gene_dict)
    fraction_5_sp=total_genes_with_5_sp/len(gene_dict)
    fraction_10_sp=total_genes_with_10_sp/len(gene_dict)

    fraction_1_both=total_genes_with_1_both/len(gene_dict)
    fraction_2_both=total_genes_with_2_both/len(gene_dict)
    fraction_5_both=total_genes_with_5_both/len(gene_dict)
    fraction_10_both=total_genes_with_10_both/len(gene_dict)

    print("coverage at 1, 2, 5, and 10 inserts each:")
    print("sc: "+str(100*fraction_1_sc)+"%, "+str(100*fraction_2_sc)+"%, "+str(100*fraction_5_sc)+"%, "+str(100*fraction_10_sc)+"%")
    print("sp: "+str(100*fraction_1_sp)+"%, "+str(100*fraction_2_sp)+"%, "+str(100*fraction_5_sp)+"%, "+str(100*fraction_10_sp)+"%")
    print("both: "+str(100*fraction_1_both)+"%, "+str(100*fraction_2_both)+"%, "+str(100*fraction_5_both)+"%, "+str(100*fraction_10_both)+"%")

    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16])
    
    x = [1,2,5,10]
    plt.xlabel("# inserts per allele")
    plt.ylabel("%"+" genes covered")
    plt.ylim(0,100)
    plt.xticks(ticks=x, labels=x)
    
    plt.plot(x,[100*fraction_1_sc, 100*fraction_2_sc, 100*fraction_5_sc, 100*fraction_10_sc], marker='o',color = '#3288bd', label = "Scer")
    plt.plot(x,[100*fraction_1_sp, 100*fraction_2_sp, 100*fraction_5_sp, 100*fraction_10_sp], marker='o',color = '#fdae61', label = "Spar")
    plt.plot(x,[100*fraction_1_both, 100*fraction_2_both, 100*fraction_5_both, 100*fraction_10_both], marker='o',color = "#2ca02c", label = "both") #green

    plt.legend(frameon=False)

    pooled_readfile_id=pooled_readfile.split('.')[0]
    t = time.localtime()
    timestamp = time.strftime('.%H.%M', t)

    plt.savefig('cov_insert_plot_'+pooled_readfile_id+".pdf", format='pdf')






def histogram_floor(gene_dict, pooled_readfile):
    floor=[]#floor is the lower of the two insert counts
    for gene in gene_dict:
        floor.append(min(gene_dict[gene]))

    num_bins=50

    n, bins, patches = plt.hist(floor, num_bins, facecolor='blue', density=1, alpha=0.5, range=(0, 50))

    plt.xlabel('#inserts in both species (min of #sc, #sp)')
    plt.ylabel('proportion of inserts')
    plt.title('Histogram of coverage for '+pooled_readfile)

    plt.subplots_adjust(left=0.15)

    pooled_readfile_id=pooled_readfile.split('.')[0]

    t = time.localtime()
    timestamp = time.strftime('.%H.%M', t)
    plt.savefig('coverage_propor_hist_for_'+pooled_readfile_id+timestamp+".pdf", format='pdf')

def histogram_logplot(gene_dict, pooled_readfile):
    floor=[]#floor is the lower of the two insert counts
    for gene in gene_dict:
        floor.append(min(gene_dict[gene]))

    num_bins=50

    n, bins, patches = plt.hist(floor, num_bins, facecolor='orange', density=1, log=1, alpha=0.5, range=(0, 50))

    plt.xlabel('#inserts in both species')
    plt.ylabel('log proportion of inserts')
    plt.title('Histogram of coverage for '+pooled_readfile)

    plt.subplots_adjust(left=0.15)

    pooled_readfile_id=pooled_readfile.split('.')[0]

    t = time.localtime()
    timestamp = time.strftime('.%H.%M', t)
    plt.savefig('coverage_proporlog_hist_for_'+pooled_readfile_id+timestamp+".pdf", format='pdf')


##START##
gene_dict=parse_pooled_outfile(pooled_readfile)
get_quick_stats(gene_dict)


#histogram_floor(gene_dict, pooled_readfile)
#histogram_logplot(gene_dict, pooled_readfile)