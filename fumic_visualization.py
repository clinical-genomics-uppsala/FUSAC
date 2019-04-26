
# !/usr/bin/env python3

# FUMIC - FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data
# By Hugo Swenson, with assistance from Patrik Smeds and Claes Edenvall
# Made for Klinisk Genetik, Uppsala Akademiska Sjukhus 2019

# Imports modules
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def vcf_count(vcf_file):
    ffpe_ind = 0
    ind = 0
    var_lst = []
    ffpe_lst = []
    perc_lst = []
    pos_lst = []
    change_lst = []
    for record in vcf_file.fetch():
        try:
            r_f = record.filter
            for f_val in r_f:
                if f_val == "FFPE":
                    # ref, mut, ffpe, n, del, ref_paired, var_paired, ref_single, var_single]
                    ffpe_ind += 1
                    rec_samp = record.samples
                    pos_lst.append(str(record.pos))
                    change_lst.append(str(record.ref) + ">" + str(list(record.alts)[0]))
                    for sample in rec_samp:
                        samp_dat = list(record.samples[sample]['UMI'])
                        rec_str = samp_dat[0]
                        rec_splt = rec_str.split(';')
                        rec_var = int(rec_splt[1])
                        rec_ffpe = int(rec_splt[2])
                        rec_ref = int(rec_splt[0])
                        rec_n = int(rec_splt[3])
                        rec_del = int(rec_splt[4])
                        var_lst.append(rec_var)
                        ffpe_lst.append(rec_ffpe)
                        if rec_var == 0:
                            perc_rat = 0
                        else:
                            perc_rat = (rec_ffpe/(rec_ref + rec_var + rec_ffpe + rec_n + rec_del)) * 100
                        perc_lst.append(perc_rat)
            else:
                ind += 1
        except KeyError as e:
            print("ERROR: The requested filter tag " + str(e) + " does not exist")
    fig, sub_plt = plt.subplots()
    x_ind = np.arange(len(var_lst))
    bar_width = 0.3
    sub_plt.bar(x_ind - bar_width / 2, var_lst, width=0.3, color='b', align='center', label="Var")
    sub_plt.bar(x_ind + bar_width/2, ffpe_lst, width=0.3, color='r', align='center', label="FFPE")
    sub_plt.set_xlabel('Called Var')
    sub_plt.set_ylabel('Counts')
    sub_plt.set_xticks(x_ind)
    sub_plt.legend()
    fig.savefig("performance_plot.png")
    plt.close()

    sns.set(style="ticks", color_codes=True)
    # Percentage line-plot of hte read-depth
    perc_df = pd.DataFrame({'Variant': x_ind, 'FFPE': perc_lst})
    perc_plot = sns.lineplot(x='Variant', y='FFPE', data=perc_df)
    perc_plot.set(xlabel='Variant', ylabel='FFPE/Depth', title="Percentage plot of FFPE vs. Depth")
    perc_plot.figure.savefig("percentage_plot.pdf", format='pdf', dpi=1000)
    plt.close()

    # Scatter plot of read depth for each variant
    var_type_lst = ['Var']*len(var_lst)
    ffpe_type_lst = ['FFPE']*len(var_lst)
    type_lst = var_type_lst + ffpe_type_lst
    depth_lst = var_lst + ffpe_lst
    pos2_lst = list(range(1, len(pos_lst)+1)) + list(range(1, len(pos_lst)+1))
    sca1_df1 = pd.DataFrame({'Variant': pos2_lst, 'Depth': depth_lst, 'Type': type_lst})
    sns.scatterplot(x='Variant', y='Depth', hue='Type', data=sca1_df1)
    plt.savefig("sns_scatter_1.pdf", format='pdf', dpi=1000)

    # VS. Scatter plot of read depth for each variant
    sca2_df = pd.DataFrame({'Var_Depth': var_lst, 'FFPE_Depth': ffpe_lst, 'Variant': pos_lst})
    sca_sns_plot_2 = sns.catplot(x="Var_Depth", y="FFPE_Depth", hue='Variant', data=sca2_df, height=4, aspect=2)
    for ax in sca_sns_plot_2.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)
    sns_leg = sca_sns_plot_2._legend
    sns_leg.set_title('Variant')
    plt.setp(sns_leg.get_texts(), fontsize=8)
    for lh in sns_leg.legendHandles:
        lh._sizes = [5]
    sca_sns_plot_2.savefig("sns_scatter_2.pdf", format='pdf', dpi=1000)

    # VS. Scatter plot of read depth for each variant
    sca3_df = pd.DataFrame({'Var_Depth': var_lst, 'FFPE_Depth': ffpe_lst, 'Base_Change': change_lst})
    sca_sns_plot_3 = sns.catplot(x="Var_Depth", y="FFPE_Depth", hue='Base_Change', data=sca3_df, height=4, aspect=2)
    for ax in sca_sns_plot_3.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)
    sns_leg_2 = sca_sns_plot_3._legend
    sns_leg_2.set_title('Base Change')
    plt.setp(sns_leg_2.get_texts(), fontsize=8)
    for lh in sns_leg_2.legendHandles:
        lh._sizes = [5]
    sca_sns_plot_3.savefig("sns_scatter_3.pdf", format='pdf', dpi=1000)


        # out_file = open("fumic_count.txt", "w")
        # out_file.write("Total number of reads: " + str(ind) + "\n")
        # out_file.write("Total number of FFPE-artefacts found; " + str(ffpe_ind) + "\n")
        # out_file.write("Fraction of FFPE-artefacts in sample: " + str(ffpe_ind/ind) + "\n")
        # out_file.close()


def main():
    # Reads in the bam file
    vcf_file = pysam.VariantFile("fumic_output.vcf", "r")
    vcf_count(vcf_file)

if __name__ == "__main__":
    main()