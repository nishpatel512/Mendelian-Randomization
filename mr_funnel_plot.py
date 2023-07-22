import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from blank_plot import blank_plot
from mr_singlesnp import mr_singlesnp
from extract_instruments import extract_instruments
from extract_outcome_data import extract_outcome_data
from harmonise_data import harmonise_data


def mr_funnel_plot(singlesnp_results):
    res = []

    # Define a color map for SNPs
    cmap = plt.get_cmap("tab10")

    for group, d in singlesnp_results.groupby(["id.exposure", "id.outcome"]):
        if sum(~d["SNP"].str.contains("All")) < 2:
            res.append(blank_plot("Insufficient number of SNPs"))
            continue

        am = d.loc[d["SNP"].str.contains("All"), "SNP"]
        d["SNP"] = d["SNP"].str.replace("All - ", "")
        am = am.str.replace("All - ", "")

        plt.figure()
        # Assign unique colors to SNPs based on the color map
        colors = [cmap(i % 10) for i in range(len(am))]
        plt.scatter(x=d.loc[~d["SNP"].isin(am), "b"], y=1 / d.loc[~d["SNP"].isin(am), "se"], color="blue")
        for snp, color in zip(am, colors):
            plt.axvline(x=d.loc[d["SNP"] == snp, "b"].values[0], color=color, label=snp)

        plt.xlabel(r"$\beta_{IV}$")
        plt.ylabel(r"$1/$SE_{IV}$")
        plt.legend(title="MR Method")
        plt.tight_layout()

        res.append(plt.gca())

    plt.show()

# exp = extract_instruments("ieu-a-2")
# out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

# dat = harmonise_data(exp, out)

# res = mr_singlesnp(dat)

# mr_funnel_plot(res)
