import pandas as pd
from qmplot import qqplot
import matplotlib.pyplot as plt
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--infile', help='tsv文件,包含header')
@click.option('--val-col', help='y轴值列名')
@click.option('--title', default="Q-Q plot", help='输入QQ图的标题')
@click.option('--outfile', help='png格式文件名')
def main(infile, val_col, title, outfile):
    df = pd.read_table(infile, sep="\t", usecols=[val_col], dtype={val_col: float})
    df = df.dropna(how="any", axis=0)  # If any NA values are present, drop that row
    # Create a Q-Q plot
    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(data=df[val_col], marker="o", title=title, xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", dpi=600, figname=outfile, s=2, ax=ax)
if __name__ == '__main__':
    main()
