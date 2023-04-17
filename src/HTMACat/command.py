import os
from HTMACat.model.Construct_adsorption_yaml import *
from pathlib import *
import shutil
import typer

htmat = typer.Typer()
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

def main():
    htmat()


@htmat.callback(invoke_without_command=True, context_settings=CONTEXT_SETTINGS)
def main_command():
    print(f"HTMACat-Kit Version: 1.0.3")
    print('A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes.')


@htmat.command()
def ads(in_dir: str = typer.Option('./', '-i', '--inputdir',
                                   help="relative directory of input file"),
        out_dir: str = typer.Option('./', '-o', '--outputdir',
                                    help="relative directory of output file")):
    """
    Construct adsorption configuration
    """
    wordir = Path(in_dir).resolve()
    outdir = Path(out_dir).resolve()
    StrucInfo = 'config.yaml'
    if not outdir == wordir:
        outdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(wordir / StrucInfo, outdir)
        os.chdir(outdir)
    Construct_adsorption_yaml(StrucInfo)

# def ads():
#     parser = argparse.ArgumentParser(description='High-throughput single adsorption modeling')
#     parser.add_argument('-i', '--inputdir', type=str, default="./", help="The folder path of input files")
#     parser.add_argument('-o', '--outputdir', type=str, default="./", help="The folder path of output files")
#     # 添加了定向输出到指定路径 by Yuxiao-Lan 2023/03/03
#     args = parser.parse_args()
#     wordir = Path(args.inputdir).resolve()
#     outdir = Path(args.outputdir).resolve()
#     StrucInfo = 'StrucInfo'
#     ModelInfo = 'Model'
#     if not outdir == wordir:
#         outdir.mkdir(parents=True, exist_ok=True)
#         shutil.copy(wordir / StrucInfo, outdir)
#         shutil.copy(wordir / ModelInfo, outdir)
#         os.chdir(outdir)
#     with open(StrucInfo, 'r+') as f:
#         for i, item in enumerate(f):
#             Construct_adsorption(item, ModelInfo)
#
#
# def ads_yaml():
#     parser = argparse.ArgumentParser(description='High-throughput single adsorption modeling')
#     parser.add_argument('-i', '--inputdir', type=str, default="./", help="The folder path of input files")
#     parser.add_argument('-o', '--outputdir', type=str, default="./", help="The folder path of output files")
#     # 添加了定向输出到指定路径 by Yuxiao-Lan 2023/04/13
#     args = parser.parse_args()
#     wordir = Path(args.inputdir).resolve()
#     outdir = Path(args.outputdir).resolve()
#     StrucInfo = 'config.yaml'
#     if not outdir == wordir:
#         outdir.mkdir(parents=True, exist_ok=True)
#         shutil.copy(wordir / StrucInfo, outdir)
#         os.chdir(outdir)
#     Construct_adsorption_yaml(StrucInfo)
#
#
# def coads():
#     parser = argparse.ArgumentParser(description='High-throughput co-adsorption modeling')
#     parser.add_argument('-i', '--inputdir', type=str, default="./", help="The folder path of input files")
#     parser.add_argument('-o', '--outputdir', type=str, default="./", help="The folder path of output files")
#     # 添加了定向输出到指定路径 by Yuxiao-Lan 2023/03/03
#     args = parser.parse_args()
#     args = parser.parse_args()
#     wordir = Path(args.inputdir).resolve()
#     outdir = Path(args.outputdir).resolve()
#     StrucInfo = 'StrucInfo'
#     ModelInfo = 'Model'
#     if not outdir == wordir:
#         outdir.mkdir(parents=True, exist_ok=True)
#         shutil.copy(wordir / StrucInfo, outdir)
#         shutil.copy(wordir / ModelInfo, outdir)
#         os.chdir(outdir)
#     # spec_ads,spec_ads_stable=get_site_stable(Efile,Ecut=-0.1)
#     # spec_ads_stable={'NH3':[1,2],'NH2':[2],'NH':[2,4],'N':[2,4],'O':[2,4],'OH':[2,4],'NO':[2,4],'H2O':[1],'H':[2,4]}
#     spec_ads_stable = {'NH3': [1], 'NH2': [2], 'NH': [2, 4], 'N': [2, 4], 'O': [2, 4], 'OH': [2, 4], 'NO': [2, 4],
#                        'H2O': [1], 'H': [2, 4]}
#     # spec_ads_stable={'NH2':[2],'NH':[3],'NO':[3],'NH3':[1],'N':[3],'O':[3],'OH':[3]}
#     Construct_coadsorption(StrucInfo, ModelInfo, spec_ads_stable)
