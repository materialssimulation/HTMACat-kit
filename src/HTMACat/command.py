import os
from HTMACat.model.Construct_adsorption_yaml import *
from HTMACat.IO import print_templator, out_templator_file
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


@htmat.command()
def templator():
    """
    Print out input templator
    """
    print_templator()
    out_templator_file()
