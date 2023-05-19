import os
from HTMACat.model.Construct_adsorption_yaml import *
from HTMACat.IO import print_templator, out_templator_file
from pathlib import *
import shutil
import typer
from rich import print

htmat = typer.Typer(add_completion=False)
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def main():
    htmat()


@htmat.callback(
    invoke_without_command=True,
    no_args_is_help=True,
    epilog="""+--------------------------------------------------------------------------------------------------------+\n
|                                           HTMACat-Kit                                                  |\n
|                                          Version: 1.0.4                                                |\n
|    A high-throughput modeling, calculation, and analysis framework for catalytic reaction processes.   |\n
|              More information, please visit https://stanfordbshan.github.io/HTMACat-kit/.              |\n
+--------------------------------------------------------------------------------------------------------+""",
    context_settings=CONTEXT_SETTINGS,
)
def main_command():
    pass


@htmat.command(context_settings=CONTEXT_SETTINGS)
def ads(
    in_dir: str = typer.Option("./", "-i", "--inputdir", help="relative directory of input file"),
    out_dir: str = typer.Option(
        "./", "-o", "--outputdir", help="relative directory of output file"
    ),
):
    """Construct adsorption configuration."""
    print("Construct adsorption configuration ... ...")
    wordir = Path(in_dir).resolve()
    outdir = Path(out_dir).resolve()
    StrucInfo = "config.yaml"
    if not outdir == wordir:
        outdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(wordir / StrucInfo, outdir)
        os.chdir(outdir)
    Construct_adsorption_yaml(StrucInfo)


@htmat.command(context_settings=CONTEXT_SETTINGS)
def templator():
    """Print out input templator."""
    print_templator()
    out_templator_file()
