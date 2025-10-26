import os
from HTMACat.model.Construct_adsorption_yaml import *
from HTMACat.IO import print_templator, out_templator_file, yaml2dict
from HTMACat.CRN import runCRN_net
from HTMACat.CRN import run_crnconfiggen
from HTMACat.Split import coads_split
from HTMACat.Show_net import draw_net
from HTMACat.__version__ import __title__, __version__
from pathlib import *
import shutil
import typer
import json
from rich import print

htmat = typer.Typer(add_completion=False)
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def main():
    htmat()


@htmat.callback(
    invoke_without_command=True,
    no_args_is_help=True,
    epilog=f"""+--------------------------------------------------------------------------------------------------------+\n
|                                           HTMACat-Kit                                                  |\n
|                                          Version: {__version__}                                                |\n
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
):
    """Construct adsorption configuration."""
    import os
    from HTMACat.api import construct_adsorption
    from rich import print

    print("[bold green]Construct adsorption configuration via API...[/bold green]")

    # 找到配置文件路径
    config_path = os.path.join(in_dir, "config.yaml")
    if not os.path.exists(config_path):
        print(f"[red]Error:[/red] config.yaml not found in {in_dir}")
        raise typer.Exit(code=1)

    # 直接传入路径，不需要读取内容
    construct_adsorption(config_yaml=config_path)

    print("[bold green]✅ Adsorption configuration generated successfully![/bold green]")


@htmat.command(context_settings=CONTEXT_SETTINGS)
def templator():
    """Print out input templator."""
    print_templator()
    out_templator_file()

@htmat.command(context_settings=CONTEXT_SETTINGS)
def complete(in_dir: str = typer.Option("./", "-i", "--inputdir", help="relative directory of input file")):
    """Complete config"""
    StrucInfo = "config.yaml"
    os.chdir(in_dir)
    substrate_dict,species_dict,ads_dict = yaml2dict(StrucInfo)
    config = {'StrucInfo':substrate_dict, 'Species':species_dict, 'Model':ads_dict}
    config_str = json.dumps(config)
    with open('./complete_config.json','w',encoding='utf-8') as f:
        f.write(config_str)

@htmat.command(context_settings=CONTEXT_SETTINGS)
def crn():
    """Generate the Chemical Reaction Network."""
    try:
        log_content = runCRN_net()
        with open('CRNGenerator_log.txt', 'w', encoding='utf-8') as f:
            f.write(log_content)
    except Exception as e:
        print(f"Error generating CRN: {e}")


@htmat.command(context_settings=CONTEXT_SETTINGS)
def crngen():
    """Generate structured directories and input files based on CRNGenerator_log.txt"""
    run_crnconfiggen()

@htmat.command(context_settings=CONTEXT_SETTINGS)#lbx
def split(filename,key_atom):
    """split configuration."""
    print("split ... ...")
    coads_split(filename,key_atom)
@htmat.command(context_settings=CONTEXT_SETTINGS)
def drawnet():
    """Draw the Chemical Reaction Network."""
    draw_net()