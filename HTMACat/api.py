import os
import yaml
import tempfile
from pathlib import Path
from HTMACat.model.Construct_adsorption_yaml import Construct_adsorption_yaml
from rich import print


def construct_adsorption(
    config_yaml: str = None,
    StrucInfo=None,
    Species=None,
    Model=None,
    workdir="./"
):
    """
    构建吸附构型。
    支持两种调用方式：
      1. construct_adsorption(config_yaml="config.yaml")
      2. construct_adsorption(StrucInfo=dict, Model=dict, [Species=dict])
    """

    print("[HTMACat] Construct adsorption configuration ...")
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)

    # ============= 模式1：直接读取 config.yaml 文件 =============
    if config_yaml and os.path.exists(config_yaml):
        Construct_adsorption_yaml(config_yaml)
        print("✅ Adsorption configuration generated successfully!")
        return

    # ============= 模式2：使用 Python 字典输入 =============
    if StrucInfo and Model:
        config_data = {"StrucInfo": StrucInfo, "Model": Model}
        if Species:
            config_data["Species"] = Species

        # 临时 YAML 文件（不会污染用户目录）
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".yaml") as tmp:
            yaml.dump(config_data, tmp)
            tmp_path = tmp.name

        # 调用原有逻辑
        Construct_adsorption_yaml(tmp_path)

        # 清理临时文件
        os.remove(tmp_path)
        print("✅ Adsorption configuration generated successfully!")
        return

    # ============= 参数错误 =============
    raise ValueError("❌ You must provide either config_yaml path or StrucInfo+Model dicts.")


