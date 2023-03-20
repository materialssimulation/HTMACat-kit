import pytest
import os
import filecmp

@pytest.mark.parametrize("dir1, dir2", [
    ("./example/coads", "./tests/results/coads"),
    ("./example/coads-SML", "./tests/results/coads-SML")
])

def test_coads(dir1, dir2):
    cwd = os.getcwd()
    os.chdir(dir1)
    os.system('coads')
    # 获取文件夹中的文件列表
    files1 = [os.path.join(dp, f) for dp, dn, filenames in os.walk(dir1) for f in filenames if os.path.isfile(os.path.join(dp, f)) and f.endswith(".vasp")]
    files2 = [os.path.join(dp, f) for dp, dn, filenames in os.walk(dir2) for f in filenames if os.path.isfile(os.path.join(dp, f)) and f.endswith(".vasp")]

# 对比文件夹中的文本文件
    for file1, file2 in zip(files1, files2):
        assert filecmp.cmp(file1, file2, shallow=False)
    os.chdir(cwd)