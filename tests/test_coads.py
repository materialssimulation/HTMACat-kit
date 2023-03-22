import pytest
import os
import filecmp

@pytest.mark.parametrize("dir1, dir2", [
    ("./example/coads", "./tests/results/coads"),
    ("./example/coads-SML", "./tests/results/coads-SML")
])

def test_ads(dir1, dir2):
    cwd = os.getcwd()
    os.chdir(dir1)
    os.system('coads')
    os.chdir(cwd)
    # 获取文件夹中的文件列表
    assert compare_folders(dir1, dir2)

def compare_folders(folder1, folder2):
    # 获取两个文件夹下的所有文件列表
    files1 = set(os.listdir(folder1))
    files2 = set(os.listdir(folder2))

    # 判断文件夹中的文件是否完全相同
    if not files1 == files2:
        print(folder1+"缺少文件")
        return False

    # 逐个比较文件内容
    for file in files1:
        file1_path = os.path.join(folder1, file)
        file2_path = os.path.join(folder2, file)
        if not filecmp.cmp(file1_path, file2_path, shallow=False):
            print(f"文件 {file} 不同")
            return False

    # 文件夹中的文件全部相同
    return True