import pytest
import os
import filecmp

@pytest.mark.parametrize("dir1, dir2", [
    ("./example/Formic_acid_ads", "./tests/results/Formic_acid_ads"),
    ("./example/Formic_acid_ads-SML", "./tests/results/Formic_acid_ads-SML")
])

def test_ads(dir1, dir2):
    cwd = os.getcwd()
    os.chdir(dir1)
    os.system('ads')
    os.chdir(cwd)
    # 比较文件夹的文件
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
            with open(file1_path.'r') as f:
                print('file1')
                print(f)
            with open(file_path.'r') as f:
                print('file2')
                print(f)
            print(f"文件 {file} 不同")
            return False

    # 文件夹中的文件全部相同
    return True