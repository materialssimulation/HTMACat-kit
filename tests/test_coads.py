import pytest
import os
import filecmp
import hashlib

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
        # if not filecmp.cmp(file1_path, file2_path, shallow=False):
        #     with open(file1_path,'r') as f:
        #         print('file1')
        #         print(f.read())
        #     with open(file2_path,'r') as f:
        #         print('file2')
        #         print(f.read())
        #     print(f"文件 {file} 不同")
        #     return False

        file1_md5 = calculate_md5(file1_path)
        file2_md5 = calculate_md5(file2_path)
        if file1_md5 != file2_md5:
            with open(file1_path,'r') as f:
                print('file1')
                print(f.read())
            with open(file2_path,'r') as f:
                print('file2')
                print(f.read())
            print(f"文件 {file} 不同")
            return False


    # 文件夹中的文件全部相同
    return True

def calculate_md5(file_path):
    with open(file_path, 'rb') as f:
        content = f.read()
        md5_hash = hashlib.md5(content)
        return md5_hash.hexdigest()