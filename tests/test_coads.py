import pytest
import os
import filecmp
import hashlib
from ase.io.vasp import read_vasp
from ase import Atoms
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
        if file == 'Model' or file == 'StrucInfo':
            continue
        file1_path = os.path.join(folder1, file)
        file2_path = os.path.join(folder2, file)

        # 读取两个VASP文件
        atoms1 = read_vasp(file1_path)
        atoms2 = read_vasp(file2_path)

        # 将VASP结构对象转换为ASE结构对象
        atoms1 = Atoms(cell=atoms1.get_cell(),
                    positions=atoms1.get_positions(),
                    symbols=atoms1.get_chemical_symbols(),
                    pbc=True)
        atoms2 = Atoms(cell=atoms2.get_cell(),
                    positions=atoms2.get_positions(),
                    symbols=atoms2.get_chemical_symbols(),
                    pbc=True)

        # 判断两个ASE结构对象是否一样
        if atoms1 == atoms2:
            # print('The two structures are identical.')
            continue
        else:
            with open(file1_path,'r') as f:
                print('file1')
                print(f.read())
            with open(file2_path,'r') as f:
                print('file2')
                print(f.read())
            print(f"文件 {file} 不同")
            # print('The two structures are different.')
            return False
    return True

        # if not filecmp.cmp(file1_path, file2_path, shallow=False):
        #     with open(file1_path,'r') as f:
        #         print('file1')
        #         print(f.read())
        #     with open(file2_path,'r') as f:
        #         print('file2')
        #         print(f.read())
        #     print(f"文件 {file} 不同")
        #     return False

        # file1_md5 = calculate_md5(file1_path)
        # file2_md5 = calculate_md5(file2_path)
        # if file1_md5 != file2_md5:
        #     with open(file1_path,'r') as f:
        #         print('file1')
        #         print(f.read())
        #     with open(file2_path,'r') as f:
        #         print('file2')
        #         print(f.read())
        #     print(f"文件 {file} 不同")
        #     return False


    # 文件夹中的文件全部相同
    # return True

# def calculate_md5(file_path):
#     with open(file_path, 'rb') as f:
#         content = f.read()
#         content = content.strip()
#         md5_hash = hashlib.md5(content)
#         return md5_hash.hexdigest()