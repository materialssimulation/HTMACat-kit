### Function1: a is inclued in b return True
def isinclude(a, b):
    item = []
    for i in a:
        if i not in b:
            item += [False]
            break
        elif i in b:
            item += [True]
    if len(a) == item.count(True):
        return True
    else:
        return False


### Function2: a and b are equal return True
def isequal(a, b):
    item = []
    if len(a) == len(b):
        for i in range(len(a)):
            if a[i] == b[i]:
                item += [True]
            else:
                item += [False]
        if len(a) == item.count(True):
            return True
        else:
            return False
    else:
        return False

### Function3: 用0.000000代替vasp文件中的-0.000000
def replace(file):
    with open(file,'r') as f:
        content = f.read()
    content = content.replace('-0.000000',' 0.000000')
    with open(file,'w') as f:
        f.write(content)

if __name__ == '__main__':
    a = [1, 2, 3, 4]
    b = [1, 2, 3, 4, 5]
    print(isequal(a, b))
