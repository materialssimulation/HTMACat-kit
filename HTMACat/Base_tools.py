# Function1: a is inclued in b return True
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


# Function2: a and b are equal return True
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

# Function3: Find the corresponding key and value from the dictionary to form a new dictionary
def get_new_dict(list, dict):
    new_dict = {}
    for key, value in dict.items():
        if key in list:
            new_dict[key] = value
    return new_dict

if __name__ == "__main__":
    a = [1, 2, 3, 4]
    b = [1, 2, 3, 4, 5]
    print(isequal(a, b))
