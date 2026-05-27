import gzip

def comp_check(*files):
    for item in files:
        with open(item,'rb') as f:
            return f.read(3) == b'\x1f\x8b\x08'


