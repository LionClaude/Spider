def read_distances(file_name):
    with open(file_name, 'r') as data:
        y = []
        for line in data:
            y.append(float(line))
    return y

def read_tseries(file_name):
    with open (file_name, 'r') as data:
        X = []
        s = []
        for line in data:
            if line != '\n':
                s.append(float(line))
            else:
                X.append(s)
                s = []
    return X

def read_velseries(file_name):
    with open (file_name, 'r') as data:
        X = []
        Z = []
        sx = []
        sz = []
        for line in data:
            if line != '\n':
                vx, vz = line.split(" ")
                sx.append(float(vx))
                sz.append(float(vz))
            else:
                X.append(sx)
                Z.append(sz)
                sx = []
                sz = []
    return X, Z

def read_qvelseries(file_name):
    with open (file_name, 'r') as data:
        X = []
        Xold = []
        Vx = []
        sx = []
        sxold = []
        svx = []
        for line in data:
            if line != '\n':
                x, xold, vx = line.split(" ")
                sx.append(float(x))
                sxold.append(float(xold))
                svx.append(float(vx))
            else:
                X.append(sx)
                Xold.append(sxold)
                Vx.append(svx)
                sx = []
                sxold = []
                svx = []
    return X, Xold, Vx
