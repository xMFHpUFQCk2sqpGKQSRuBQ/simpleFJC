import numpy as np
from numba.decorators import jit
from numba import prange, njit
from gaussxw.gaussxw import gaussxwab  

def _int_gauss(f, a, b, N, error=0.0):
    k = np.arange(N)
    s, w = gaussxwab(N, a, b)
    result = f(s[k])
    while w.ndim != result.ndim:
        w = w[:, None]
    return np.nansum((w * result), axis=0)

int_gauss = jit(_int_gauss)

def _int_romb(f, a, b, m, error=0.0):
    N = 3 
    prevI = np.array([int_tra(f, a, b, N)])
    for i in range(1, m + 1):
        nextI = np.empty([i + 1], dtype=np.ndarray)
        nextI[0] = int_tra(f, a, b, N * 2 ** i, prevI[0])
        for j in range(1, i + 1):
            e = (nextI[j - 1] - prevI[j - 1]) / (4 ** j - 1)
            nextI[j] = nextI[j - 1] + e
        prevI = nextI
        if (abs(e) < error).all():
            break
    return nextI[-1]

int_romb = jit(_int_romb)

def _int_tra(f, a, b, N, prevI=None, error=0.0):
    if error != 0:
        startN = 10
        i = 1
        I = int_tra(f, a, b, startN)
        I2 = int_tra(f, a, b, startN * 2 ** i, I)
        
        while (((np.abs(I2 - I) / 3 > error).any()) and (startN < N)):
            i += 1
            I = I2
            I2 = int_tra(f, a, b, startN * 2 ** i, I)            
        return I2
    k = np.arange(1, N)
    h = (b - a) / N
    if prevI is None:
        s = h / 2 * (f(a) + f(b))
    else:
        s = prevI / 2
        k = k[::2]
    return s + h * np.sum(f(a + k * h), axis=0)

int_tra = jit(_int_tra)

def _int_simp(f, a, b, N, error=0):
    k = np.arange(1, N)
    h = (b - a) / N
    return h / 3 * (f(a) + f(b) + 4 * np.nansum(f(a + k[::2] * h), axis=0) + 2 * np.nansum(f(a + k[1::2] * h), axis=0))

int_simp = jit(_int_simp)

# this is the original 3D permute
# for general use
def _permutexyz(x=0, y=0, z=0):
    isx = isinstance(x, np.ndarray)
    isy = isinstance(y, np.ndarray)
    isz = isinstance(z, np.ndarray)
    
    if isx and isy and isz:
        x = x[:, None, None]
        y = y[None, :, None]
        z = z[None, None, :]
    elif isx and isy:
        x = x[:, None]
        y = y[None, :]
    elif isy and isz:
        y = y[:, None]
        z = z[None, :]
    elif isx and isz:
        x = x[:, None]
        z = z[None, :]
    return x, y, z

permutexyz = jit(_permutexyz)

# the original permute up to 6 dimensions
# for general 6D use
def _permutexyz6(x=0, y=0, z=0, xi=0, yi=0, zi=0):
    isx  = isinstance(x, np.ndarray)
    isy  = isinstance(y, np.ndarray)
    isz  = isinstance(z, np.ndarray)
    isxi = isinstance(xi, np.ndarray)
    isyi = isinstance(yi, np.ndarray)
    iszi = isinstance(zi, np.ndarray)
        
    dim   = [ x , y , z , xi, yi, zi]
    isdim = np.array([ isx , isy , isz , isxi, isyi, iszi], dtype=np.bool)
    count = np.sum(isdim)

    if count == 6:        
        dim[0]  =  dim[0][:, None, None, None, None, None]
        dim[1]  =  dim[1][None, :, None, None, None, None]
        dim[2]  =  dim[2][None, None, :, None, None, None]
        dim[3]  =  dim[3][None, None, None, :, None, None]
        dim[4]  =  dim[4][None, None, None, None, :, None]
        dim[5]  =  dim[5][None, None, None, None, None, :]
    elif count == 5:
        aindex = 0
        for index in range(6):
            if isdim[index]:
                if aindex == 0:
                    dim[index] = dim[index][:, None, None, None, None]
                elif aindex == 1:
                    dim[index] = dim[index][None, :, None, None, None]
                elif aindex == 2:
                    dim[index] = dim[index][None, None, :, None, None]
                elif aindex == 3:
                    dim[index] = dim[index][None, None, None, :, None]
                elif aindex == 4:
                    dim[index] = dim[index][None, None, None, None, :]
                aindex += 1
    elif count == 4:
        aindex = 0
        for index in range(6):
            if isdim[index]:
                if aindex == 0:
                    dim[index] = dim[index][:, None, None, None]
                elif aindex == 1:
                    dim[index] = dim[index][None, :, None, None]
                elif aindex == 2:
                    dim[index] = dim[index][None, None, :, None]
                elif aindex == 3:
                    dim[index] = dim[index][None, None, None, :]
                aindex += 1
    elif count == 3:
        aindex = 0
        for index in range(6):
            if isdim[index]:
                if aindex == 0:
                    dim[index] = dim[index][:, None, None]
                elif aindex == 1:
                    dim[index] = dim[index][None, :, None]
                elif aindex == 2:
                    dim[index] = dim[index][None, None, :]
                aindex += 1    
    elif count == 2:
        aindex = 0
        for index in range(6):
            if isdim[index]:
                if aindex == 0:
                    dim[index] = dim[index][:, None]
                elif aindex == 1:
                    dim[index] = dim[index][None, :]
                aindex += 1  

    return dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]

permutexyz6 = jit(_permutexyz6)

# strict 3D permute 
# useful for gauss/monte carlo
def _permutexyz3D(x,y,z):
    x = x[:, None, None]
    y = y[None, :, None]
    z = z[None, None, :]

    return x, y, z

permutexyz3D = jit(_permutexyz3D)

# strict 6D permute 
# useful for gauss/monte carlo
def _permutexyz6D(x,y,z,xi,yi,zi):
    x  =  x[:, None, None, None, None, None]
    y  =  y[None, :, None, None, None, None]
    z  =  z[None, None, :, None, None, None]
    xi = xi[None, None, None, :, None, None]
    yi = yi[None, None, None, None, :, None]
    zi = zi[None, None, None, None, None, :]

    return x, y, z, xi, yi, zi

permutexyz6D = jit(_permutexyz6D)

# original 3D gauss integrator function
def gauss3D(f, bounds, steps=35, error=0.0):
    def getintegrand(Y, Z):
        def integrand(X):
            x, y, z = permutexyz3D(X, Y, Z)
            return f(x, y, z)
        return integrand
    def intx(Z):
        def inner(y):
            return int_gauss(getintegrand(y, Z), bounds[0][0], bounds[0][1], steps, error=error)
        return inner
    def inty(z):
        return int_gauss(intx(z), bounds[1][0], bounds[1][1], steps, error=error)
    return int_gauss(inty, bounds[2][0], bounds[2][1], steps, error=error)

# original 3D simp integrator function
def simp3D(f, bounds, steps=35, error=0.0):
    def getintegrand(Y, Z):
        def integrand(X):
            x, y, z = permutexyz(X, Y, Z)
            return f(x, y, z)
        return integrand
    def intx(Z):
        def inner(y):
            return int_simp(getintegrand(y, Z), bounds[0][0], bounds[0][1], steps, error=error)
        return inner
    def inty(z):
        return int_simp(intx(z), bounds[1][0], bounds[1][1], steps, error=error)
    return int_simp(inty, bounds[2][0], bounds[2][1], steps, error=error)

# original 3D trap integrator function
def trap3D(f, bounds, steps=35, error=0.0):
    def getintegrand(Y, Z):
        def integrand(X):
            x, y, z = permutexyz(X, Y, Z)
            return f(x, y, z)
        return integrand
    def intx(Z):
        def inner(y):
            return int_tra(getintegrand(y, Z), bounds[0][0], bounds[0][1], steps, error=error)
        return inner
    def inty(z):
        return int_tra(intx(z), bounds[1][0], bounds[1][1], steps, error=error)
    return int_tra(inty, bounds[2][0], bounds[2][1], steps, error=error)

def _integrate(object f,double[:,:] bounds, int steps=10 ** 5, str method="trap", double error=0):
    int_func = int_tra
    method = method.lower()
    if method in ["trap", "trapezoid"]:
        int_func = int_tra
    elif method in ["romb", "romberg"]:
        int_func = int_romb
    elif method in ["simp", "simpson"]:
        int_func = int_simp
    elif method in ["gauss", "gaussian"]:
        int_func = int_gauss
    if len(bounds) == 1:
        def integrand(x):
            return f(x)
        return int_func(integrand, bounds[0][0], bounds[0][1], steps, error=error)
    elif len(bounds) == 2:
        def getintegrand(Y):
            def integrand(X):
                x, y, _ = permutexyz(X, Y)
                return f(x, y)
            return integrand 
        def intx(y):
            return int_func(getintegrand(y), bounds[0][0], bounds[0][1], steps, error=error)
        return int_func(intx, bounds[1][0], bounds[1][1], steps, error=error)
    elif len(bounds) == 3:
        def getintegrand(Y, Z):
            def integrand(X):
                x, y, z = permutexyz(X, Y, Z)
                return f(x, y, z)
            return integrand
        def intx(Z):
            def inner(y):
                return int_func(getintegrand(y, Z), bounds[0][0], bounds[0][1], steps, error=error)
            return inner
        def inty(z):
            return int_func(intx(z), bounds[1][0], bounds[1][1], steps, error=error)
        return int_func(inty, bounds[2][0], bounds[2][1], steps, error=error)

integrate = jit(_integrate)

"""
# unfinished 6D gauss integral
def gauss6D(f, bounds, steps=35, error=0.0):
    def getintegrand(Y, Z, XI, YI, ZI):
        def integrand(X):
            x, y, z, xi, yi, zi = permutexyz6D(X, Y, Z, XI, YI, ZI)
            return f(x, y, z, xi, yi, zi)
        return integrand
    
    def intx(Z):
        def inner(y):
            return int_gauss(getintegrand(y, Z, XI, YI, ZI), bounds[0][0], bounds[0][1], steps, error=error)
        return inner
    def inty(z):
        return int_gauss(intx(z), bounds[1][0], bounds[1][1], steps, error=error)
    def intz():
        def inner():
            return int_gauss(intx(z), bounds[2][0], bounds[2][1], steps, error=error)
        return inner
    def intxi(ZI):
        def inner(yi):
            return int_gauss(intx(z), bounds[3][0], bounds[3][1], steps, error=error)
        return inner
    def intyi(zi):
        return int_gauss(intxi(zi)), bounds[4][0], bounds[4][1], steps, error=error)
        pass
    return int_gauss(intyi, bounds[5][0], bounds[5][1], steps, error=error)

    """