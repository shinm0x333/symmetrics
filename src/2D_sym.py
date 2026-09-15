import cmath, math
from functools import cmp_to_key

eps = 1e-5
def _identical(a, b):
    if type(a) in [complex, float, int]:
            return abs(a-b) < eps
    raise

def _z_sym(z: complex, arg: float) -> complex:
    return cmath.exp(2j*arg) * z.conjugate()

def _check(c_i: list[complex], arg: float) -> bool:
    c2_i = [_z_sym(z, arg) for z in c_i]
    c_i.sort(key=lambda x: (x.real, x.imag))
    c2_i.sort(key=lambda x: (x.real, x.imag))
    return all(_identical(c1, c2) for c1, c2 in zip(c_i, c2_i))

def method_i(p_i: list) -> tuple[float, float, float]:
    c_i = [complex(*p) for p in p_i]
    avg = sum(c_i)/len(c_i)
    for exp in 2, 3:
        var = sum((z-avg)**exp for z in c_i)
        if _identical(var, 0): continue
        for n in range(exp):
            arg = (cmath.phase(var) + n*cmath.pi) / exp
            if _check(c_i, arg):
                return (-math.sin(arg), math.cos(arg), math.sin(arg)*avg.real - math.cos(arg)*avg.imag)
    raise

def _half(z: complex) -> int:
    x = z.real; y = z.imag
    return 0 if (y > 0 or (y == 0 and x >= 0)) else 1

def _cross_cmp(a: complex, b: complex) -> int:
    ha, hb = _half(a), _half(b)
    if ha != hb:
        return -1 if ha < hb else 1

    cross = a.real*b.imag - a.imag*b.real
    if cross > 0: return -1
    if cross < 0: return 1
    return 0

def _compress(l: list[tuple], start_cnt: int = 0) -> list[tuple]:
    l.sort()
    cnt = start_cnt
    ret = [(start_cnt, l[0][1])]
    for i in range(1, len(l)):
        if not _identical(l[i-1][0], l[i][0]):
            cnt += 1
        ret.append((cnt, l[i][1]))
    ret.sort(key=lambda x: x[1])
    return ret

MOD = 10**9+7
def _segment_palindromes(l: list[int], sz: int) -> list[int]:
    N = len(l)
    B = max(l)+1
    orig = []
    h = 0
    DEL = pow(B, sz, MOD)
    for i in range(N):
        h = (h*B + l[i])%MOD
        if i >= sz:
            h = (h - l[i-sz]*DEL)%MOD
        if i >= sz-1:
            orig.append(h)
    h = 0
    rev = []
    for i in range(N-1, -1, -1):
        h = (h*B + l[i])%MOD
        if N-i-1 >= sz:
            h = (h - l[i+sz]*DEL)%MOD
        if N-i-1 >= sz-1:
            rev.append(h)
    rev.reverse()
    return [i for i in range(N-sz+1) if orig[i] == rev[i]]
        
def _line_eq(a: complex, b: complex) -> tuple[float, float, float]:
    return (a-b).imag, (b-a).real, a.real*b.imag - b.real*a.imag

def _bisector_eq(a: complex, b: complex) -> tuple[float, float, float]:
    return ((a-b).real, (a-b).imag, (abs(b)**2-abs(a)**2)/2)

def method_ii(p_i: list) -> list[tuple[float, float, float]]:
    c1_i = [complex(*p) for p in p_i]
    c = sum(c1_i)/len(p_i)
    c1_i = [z-c for z in c1_i]
    c_i = [z for z in c1_i if not _identical(0, z)]; N = len(c_i)
    c_i.sort(key=cmp_to_key(_cross_cmp))
    l1 = [(abs(c_i[i]), i) for i in range(N)]
    l2 = []
    for i in range(N):
        j = (i+1)%N
        l2.append(((cmath.phase(c_i[i])-cmath.phase(c_i[j])) % (2*math.pi), i))
    l1 = _compress(l1)
    l2 = _compress(l2, len(l1))
    l_star = [l2[i>>1][0] if i&1 else l1[i>>1][0] for i in range(2*N)] * 2
    l_star.pop(0); l_star.pop()
    seg = _segment_palindromes(l_star, 2*N-1)
    ret = []
    for s in seg:
        m = s + N  
        k = (m % (N*2)) // 2  
        if m&1:
            ret.append(_bisector_eq(c+c_i[k], c+c_i[(k+1)%N]))
        else:
            ret.append(_line_eq(c, c+c_i[k]))
    return ret

def _get_c_in_radius(c_i: list[complex], r: float) -> complex:
    circ = [z for z in c_i if abs(z) < r or _identical(abs(z), r)]
    return sum(circ)/len(circ)

def method_iii(p_i: list) -> tuple[float, float, float]:
    N = len(p_i)
    c_i = [complex(*p) for p in p_i]
    c = sum(c_i)/N
    c_i = [z-c for z in c_i]
    r_i = [abs(z) for z in c_i]; r_i.sort()
    for dx in 0, -1, 1, -2, 2:
        if not 0 <= N//2+dx < N: continue
        c2 = _get_c_in_radius(c_i, r_i[N//2+dx])+c
        if not _identical(c, c2):
            return _line_eq(c, c2)
    raise

# TODO: 결과값 약분하기?
if __name__ == "__main__":
    p_i = []
    tc = int(input())
    for _ in range(tc):
        x, y = map(int, input().split())
        p_i.append((x, y))
    print(method_i(p_i))
    print(method_ii(p_i))
    print(method_iii(p_i))
