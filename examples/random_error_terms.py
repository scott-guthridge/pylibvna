from libvna.cal import CalType, Parameter, Solver
import math
import cmath
import numpy as np
import sys

sqrt_2 = math.sqrt(2.0)


def random_complex(rng, σ, size=None):
    σ /= sqrt_2
    return (rng.normal(0, σ, size) + rng.normal(0, σ, size) * 1j)


class RandomSParameter:
    '''
    Make a random S parameter that varies smoothly with frequency.
    '''
    def __init__(self, rng, fmin, fmax):
        self.a = rng.normal(0, 0.5)
        self.b = rng.normal(0, 0.5)
        self.c = rng.normal(0, 0.5)
        # Choose random undamped natural frequency between
        # fmin / 10 and fmax * 10.
        fmin /= 10.0
        fmax *= 10.0
        f = fmin * 10 ** (math.log10(fmax / fmin) * rng.random())
        self.w0 = 2.0 * math.pi * f
        self.z = 0.05 + abs(rng.normal(0, 0.95))

    def evaluate(self, f):
        s = 2.0j * math.pi * f
        return ((self.a + self.b * 2.0 * self.z * self.w0 * s + self.c * s**2) /
                (self.w0**2 + 2.0 * self.z * self.w0 * s + s**2))


class RandomTL2x2:
    '''
    Make a random length of RG-175 transmission line.
    '''
    def __init__(self, rng):
        self.vf = 2.0 / 3.0             # velocity factor
        self.lm = 2.50e-7               # metal loss: Np/m/Hz**(1/2)
        self.ld = 5.86e-12              # dielectric loss: Np/m/Hz
        self.z0 = 50.0
        self.length = 0.1 + 0.25 * rng.random()   # in meters

    def get_S(self, f: float):
        '''
        Return an S-parameter representing transmission at frequency f
        assuming the characteristic impedance is 50 ohms real.
        '''
        c = 2.9979246e+08               # speed of light (m/s)
        gl = (self.lm * math.sqrt(f) + self.ld * f
              + 2.0j * math.pi * f / (self.vf * c)) * self.length
        return cmath.exp(-gl)


class RandomErrorTerms:
    class System:
        def __init__(self, rng, ctype: CalType, m_rows: int, m_cols: int,
                     s_rows: int, s_cols: int, fmin: float, fmax: float):
            # Find the nominal (for 2x2) number of error terms
            # based on the calibration type.
            if ctype == CalType.T8 or ctype == CalType.U8:
                terms = 8
            elif (ctype == CalType.TE10 or ctype == CalType.UE10
                    or ctype == CalType.UE14 or ctype == CalType.E12):
                terms = 10
            elif ctype == CalType.T16 or ctype == CalType.U16:
                terms = 16
            else:
                assert not "type not handled"

            self.transmission = []
            for port in range(max(s_rows, s_cols)):
                self.transmission.append(RandomTL2x2(rng))
            self.el = np.full((m_rows, m_cols), None, dtype=object)
            self.er = np.full((m_rows, s_rows), None, dtype=object)
            self.et = np.full((s_cols, m_cols), None, dtype=object)
            self.em = np.full((s_cols, s_rows), None, dtype=object)

            # init El terms
            for i in range(m_rows):
                for j in range(m_cols):
                    if i == j or terms > 8:
                        self.el[i, j] = RandomSParameter(rng, fmin, fmax)
            # init Er and Ef leakages if using 16-term
            if terms == 16:
                for i in range(m_rows):
                    for j in range(s_rows):
                        if i != j:
                            self.er[i, j] = RandomSParameter(rng, fmin, fmax)
                for i in range(s_cols):
                    for j in range(m_cols):
                        if i != j:
                            self.et[i, j] = RandomSParameter(rng, fmin, fmax)
            # init Em terms
            for i in range(s_cols):
                for j in range(s_rows):
                    if i == j or terms == 16:
                        self.em[i, j] = RandomSParameter(rng, fmin, fmax)

    def __init__(self, rng, ctype: CalType, rows: int, cols: int,
                 fmin: float, fmax: float):
        has_colsys = ctype == CalType.UE14 or ctype == CalType.E12
        m_rows = rows
        if has_colsys:
            m_cols = 1
            n_systems = cols
        else:
            m_cols = cols
            n_systems = 1
        s_rows = max(rows, cols)
        s_cols = s_rows

        self.rng = rng
        self.has_colsys = has_colsys
        self.m_rows = m_rows
        self.m_cols = m_cols
        self.s_rows = s_rows
        self.s_cols = s_cols
        self.fmin = fmin
        self.fmax = fmax
        self.systems = []
        for sindex in range(n_systems):
            self.systems.append(RandomErrorTerms.System(rng, ctype,
                                m_rows, m_cols, s_rows, s_cols, fmin, fmax))

    def evaluate(self, calset, f_vector, s_matrix, ab=False):
        m_rows = self.m_rows
        m_cols = self.m_cols
        s_rows = self.s_rows
        s_cols = self.s_cols
        s_matrix = [[Parameter.from_value(calset, s_matrix[i][j])
                    for i in range(s_rows)] for j in range(s_cols)]
        n_systems = len(self.systems)
        frequencies = len(f_vector)
        result_cols = max(m_cols, n_systems)
        result = np.empty((frequencies, m_rows, result_cols),
                               dtype=complex)
        for sindex, system in enumerate(self.systems):
            for findex, f in enumerate(f_vector):
                s = np.empty(shape=(s_rows, s_cols), dtype=complex)
                for i in range(s_rows):
                    for j in range(s_cols):
                        s[i, j] = s_matrix[i][j].get_value(f)
                el = np.zeros((m_rows, m_cols), dtype=complex)
                er = np.zeros((m_rows, s_rows), dtype=complex)
                et = np.zeros((s_cols, m_cols), dtype=complex)
                em = np.zeros((s_cols, s_rows), dtype=complex)

                # Add random transmission lines between VNA ports
                # and corresponding DUT ports.
                for port in range(max(m_rows, m_cols)):
                    temp = system.transmission[port].get_S(f)
                    if port < m_rows:
                        er[port, port] = temp
                    if self.has_colsys:
                        if port == sindex:
                            et[port, 0] = temp
                    elif port < m_cols:
                        et[port, port] = temp

                # Add el directivity/leakage error terms
                for i in range(m_rows):
                    for j in range(m_cols):
                        term = system.el[i, j]
                        if term is not None:
                            el[i, j] = term.evaluate(f)

                # Add reflection leakage terms
                for i in range(m_rows):
                    for j in range(s_rows):
                        term = system.er[i, j]
                        if term is not None:
                            er[i, j] = term.evaluate(f)

                # Add transmission leakage terms
                for i in range(s_cols):
                    for j in range(m_cols):
                        term = system.et[i, j]
                        if term is not None:
                            et[i, j] = term.evaluate(f)

                # Add DUT port crosstalk leakage terms
                for i in range(s_cols):
                    for j in range(s_rows):
                        term = system.em[i, j]
                        if term is not None:
                            em[i, j] = term.evaluate(f)

                # Compute the measured values and place them into
                # the output array.
                m = (el + er @ np.linalg.inv(np.identity(s_rows) - s @ em)
                     @ s @ et)
                if self.has_colsys:
                    assert m.shape[1] == 1
                    result[findex, :, sindex] = m[:, 0]
                else:
                    result[findex, :, :] = m

        if ab:
            if self.has_colsys:
                # Construct an arbitrary A matrix consisting of a
                # systems wide vector of 1x1 matrices, giving the
                # reference values for each system (column).
                a = np.empty((frequencies, 1, self.systems),
                             dtype=complex)
                b = np.empty((frequencies, self.m_rows, self.systems),
                             dtype=complex)
                ω0 = 2.0 * math.pi * self.fmax / 10.0
                zeta = 1 / math.sqrt(2.0)
                for findex, f in enumerate(f_vector):
                    jω = 2.0j * math.pi * f
                    d = ω0**2 + 2.0 * zeta * ω0 * jω + jω**2
                    lp = ω0**2 / d
                    for j in range(self.systems):
                        a[findex, 0, j] = lp
                    b[findex, ...] = (result[findex, ...]
                                      @ np.diag(a[findex, 0, :]))

                return a, b

            else:
                # Construct an arbitrary A matrix where the diagonal
                # has a low-pass response and off-diagonals have high
                # pass responses with alternating signs.
                a = np.empty((frequencies, self.m_cols, self.m_cols),
                             dtype=complex)
                b = np.empty((frequencies, self.m_rows, self.m_cols),
                             dtype=complex)
                ω0 = 2.0 * math.pi * self.fmax / 10.0
                zeta = 1 / math.sqrt(2.0)
                for findex, f in enumerate(f_vector):
                    jω = 2.0j * math.pi * f
                    d = ω0**2 + 2.0 * zeta * ω0 * jω + jω**2
                    lp = ω0**2 / d
                    hp = jω**2 / d
                    for i in range(self.m_cols):
                        for j in range(self.m_cols):
                            if i == j:
                                a[findex, i, j] = lp
                            elif (j & 1) == 1:
                                a[findex, i, j] = hp
                            else:
                                a[findex, i, j] = -hp
                    b[findex, ...] = result[findex, ...] @ a[findex, ...]

                return a, b

        return result


def print_matrix(m, file=None, name="m", indent=0, asarray=False):
    (frequencies, rows, cols) = m.shape
    print(f'{" " * indent}{name} = ', end='', file=file)
    if asarray:
        print('np.asarray(', end='', file=file)
    print(f'[', end='', file=file)
    for findex in range(frequencies):
        if findex != 0:
            print(',', end='', file=file)
        print(f'\n{" " * indent}    [', end='', file=file)
        for row in range(rows):
            if row != 0:
                print(f',\n{" " * indent}     ', end='', file=file)
            print('[', end='', file=file)
            for col in range(cols):
                if col != 0:
                    print(', ', end='', file=file)
                v = m[findex, row, col]
                print(f'{v.real:+8.5f}{v.imag:+8.5f}j', end='', file=file)
            print(']', end='', file=file)
        print(']', end='', file=file)
    print(f'\n{" " * indent}]', end='', file=file)
    if asarray:
        print(')', end='', file=file)
    print(file=file)
