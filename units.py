"""units.py — 단위 변환 헬퍼.

원칙: 물리 계산은 전부 SI로 수행하고, 변환은 이 모듈에서만 한다.
근거 상수/환산: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019) 및 표준 물리상수.
"""

# ---------------------------------------------------------------------------
# 물리 상수 (SI, CODATA)
# ---------------------------------------------------------------------------
K_B = 1.380649e-23        # 볼츠만 상수 [J/K]
N_A = 6.02214076e23       # 아보가드로 수 [1/mol]

# ---------------------------------------------------------------------------
# 압력  ?  ->  Pa
# ---------------------------------------------------------------------------
TORR_TO_PA = 133.322368421  # 1 Torr = 101325/760 Pa
PRESSURE_TO_PA = {
    "Pa": 1.0,
    "kPa": 1e3,
    "Torr": TORR_TO_PA,
    "mTorr": TORR_TO_PA * 1e-3,
    "mbar": 100.0,
    "bar": 1e5,
}

# ---------------------------------------------------------------------------
# 노출량(exposure)  ?  ->  Pa·s
#   1 Langmuir = 1e-6 Torr·s = 1.33322e-4 Pa·s   (즉 1 Pa·s ≈ 7500 L)
#   ※ 논문 본문의 "1 L = 7500 Pa s" 표기는 혼동을 부르는 서술이며,
#     정확히는 1 Pa·s ≈ 7500 L 이다.
# ---------------------------------------------------------------------------
LANGMUIR_TO_PA_S = 1e-6 * TORR_TO_PA   # ≈ 1.33322e-4
EXPOSURE_TO_PA_S = {
    "Pa·s": 1.0,
    "L": LANGMUIR_TO_PA_S,
    "Langmuir": LANGMUIR_TO_PA_S,
    "Torr·s": TORR_TO_PA,
}

# ---------------------------------------------------------------------------
# 길이  ?  ->  m
# ---------------------------------------------------------------------------
LENGTH_TO_M = {
    "m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6,
    "nm": 1e-9, "A": 1e-10, "Å": 1e-10,
}


# ---------------------------------------------------------------------------
# 변환 함수 (to_* : 단위->SI,  from_* : SI->단위)
# ---------------------------------------------------------------------------
def to_pa(value, unit):       return value * PRESSURE_TO_PA[unit]
def from_pa(value, unit):     return value / PRESSURE_TO_PA[unit]

def to_pa_s(value, unit):     return value * EXPOSURE_TO_PA_S[unit]
def from_pa_s(value, unit):   return value / EXPOSURE_TO_PA_S[unit]

def to_m(value, unit):        return value * LENGTH_TO_M[unit]
def from_m(value, unit):      return value / LENGTH_TO_M[unit]

def celsius_to_kelvin(t_c):   return t_c + 273.15
def kelvin_to_celsius(t_k):   return t_k - 273.15


def molar_mass_to_molecule_mass(M_g_per_mol):
    """분자량 M [g/mol] -> 분자 1개 질량 m [kg]."""
    return (M_g_per_mol / 1000.0) / N_A


def areal_density_per_nm2_to_per_m2(n):  return n * 1e18
def areal_density_per_m2_to_per_nm2(n):  return n * 1e-18
