"""validation.py — 논문 레퍼런스 케이스 대조 검증.

실행:  python validation.py

목적: K_max 가정의 불확실성 때문에 절대값 정밀 일치는 기대하지 않는다.
      '차수(order of magnitude)·추세·자기일관성' 검증이 목표.
출처: Cremers et al., Appl. Phys. Rev. 6, 021302 (2019)
      (HfO2: 평탄면 3–43 L, 홀 EAR≈36–43 → ~9000 L, 홀/트렌치 노출비 ≈4 등)
"""
import math
import physics as phys
import units as u


def kmax_from_gpc_practical(gpc_nm, rho_g_cm3, M_g_mol):
    """편의 래퍼: GPC[nm], ρ[g/cm³], M[g/mol] -> K_max[/m²]."""
    return phys.kmax_from_gpc(gpc_nm * 1e-9, rho_g_cm3 * 1000.0, M_g_mol / 1000.0)


results = []
def check(name, value, lo, hi, unit=""):
    ok = lo <= value <= hi
    results.append((name, value, lo, hi, unit, ok))
    return ok


# --- 공통 입력: Hf(NMe2)4 / HfO2, 200°C ------------------------------------
T = u.celsius_to_kelvin(200.0)
M_HfNMe2 = 354.8                                  # g/mol (Hf 178.49 + 4×NMe2 ≈44.08)
m = u.molar_mass_to_molecule_mass(M_HfNMe2)
# HfO2: GPC≈0.10 nm/cycle, ρ≈9.68 g/cm³, M(HfO2)=210.49 g/mol
K_max = kmax_from_gpc_practical(0.10, 9.68, 210.49)

# --- Case 1: 평탄면 포화 노출량 → 논문 3–43 L -------------------------------
flat_pa_s = phys.flat_saturation_exposure(K_max, m, T)
flat_L = u.from_pa_s(flat_pa_s, "L")
check("C1 평탄면 (Pt)_flat", flat_L, 1.0, 60.0, "L")          # 3–43 L 밴드(여유 포함)

# --- Case 2: HfO2 홀 EAR≈36–43 필요 노출량 → 논문 ~9000 L -------------------
for a in (36.0, 43.0):
    Pt = phys.required_exposure(a, K_max, m, T)
    check(f"C2 홀 EAR={a:.0f} 필요노출", u.from_pa_s(Pt, "L"), 3e3, 2e4, "L")

# --- Case 3: 같은 L/w 에서 홀/트렌치 노출비 → 고AR에서 ≈4 -------------------
Lw, w = 50.0, 1.0
a_hole   = phys.aspect_ratio("circular_hole", Lw, w)
a_trench = phys.aspect_ratio("trench", Lw, w)
ratio = phys.required_exposure(a_hole, K_max, m, T) / phys.required_exposure(a_trench, K_max, m, T)
check("C3 홀/트렌치 노출비(L/w=50)", ratio, 3.6, 4.0, "x")

# --- Case 4: 기하 EAR 환산식 직접 검증 -------------------------------------
check("C4 trench=(L/w)/2", phys.aspect_ratio("trench", 100, 1), 49.99, 50.01)
check("C4 elong(z=w)=L/w", phys.aspect_ratio("elongated_hole", 100, 2, z=2), 49.99, 50.01)
exp_pillar = 100 / (2 * math.sqrt(2) * 1)
check("C4 pillar=L/(2v2 w)", phys.aspect_ratio("square_pillar", 100, 1),
      exp_pillar - 1e-6, exp_pillar + 1e-6)

# --- Case 5: 정방향↔역방향 자기일관성 (required ↔ max_AR 왕복) --------------
a_in = 20.0
Pt5 = phys.required_exposure(a_in, K_max, m, T)
a_out = phys.max_aspect_ratio(Pt5, K_max, m, T)
check("C5 정/역 왕복 a", a_out, a_in - 1e-3, a_in + 1e-3)

# --- Case 6: 펄스 시간 = Pt/P ---------------------------------------------
check("C6 펄스시간 Pt/P", phys.pulse_time(0.4, 0.2), 1.999, 2.001, "s")

# --- Case 7: 분자 흐름 점검 (저압 펌프형, dp=100nm 이면 Kn≫1) --------------
lmbda = phys.mean_free_path(T, P=0.27, d=7e-10)   # P=0.27 Pa, d≈7 Å
Kn = phys.knudsen_number(lmbda, d_p=100e-9)
check("C7 Kn (저압,dp=100nm)", Kn, 10.0, 1e12)     # 분자 흐름이면 Kn>10


# ---------------------------------------------------------------------------
# 출력
# ---------------------------------------------------------------------------
print(f"{'검증 항목':<26}{'계산값':>13}{'허용범위':>26}  판정")
print("-" * 74)
allok = True
for name, val, lo, hi, unit, ok in results:
    allok &= ok
    rng = f"[{lo:.3g}, {hi:.3g}]"
    print(f"{name:<26}{val:>11.3g} {unit:<2}{rng:>24}  {'PASS' if ok else 'FAIL'}")
print("-" * 74)
print(f"K_max = {K_max:.3e} /m²  = {u.areal_density_per_m2_to_per_nm2(K_max):.2f} /nm²")
print(f"(Pt)_flat = {flat_L:.2f} L,   λ(0.27Pa) = {lmbda*1e3:.2f} mm,   Kn = {Kn:.3g}")
print("전체:", "ALL PASS" if allok else "일부 FAIL")
