"""physics.py — Gordon 모델 물리 코어.

모든 함수는 SI 단위 입출력 (Streamlit 등 UI에 의존하지 않는 순수 함수).

근거 문헌:
  V. Cremers, R. L. Puurunen, J. Dendooven,
  "Conformality in atomic layer deposition: Current status overview of
   analysis and modelling," Appl. Phys. Rev. 6, 021302 (2019).
원 모델:
  R. G. Gordon, D. Hausmann, E. Kim, J. Shepard,
  Chem. Vap. Depos. 9(2), 73-78 (2003).

핵심 가정 (적용 영역):
  - sticking probability s = 1
  - 분자 흐름(molecular flow) 영역
  - 확산 제한(diffusion-limited) 성장
  - 비가역 반응, 반응기 고갈 무시
"""
import math
from units import K_B, N_A


# ===========================================================================
# 1. 기하: (등가)종횡비
# ===========================================================================
def aspect_ratio(geometry, L, w, z=None):
    """Gordon a = L·p/(4A) 에 기반한 (등가)종횡비 EAR.

    geometry : 'circular_hole' | 'square_hole' | 'trench'
               | 'elongated_hole' | 'square_pillar'
    L : 깊이/높이 [m]
    w : 폭/갭 [m]   (pillar: 인접 기둥 사이 갭)
    z : elongated_hole 의 홀 길이 [m]

    반환: a (= EAR, 무차원)

    환산식 (리뷰 논문 II G):
      circular/square hole : a = L/w
      trench               : a = L/(2w)
      elongated_hole       : a = L(w+z)/(2 w z)
      square_pillar        : a = L/(2√2 · w)
                             ※ MC 기반(분자흐름), w/wpillar=3, L/w 5–50 에서 유효
    """
    g = geometry.lower()
    if g in ("circular_hole", "square_hole"):
        return L / w
    if g == "trench":
        return L / (2.0 * w)
    if g == "elongated_hole":
        if z is None:
            raise ValueError("elongated_hole 에는 z(홀 길이)가 필요합니다.")
        return L * (w + z) / (2.0 * w * z)
    if g == "square_pillar":
        return L / (2.0 * math.sqrt(2.0) * w)
    raise ValueError(f"알 수 없는 geometry: {geometry!r}")


def aspect_ratio_simple(L, w):
    """단순 AR = L/w (단면 기준, 3D 기하 무시). 참고 표시용."""
    return L / w


# ===========================================================================
# 2. 노출량 / 침투 깊이 (Gordon 모델 본체)
# ===========================================================================
def flat_saturation_exposure(K_max, m, T):
    """평탄면 포화 노출량 (Pt)_flat [Pa·s].

    (Pt)_flat = K_max · sqrt(2π·m·k_B·T)

    K_max : 포화 면적당 reactant 분자 수 [molecules/m²]
    m     : reactant 분자 1개 질량 [kg]
    T     : 온도 [K]
    """
    return K_max * math.sqrt(2.0 * math.pi * m * K_B * T)


def required_exposure(a, K_max, m, T):
    """리뷰 Eq.(14): 종횡비 a 컨포멀 코팅 필요 노출량 [Pa·s].

    Pt = (Pt)_flat · (1 + 19/4·a + 3/2·a²)
    대(大) a 에서 (3/2)a² 지배 → 노출량 ∝ a².
    """
    flat = flat_saturation_exposure(K_max, m, T)
    return flat * (1.0 + (19.0 / 4.0) * a + (3.0 / 2.0) * a * a)


def penetration_depth(Pt, w, K_max, m, T):
    """리뷰 Eq.(24): 비포화 노출량 Pt 에서 완전 피복 깊이 l [m].

    l = (4w/3) · ( sqrt(1 + 3/8·E*) − 1 ),   E* = Pt / (Pt)_flat

    ⚠️ 계수 4/3 은 **원형 홀 기준 유도**임. 트렌치·square pillar 에 대해서도 동일 식
       (폭 w 사용)을 적용하므로 그 기하에서는 침투 깊이가 근사임.
       (EAR 기반인 required_exposure(Eq.14)는 기하 일반식이라 영향 없음.)

    ※ Eq.(14)와 완전한 역함수가 아님(저~중 AR에서 차이 ≈ 1 + 0.75a).
      따라서 '필요 노출량'은 required_exposure, '침투 깊이'는 본 함수로
      각각 산출하고, 어느 식을 썼는지 결과에 명시할 것.
    """
    flat = flat_saturation_exposure(K_max, m, T)
    E_star = Pt / flat
    return (4.0 * w / 3.0) * (math.sqrt(1.0 + (3.0 / 8.0) * E_star) - 1.0)


def max_aspect_ratio(Pt, K_max, m, T):
    """노출량 예산 Pt 로 컨포멀 코팅 가능한 최대 a (Eq.14 역산).

    1.5·a² + 4.75·a + (1 − E*) = 0 의 양의 근.  E* = Pt/(Pt)_flat.
    E* ≤ 1 이면 평탄면조차 미포화 → 0 반환.
    """
    flat = flat_saturation_exposure(K_max, m, T)
    E_star = Pt / flat
    if E_star <= 1.0:
        return 0.0
    disc = 4.75 ** 2 - 4.0 * 1.5 * (1.0 - E_star)
    a = (-4.75 + math.sqrt(disc)) / (2.0 * 1.5)
    return max(a, 0.0)


def pulse_time(Pt, P):
    """노출량 Pt [Pa·s] 와 reactant 분압 P [Pa] 로부터 펄스 시간 t [s] = Pt/P."""
    return Pt / P


def penetration_depth_vs_time(t, P, w, K_max, m, T):
    """고정 분압 P 에서 feeding(펄스) 시간 t [s] → 완전 피복 침투 깊이 l [m].

    Eq.(24) 를 Pt = P·t 로 치환한 것. (분압 P 가 펄스 동안 일정하다는 가정)
    """
    return penetration_depth(P * t, w, K_max, m, T)


def feeding_time_for_depth(l_target, P, w, K_max, m, T):
    """고정 분압 P 에서 목표 침투 깊이 l_target [m] 도달에 필요한 feeding 시간 t [s].

    Eq.(24) 역산:  E* = (8/3)[(3·l/4w + 1)² − 1],  t = E*·(Pt)_flat / P.
    ⚠️ 반응기 고갈·펄스 rise/fall 무시. 실제 필요한 시간은 이보다 클 수 있음.
    """
    flat = flat_saturation_exposure(K_max, m, T)
    E_star = (8.0 / 3.0) * ((3.0 * l_target / (4.0 * w) + 1.0) ** 2 - 1.0)
    return E_star * flat / P


# ===========================================================================
# 3. 분자 흐름 가정 점검
# ===========================================================================
def mean_free_path(T, P, d):
    """리뷰 Eq.(2): 평균자유행로 λ [m] = k_B·T / (√2·π·d²·P).  d: 분자 지름 [m]."""
    return K_B * T / (math.sqrt(2.0) * math.pi * d * d * P)


def knudsen_number(lmbda, d_p):
    """Kn = λ / d_p.  Kn ≫ 1 (대략 >10) 이면 분자 흐름 영역."""
    return lmbda / d_p


# ===========================================================================
# 4. K_max 보조 산정 (GPC 기반)  — ⚠️ 단순화 근사
# ===========================================================================
def kmax_from_gpc(gpc, rho, M):
    """GPC 기반 K_max 추정 [molecules/m²].

    K_max ≈ GPC · ρ · N_A / M
    gpc : [m/cycle],  rho : 막 밀도 [kg/m³],  M : 막 화학식단위 몰질량 [kg/mol]

    ⚠️ 한계(논문 직접 지적): 반응성 흡착 사이트 수 ≠ 사이클당 증착 원자 수.
       (예: TMA/H2O 에서 OH ~7–9/nm² vs 증착 Al ~4.5/nm²)
       본 추정은 본질적으로 단순화된 근사이므로 결과는 차수 수준으로 해석할 것.
    """
    return gpc * rho * N_A / M
