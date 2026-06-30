"""multicycle.py — 다중 사이클 EAR 변화 (확장 모델).

배경(논문 Sec V C 6; Gordon HfO2 홀 EAR 36→43; Perez et al.):
  매 ALD 사이클마다 벽에 GPC 두께가 쌓여 피처 폭이 줄고 EAR 이 증가한다.
  고정(비포화) 노출에서는 사이클이 진행될수록 침투 깊이가 감소
  (→ 경사진 두께 프로파일). Gordon 모델을 좁아지는 피처에 반복 적용한다.

가정/한계:
  - 벽 성장은 등각(conformal)·일정 GPC 가정. 실제 GPC 는 깊이/사이클에 따라 변할 수 있음.
  - 폭만 감소(w_n = w0 − walls·GPC·n), 깊이 L 은 일정으로 근사(L ≫ w).
  - 침투 깊이는 Eq.(24)(s=1, 분자 흐름) 기반. s<1 효과는 별도(Monte Carlo) 모듈.
"""
import numpy as np
import physics as phys


def pore_width(w0, gpc, n, walls=2):
    """n 사이클 후 폭 [m]. walls=2: 마주보는 두 벽이 각각 GPC 만큼 성장."""
    return w0 - walls * gpc * n


def cycles_to_close(w0, gpc, walls=2):
    """폭이 0 이 되는(클로깅) 사이클 수."""
    return w0 / (walls * gpc)


def evolution(geom, L, w0, gpc, n_max, z=None, walls=2):
    """사이클별 (n, w_n[m], EAR_n, film_thickness[m]). 폭 > 0 구간만 반환."""
    n = np.arange(0, int(n_max) + 1)
    w_n = pore_width(w0, gpc, n, walls)
    valid = w_n > 0
    n, w_n = n[valid], w_n[valid]
    ear = np.array([phys.aspect_ratio(geom, L, w, z) for w in w_n])
    film = gpc * n
    return n, w_n, ear, film


def penetration_per_cycle(geom, L, w0, gpc, Pt_cycle, K_max, m, T,
                          n_max, z=None, walls=2):
    """고정 per-cycle 노출 Pt_cycle [Pa·s] 에서 사이클별 침투 깊이 l_n [m].

    Eq.(24) 를 좁아지는 폭 w_n 에 반복 적용. l_n ∝ w_n 이므로 사이클이 진행될수록
    절대 침투 깊이는 감소(Perez/Gordon 과 일치).
    """
    n, w_n, ear, film = evolution(geom, L, w0, gpc, n_max, z, walls)
    l = np.array([phys.penetration_depth(Pt_cycle, w, K_max, m, T) for w in w_n])
    return n, l


# ---------------------------------------------------------------------------
# 자가검증: Gordon 홀 EAR 36→43 재현
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    L = 7200e-9      # 7.2 µm
    w0 = 200e-9      # 0.20 µm  → EAR0 = 36
    gpc = 0.10e-9    # 0.10 nm/cycle
    N = 163          # Δw = 2·gpc·N = 32.6 nm → w_N = 167.4 nm → EAR ≈ 43

    n, w_n, ear, film = evolution("circular_hole", L, w0, gpc, N)
    print(f"EAR0 = {ear[0]:.1f}  (기대 36)")
    print(f"EAR@N={n[-1]} = {ear[-1]:.1f}  (기대 ≈43)")
    print(f"폭: {w_n[0]*1e9:.1f} → {w_n[-1]*1e9:.1f} nm,  막두께 {film[-1]*1e9:.1f} nm")
    print(f"클로깅 사이클 = {cycles_to_close(w0, gpc):,.0f}")
    ok = abs(ear[0] - 36) < 0.5 and abs(ear[-1] - 43) < 0.5
    print("자가검증:", "PASS" if ok else "FAIL")
