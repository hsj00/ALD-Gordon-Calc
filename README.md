
# ALD Gordon Model Calculator v5 — 사용자 매뉴얼

> Gordon et al. *Chem. Vap. Deposition* 9, 73 (2003) · Cremers et al. *Appl. Phys. Rev.* 6, 021302 (2019)

---

## 목차

1. [소개](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#1-%EC%86%8C%EA%B0%9C)
2. [설치 및 실행](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#2-%EC%84%A4%EC%B9%98-%EB%B0%8F-%EC%8B%A4%ED%96%89)
3. [핵심 물리 모델](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#3-%ED%95%B5%EC%8B%AC-%EB%AC%BC%EB%A6%AC-%EB%AA%A8%EB%8D%B8)
4. [입력 파라미터 선정 가이드](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#4-%EC%9E%85%EB%A0%A5-%ED%8C%8C%EB%9D%BC%EB%AF%B8%ED%84%B0-%EC%84%A0%EC%A0%95-%EA%B0%80%EC%9D%B4%EB%93%9C)
5. [실전 사용 예시](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#5-%EC%8B%A4%EC%A0%84-%EC%82%AC%EC%9A%A9-%EC%98%88%EC%8B%9C)
6. [결과 해석 방법](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#6-%EA%B2%B0%EA%B3%BC-%ED%95%B4%EC%84%9D-%EB%B0%A9%EB%B2%95)
7. [그래프 탭 활용법](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#7-%EA%B7%B8%EB%9E%98%ED%94%84-%ED%83%AD-%ED%99%9C%EC%9A%A9%EB%B2%95)
8. [Fill Tank 모델 상세](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#8-fill-tank-%EB%AA%A8%EB%8D%B8-%EC%83%81%EC%84%B8)
9. [단위 체계 및 검증](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#9-%EB%8B%A8%EC%9C%84-%EC%B2%B4%EA%B3%84-%EB%B0%8F-%EA%B2%80%EC%A6%9D)
10. [단위 테스트](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#10-%EB%8B%A8%EC%9C%84-%ED%85%8C%EC%8A%A4%ED%8A%B8)
11. [제한사항 및 주의사항](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#11-%EC%A0%9C%ED%95%9C%EC%82%AC%ED%95%AD-%EB%B0%8F-%EC%A3%BC%EC%9D%98%EC%82%AC%ED%95%AD)
12. [참고문헌](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#12-%EC%B0%B8%EA%B3%A0%EB%AC%B8%ED%97%8C)

---

## 1. 소개

이 계산기는 **고종횡비(HAR) 구조에서 ALD precursor의 conformal 코팅에 필요한 최소 노출량과 도징 시간을 예측**하는 도구입니다.

 **주요 기능** :

* Gordon 모델 기반 필요 노출량 및 도징 시간 계산
* Fill Tank ODE 적분을 통한 실제 노출량 산출
* 시간별 침투 깊이(penetration depth) 추적
* 평탄 기판 대비 배수(Dose Multiple) 산출
* 피복률(coverage) 프로파일 θ(z) 시각화
* 5종 프리셋 (DRAM, 3D NAND, FinFET, Flash WL)
* 압력 단위 Torr/mTorr 선택 가능

---

## 2. 설치 및 실행

```bash
pip install streamlit numpy pandas plotly scipy
streamlit run ald_gordon_calculator_v5.py
```

---

## 3. 핵심 물리 모델

### 3.1 Gordon 필요 노출량 (Eq.14)

구조 바닥까지 완전한 conformal 코팅에 필요한 최소 precursor 노출량:

```
E_required = K_max × √(2π·m·k_B·T) × (1/19 + 4a/3 + 2a²)
             ───────────────────────   ──────────────────────
                    scale [Pa·s]           f(a) [무차원]
```

각 항의 역할:

* **K_max** [molecules/m²]: 벽면 1 m²를 포화시키는 데 필요한 분자 수 — "벽면이 얼마나 많이 소비하는가"
* **√(2πmk_BT)** [kg·m/s]: 분자의 운동량 스케일. 분자 충돌률(flux)을 노출량(P·t)으로 변환하는 계수
* **f(a)** : 구조의 기하학적 확산 저항 — "얼마나 깊이 밀어넣어야 하는가"

### 3.2 일반화 종횡비 (EAR)

| 구조                | EAR 공식      | 예시 (L=5μm, w=100nm) |
| ------------------- | ------------- | --------------------- |
| Cylindrical Hole    | a = L/w       | 50                    |
| Infinite Trench     | a = L/(2w)    | 25                    |
| Square Pillar Array | a = L/(2√2·w) | 17.7                  |

같은 AR이라도 Trench는 양면이 열려있으므로 EAR = AR/2. Hole보다 코팅이 쉽습니다.

### 3.3 K_max — Film property와 Precursor property의 구분

K_max는 **"한 ALD half-cycle에서 표면 단위 면적당 흡착되는 precursor 분자의 최대 개수"**입니다.

계산 공식:

```
K_max = GPC [m] × ρ_film [kg/m³] × N_A [/mol] / MW_film [g/mol]
```

이 공식의 논리: "GPC 두께의 film 안에 원자가 몇 개인가?" = "한 cycle에 소비된 precursor 분자 수"

 **Film property vs Precursor property 사용처** :

| 파라미터            | 무엇의 성질?         | 이유                       |
| ------------------- | -------------------- | -------------------------- |
| MW (분자량)         | **Precursor**        | 기체상 확산 속도 결정      |
| d (분자 직경)       | **Precursor**        | Knudsen 수, 평균 자유 경로 |
| s₀ (sticking coeff) | **Precursor + 표면** | 표면 반응 확률             |
| GPC                 | **공정 측정값**      | 한 cycle의 성장 두께       |
| ρ_film, MW_film     | **Film**             | K_max 추정에 사용          |

예: MoO₂Cl₂로 Mo metal 증착 시

* Precursor MW = 198.85 g/mol (MoO₂Cl₂) → 확산 속도 계산에 사용
* Film ρ = 10.2 g/cm³, MW = 95.95 g/mol (Mo) → K_max 계산에 사용

### 3.4 침투 깊이 (Eq.24)

```
l(t) = (4w/3) × [√(1 + (3/8) × E(t)/scale) − 1]
```

* 초기: l ∝ E (선형)
* 후기: l ∝ √E (확산 지배) — 노출량 4배 → 침투 깊이 2배

### 3.5 Knudsen 수

| Kn 범위  | 유동 영역                                    | Gordon 모델 |
| -------- | -------------------------------------------- | ----------- |
| Kn > 10  | Molecular Flow — 분자가 벽에 탁구공처럼 튕김 | 정확        |
| 0.1~10   | Transition — 벽-분자, 분자-분자 충돌 공존    | 근사 적용   |
| Kn < 0.1 | Viscous Flow — 유체 거동                     | 부적합      |

---

## 4. 입력 파라미터 선정 가이드

### 4.1 Depth (L)과 Width (w)

 **측정 방법** : SEM/TEM 단면 분석, CD-SEM (bottom CD 사용 권장)

| 구조                       | Depth (L)        | Width (w)              | AR        |
| -------------------------- | ---------------- | ---------------------- | --------- |
| DRAM Capacitor (DDR4)      | 3~5 μm           | 60~100 nm              | 40~80:1   |
| DRAM Capacitor (DDR5+)     | 5~8 μm           | 40~70 nm               | 80~150:1  |
| 3D NAND 128L Channel Hole  | 8~9 μm           | 100~120 nm             | 70~90:1   |
| 3D NAND 200L+ Channel Hole | 12~15 μm         | 90~110 nm              | 120~160:1 |
| 3D NAND WL Lateral Cavity  | 2~4 μm (lateral) | 15~20 nm (gate height) | 100~200:1 |
| FinFET Gate Trench         | 40~60 nm         | 7~20 nm                | 2~8:1     |

 **3D NAND Word Line 구조 참고** : WL 금속 충전은 수직 slit을 통해 precursor가 진입한 후 **수평 방향의 lateral recess cavity**를 채우는 구조입니다. Cylindrical Hole이 아닌 **Infinite Trench**에 가깝습니다.

* L = slit에서 channel hole까지의 수평 거리 (~2~4 μm)
* w = WL gate height (~15~20 nm, z-pitch ~40 nm의 약 절반)
* 현재 z-pitch: ~40 nm, 차세대: 25~30 nm (gate length 10~15 nm)

### 4.2 Sticking Coefficient (s₀)

 **측정 방법** : Lateral HARS 테스트 구조 + TEM 분석, QCM, Monte Carlo 피팅

| Precursor | Film         | s₀      | 출처                   |
| --------- | ------------ | ------- | ---------------------- |
| TMA       | Al₂O₃        | ~0.01   | Cremers (2019)         |
| TiCl₄     | TiO₂/TiN     | ~0.006  | Cremers (2019)         |
| TEMAHf    | HfO₂         | ~0.1    | Gordon (2003)          |
| DEZ       | ZnO          | ~0.007  | Gordon (2003)          |
| BDEAS     | SiO₂ (PEALD) | ~3×10⁻⁵ | 문헌 추정              |
| MoCl₅     | Mo           | ~0.05   | 문헌 제한적, 실측 필요 |
| MoO₂Cl₂   | Mo           | ~0.04   | 문헌 제한적, 실측 필요 |

### 4.3 GPC

평탄 기판 위 **포화 조건**에서 Ellipsometry 또는 XRR로 측정한 값을 사용합니다.

### 4.4 K_max

자동 추정: `K_max = GPC × ρ_film × N_A / MW_film`. 오차 ±20~50%. 정밀값은 in-situ QCM 측정 권장.

### 4.5 압력 (P)

계산기에서 **Torr 또는 mTorr 중 선택 가능**합니다. 기본값은 Torr입니다. 장비 게이지 단위에 맞춰 선택하세요.

모든 내부 계산은 `Torr_to_Pa()` 함수를 통해 Pa로 변환 후 SI 단위로 수행됩니다. mTorr 입력 시에도 mTorr → Torr → Pa 순서로 함수를 통해 변환됩니다.

### 4.6 Fill Tank 파라미터

| 파라미터        | 결정 방법                                       | 일반 범위         |
| --------------- | ----------------------------------------------- | ----------------- |
| V_fill (cc)     | 장비 매뉴얼 'canister volume' + 배관 사체적     | 10~200 cc         |
| P_before (Torr) | Fill tank 게이지 판독 (충전 후)                 | 1~50 Torr         |
| P_after (Torr)  | Fill tank 게이지 판독 (방출 후)                 | P_before의 50~90% |
| V_chamber (L)   | 장비 매뉴얼 또는 N₂ purge로 실측                | 1~50 L            |
| S_pump (L/s)    | Spec의 50~80% (유효 S). dP/dt 실측 권장         | 50~500 L/s        |
| C_valve (L/s)   | Fast-fill이면 생략. Flow restrictor 있으면 입력 | 10~100 L/s        |

---

## 5. 실전 사용 예시

### 예시 1: DRAM 커패시터 HfO₂

1. 프리셋 "DRAM Capacitor" 선택
2. TEM 결과에 맞춰 L=4.2μm, w=75nm으로 수정
3. 결과: EAR=56, Required ~350 L, Pulse ~0.46s @ 0.1 Torr
4. Fill Tank 확인: V_fill=50cc, ΔP=5T → E_actual=2500 L → 포화도 충분

### 예시 2: 3D NAND Channel Hole Al₂O₃

1. 프리셋 "3D NAND Channel Hole" 선택
2. L=8μm, w=110nm, TMA, 0.2 Torr
3. EAR=72.7, Required ~1,200 L → 도전적
4. Fill Tank 최적화: P_before=20T, V_fill=100cc → 포화도 110%

### 예시 3: Flash WL Mo Fill

1. 프리셋 "Flash WL Mo Fill" 선택
2. 구조: **Infinite Trench** (lateral recess cavity)
3. L=3μm (slit→hole lateral 거리), w=20nm (gate height)
4. EAR=75:1, Dose Multiple ~수백 배 → multi-pulse 전략 검토

---

## 6. 결과 해석 방법

### 6.1 결과 카드 (6개)

| 카드              | 의미                                              |
| ----------------- | ------------------------------------------------- |
| AR                | 기하학적 종횡비                                   |
| EAR               | 일반화 종횡비 (모델 입력)                         |
| Required Exposure | 100% step coverage에 필요한 최소 노출량 [L]       |
| Required Pulse    | 기준 압력에서의 최소 도징 시간 [s]                |
| Knudsen Kn        | 유동 영역 판정 (>10이면 모델 신뢰)                |
| **Dose Multiple** | 평탄 기판 포화 노출량 대비 HAR 필요 노출량의 배수 |

### 6.2 신호등 색상 해석

 **EAR** : 🟢 < 10 (쉬움) / 🟡 10~50 (DRAM급) / 🔴 > 50 (3D NAND급)

 **필요 노출량** : 🟢 < 100 L / 🟡 100~1,000 L / 🔴 > 1,000 L

 **Knudsen** : 🟢 Kn > 10 (molecular) / 🟡 0.1~10 (transition) / 🔴 < 0.1 (viscous)

### 6.3 Dose Multiple — "평탄 기판 대비 100~1000배"의 정체

ALD에서 "평탄 대비 수백 배 노출이 필요하다"는 표현이 자주 등장합니다. 이 계산기의 **Dose Multiple**이 바로 이 비율입니다.

```
Dose Multiple = E_required(HAR) / E_flat(평탄 기판)
```

 **문헌 근거** : Cremers (2019)에 따르면, HfO₂ ALD에서 평탄 기판 포화에 3~43 L이 필요한 반면, AR≈43 홀에서는 9,000 L이 필요하여 약 200~3,000배에 달합니다.

 **핵심 구분 — 이 두 비율을 혼동하지 마세요** :

| 비율                 | 정의                  | 의미                                                        |
| -------------------- | --------------------- | ----------------------------------------------------------- |
| **Dose Multiple**    | E_req(HAR) / E_flat   | "평탄 대비 몇 배 어려운가" →**Gordon Required에 이미 포함** |
| **Saturation Ratio** | E_actual / E_req(HAR) | "현재 공급량이 충분한가" →**≥100%이면 conformal 달성**      |

Dose Multiple이 500×라는 것은 "이 구조는 평탄 기판보다 500배 어렵다"는 뜻이지, Required Exposure에 다시 500을 곱해야 한다는 뜻이 아닙니다. Required Exposure에 이미 이 배수가 포함되어 있습니다.

### 6.4 포화도 (Saturation Ratio) 해석

| 포화도 | 코팅 상태             | 조치                     |
| ------ | --------------------- | ------------------------ |
| ≥ 100% | 바닥까지 conformal    | 현 조건 유지             |
| 90~99% | 하단 약간 얇아짐      | 도징 시간 소폭 증가      |
| 70~89% | 하단 상당 부분 미코팅 | 시간/압력 1.2~1.5배 증가 |
| < 50%  | 유의미한 코팅 불가    | 근본적 조건 변경 필요    |

### 6.5 Fill Tank Model A vs B

| Model            | 정의                     | 의미              |
| ---------------- | ------------------------ | ----------------- |
| A: P_eq × t_dose | 상한 (P가 일정하다 가정) | 낙관적            |
| B: ODE 수치적분  | 실제 압력 변화 반영      | **현실적 (권장)** |

차이 > 20%이면: t_dose < τ → 도징 시간 연장 효과적 / t_dose > τ → ΔP 또는 V_fill 증대 필요

---

## 7. 그래프 탭 활용법

| 탭                       | 답해주는 질문                                            |
| ------------------------ | -------------------------------------------------------- |
| Exposure vs EAR          | 구조 shrink 시 노출량이 얼마나 증가하나? (EAR² 스케일링) |
| Pulse Time vs P          | 어떤 압력에서 목표 도징 시간을 달성하나?                 |
| Pulse Time vs T          | 공정 온도 변경이 도징 시간에 미치는 영향?                |
| EAR Structure Comparison | 같은 AR의 Hole vs Trench 코팅 난이도 차이?               |

---

## 8. Fill Tank 모델 상세

### ODE 수식

```
dP_f/dt = −C(P_f − P_c) / V_f        (fill tank 압력 감소)
dP_c/dt =  C(P_f − P_c) / V_c − S·P_c / V_c   (챔버 압력)
```

### 핵심 파라미터

* **P_eq = ΔP × V_fill / V_chamber** : 평형 챔버 압력
* **τ = V_chamber / S_pump** : 시상수
* **E_max = P_eq × τ** : 달성 가능 최대 노출량

### 판단 기준

| 조건               | 대응                              |
| ------------------ | --------------------------------- |
| E_max > E_required | t_dose ≥ t_full로 설정            |
| E_max < E_required | ΔP↑, V_fill↑, 또는 S_pump↓        |
| t_dose/τ < 0.5     | 도징 시간 연장이 효과적           |
| t_dose/τ > 3       | 추가 시간 효과 미미, ΔP 증대 필요 |

---

## 9. 단위 체계 및 검증

### 9.1 압력 단위 처리

계산기는 **Torr와 mTorr 중 선택 가능**합니다 (사이드바의 "압력 단위" 라디오 버튼).

내부 변환 경로:

```
Torr 모드:  사용자 입력 [Torr]  → P_Torr          → Torr_to_Pa() → P_Pa [Pa]
mTorr 모드: 사용자 입력 [mTorr] → ÷1000 → P_Torr  → Torr_to_Pa() → P_Pa [Pa]
```

**모든 내부 계산은 SI 단위(Pa)로 수행**됩니다. 하드코딩된 변환 상수(`0.133322`, `/133.322`)는 코드 전체에서 제거되었으며, 변환은 오직 `Torr_to_Pa()`, `Pa_to_Torr()`, `Pa_s_to_L()` 함수만을 통해 이루어집니다.

### 9.2 전체 단위 흐름

```
입력 단위     SI 변환              계산                출력 단위
──────────   ──────────────────   ────────────────   ─────────
μm           → m  (×10⁻⁶)        gordon_a → a       [무차원]
nm           → m  (×10⁻⁹)
°C           → K  (+273.15)
Torr/mTorr   → Pa (Torr_to_Pa)   E = scale × f(a)   Pa·s → L
g/mol        → kg (×10⁻³/N_A)    t = E / P          s
Å            → m  (×10⁻¹⁰)       λ = kBT/(..d²P)    m → nm
nm, g/cm³    → K_max [1/m²]      Kn = λ/w           [무차원]
cc           → m³ (×10⁻⁶)        Fill Tank ODE
L            → m³ (×10⁻³)          P_c(t)            Pa → Torr
L/s          → m³/s (×10⁻³)        E_cum(t)          Pa·s → L
```

### 9.3 검증 방법

```bash
python run_tests.py     # 44개 자동 검증 테스트
```

검증 항목:

* 물리 상수 (k_B, N_A) CODATA 2018 정확값 확인
* 단위 변환 함수 왕복 검증 (Torr↔Pa, Pa·s↔L)
* K_max 차원 분석: [m]×[g/m³]×[/mol]/[g/mol] = [molecules/m²]
* Required Exposure: 두 가지 독립 경로로 계산 → 완전 일치
* Torr vs mTorr 입력 시 동일 결과 확인
* Fill Tank ODE vs 해석해: 상대 오차 < 2%

---

## 10. 단위 테스트

```bash
python run_tests.py          # Standalone (streamlit 불필요)
pytest test_gordon_calculator.py -v   # pytest 사용 시
```

44개 테스트 항목:

* 단위 변환 정확도 (3개)
* Gordon EAR 계산 (4개)
* K_max 추정 (1개)
* 평균 자유 경로 (2개)
* Fill Tank ODE vs 해석해 (22개)
* 침투 깊이 (3개)
* find_t_full (4개)
* 입력 검증 (3개)
* Gordon 노출량 공식 (2개)

---

## 11. 제한사항 및 주의사항

### 모델 한계

1. **s₀ = 1 가정** : Gordon 모델은 완전 diffusion-limited 경우를 계산. s₀ < 1이면 실제 프로파일이 더 완만하므로, 노출량은 **보수적 상한**
2. **1차 Langmuir 흡착** : s(θ) = s₀(1−θ)만 고려. 복잡한 반응 메커니즘 미반영
3. **등온 가정** : 구조 내부 온도 구배 미고려
4. **이상 기체 가정** : 저압(< 수 Torr)에서 유효

### 실무 주의사항

1. **Kn < 1** : Gordon 모델 부적합. CFD 시뮬레이션 권장
2. **K_max 오차** : GPC 기반 추정은 ±20~50%. 정밀 공정은 QCM 실측 필요
3. **Precursor 열분해** : 고온에서 CVD 성분 혼입 가능. Gordon 모델은 순수 ALD만 고려
4. **표면 화학 불균일** : 구조 내부의 표면 화학이 깊이에 따라 다르면 s₀가 위치 의존 → 정확도 저하

---

## 12. 참고문헌

1. R. G. Gordon et al., "A Kinetic Model for Step Coverage by ALD in Narrow Holes or Trenches," *Chem. Vap. Deposition*  **9** , 73–78 (2003)
2. V. Cremers et al., "Conformality in ALD: Current Status Overview of Analysis and Modelling," *Appl. Phys. Rev.*  **6** , 021302 (2019)
3. H. C. M. Knoops et al., "Conformality of Plasma-Assisted ALD: Physical Processes and Modeling," *J. Electrochem. Soc.*  **157** , G241–G249 (2010)
4. S. M. Sze, M. K. Lee,  *Semiconductor Devices: Physics and Technology* , 3rd Ed., Wiley (2012)

---

*ALD Gordon Model Calculator v5 — 2025*
