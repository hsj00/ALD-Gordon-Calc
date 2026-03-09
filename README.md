
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
9. [단위 테스트 실행](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#9-%EB%8B%A8%EC%9C%84-%ED%85%8C%EC%8A%A4%ED%8A%B8-%EC%8B%A4%ED%96%89)
10. [제한사항 및 주의사항](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#10-%EC%A0%9C%ED%95%9C%EC%82%AC%ED%95%AD-%EB%B0%8F-%EC%A3%BC%EC%9D%98%EC%82%AC%ED%95%AD)
11. [참고문헌](https://claude.ai/chat/4e97662f-04e1-4206-806c-df6c605468a1#11-%EC%B0%B8%EA%B3%A0%EB%AC%B8%ED%97%8C)

---

## 1. 소개

이 계산기는 **고종횡비(High Aspect Ratio) 구조에서 ALD precursor의 conformal 코팅에 필요한 최소 노출량(exposure)과 도징 시간(pulse time)을 예측**하는 도구입니다.

 **해결하는 문제** : "AR=50:1 실린더 홀에 TMA로 Al₂O₃를 코팅할 때, 바닥까지 100% step coverage를 달성하려면 몇 Langmuir의 노출량이 필요한가?"

 **주요 기능** :

* Gordon 모델 기반 필요 노출량 및 도징 시간 계산
* Fill Tank ODE 적분을 통한 실제 노출량 산출
* 시간별 침투 깊이(penetration depth) 추적
* 피복률(coverage) 프로파일 θ(z) 시각화
* 다양한 구조 형상(Hole, Trench, Pillar) 비교

 **대상 사용자** : ALD 공정 엔지니어, 반도체 소자 공정 개발자, 박막 연구자

---

## 2. 설치 및 실행

### 필수 패키지

```bash
pip install streamlit numpy pandas plotly scipy
```

### 실행

```bash
streamlit run ald_gordon_calculator_v5.py
```

브라우저에서 `http://localhost:8501`로 접속합니다.

### 단위 테스트

```bash
python run_tests.py
```

또는 pytest가 설치된 경우:

```bash
pytest test_gordon_calculator.py -v
```

---

## 3. 핵심 물리 모델

### 3.1 Gordon 필요 노출량 (Eq.14)

구조 바닥까지 완전한 conformal 코팅을 달성하기 위한 최소 precursor 노출량:

```
E_required = K_max × √(2π·m·k_B·T) × (1/19 + 4a/3 + 2a²)
```

여기서 `a`는 일반화 종횡비(EAR, Equivalent Aspect Ratio)입니다.

### 3.2 일반화 종횡비 (EAR)

동일한 AR이라도 구조 형상에 따라 코팅 난이도가 다릅니다:

| 구조                | EAR 공식         | 예시 (L=5μm, w=100nm) |
| ------------------- | ---------------- | ---------------------- |
| Cylindrical Hole    | a = L/w          | 50                     |
| Square Hole         | a = L/w          | 50                     |
| Infinite Trench     | a = L/(2w)       | 25                     |
| Elongated Hole      | a = L(w+z)/(2wz) | 구조 의존              |
| Square Pillar Array | a = L/(2√2·w)  | 17.7                   |

 **핵심** : Trench는 양면이 열려있으므로 EAR = AR/2. 같은 AR이라도 Trench가 Hole보다 코팅이 쉽습니다.

### 3.3 침투 깊이 (Eq.24)

주어진 시점에서 precursor가 도달한 구조 내부 깊이:

```
l(t) = (4w/3) × [√(1 + (3/8) × E(t)/scale) − 1]
```

* 초기(E 작을 때): l ∝ E (선형 증가)
* 후기(E 클 때): l ∝ √E (확산 지배)
* 노출량이 4배가 되어야 침투 깊이가 2배가 됩니다

### 3.4 Knudsen 수

```
Kn = λ/w,   λ = k_B·T / (√2·π·d²·P)
```

| Kn 범위       | 유동 영역      | 물리적 의미                          | Gordon 모델 |
| ------------- | -------------- | ------------------------------------ | ----------- |
| Kn > 10       | Molecular Flow | 분자가 벽에 부딪히며 탁구공처럼 이동 | 정확        |
| 0.1 < Kn < 10 | Transition     | 분자-벽, 분자-분자 충돌 공존         | 근사 적용   |
| Kn < 0.1      | Viscous Flow   | 유체처럼 거동                        | 부적합      |

---

## 4. 입력 파라미터 선정 가이드

### 4.1 구조 파라미터: Depth (L)과 Width (w)

**어디서 가져오는가?**

* 설계 스펙(Design Rule)에서 목표 L과 w를 확인
* 실제 에칭 결과를 TEM 또는 SEM 단면 분석으로 측정
* CD-SEM으로 top/bottom CD 측정 (w = bottom CD 사용 권장)

 **대표값 참고표** :

| 구조                       | 대표 Depth (L) | 대표 Width (w) | AR        |
| -------------------------- | -------------- | -------------- | --------- |
| DRAM Capacitor (DDR4)      | 3~5 μm        | 60~100 nm      | 40~80:1   |
| DRAM Capacitor (DDR5+)     | 5~8 μm        | 40~70 nm       | 80~150:1  |
| 3D NAND 96L Channel Hole   | 6~7 μm        | 100~140 nm     | 50~70:1   |
| 3D NAND 128L Channel Hole  | 8~9 μm        | 100~120 nm     | 70~90:1   |
| 3D NAND 200L+ Channel Hole | 12~15 μm      | 90~110 nm      | 120~160:1 |
| 3D NAND Slit (Trench)      | 6~10 μm       | 100~200 nm     | 30~100:1  |
| FinFET Gate Trench         | 40~60 nm       | 7~20 nm        | 2~8:1     |
| GAA Nanosheet Gap          | 8~12 nm        | 8~12 nm        | ~1:1      |

 **실전 팁** : w는 구조 내부에서 가장 좁은 지점(bottleneck)을 사용하세요. 경사가 있는 에칭 프로파일의 경우 top CD가 아닌 minimum CD를 입력해야 가장 보수적인(안전한) 결과를 얻습니다.

### 4.2 Sticking Coefficient (s₀)

 **물리적 의미** : Precursor 분자가 표면에 충돌했을 때 반응하는 확률. 0에서 1 사이의 값.

 **측정 방법** :

1. **Lateral HARS 구조** : 알려진 AR의 테스트 구조에 ALD 수행 후, step coverage 프로파일을 TEM으로 측정하여 Gordon 모델로 역산
2. **QCM (Quartz Crystal Microbalance)** : In-situ로 표면 반응률을 직접 측정
3. **Monte Carlo 시뮬레이션 피팅** : 실험 프로파일에 맞추어 s₀ 추출

**문헌 대표값** (Gordon 2003, Cremers 2019 기반):

| Precursor            | Film          | s₀        | 출처                   |
| -------------------- | ------------- | ---------- | ---------------------- |
| TMA [Al(CH₃)₃]     | Al₂O₃       | ~0.01      | Cremers (2019)         |
| TiCl₄               | TiO₂/TiN     | ~0.006     | Cremers (2019)         |
| TEMAHf [Hf(NEtMe)₄] | HfO₂         | ~0.1       | Gordon (2003)          |
| DEZ [Zn(C₂H₅)₂]   | ZnO           | ~0.007     | Gordon (2003)          |
| BDEAS                | SiO₂ (PEALD) | ~3×10⁻⁵ | 문헌 추정              |
| MoCl₅               | Mo            | ~0.05      | 문헌 제한적, 실측 필요 |
| MoO₂Cl₂            | Mo            | ~0.04      | 문헌 제한적, 실측 필요 |

 **주의** : s₀는 공정 온도, 표면 상태, co-reactant 종류에 크게 의존합니다. 정확한 공정 최적화를 위해서는 실측이 필수입니다. 문헌값은 초기 설계 단계의 추정용으로만 사용하세요.

### 4.3 GPC (Growth Per Cycle)

**어디서 가져오는가?**

* **평탄 기판** 위에서 포화 조건(saturation regime)으로 ALD를 수행한 뒤 측정
* 측정 방법: Ellipsometry(가장 일반적), XRR(X-ray Reflectometry), TEM 단면
* 반드시 포화 곡선(GPC vs dose time)을 확인하여 포화 영역의 GPC를 사용

 **대표값** :

| Film                              | 대표 GPC (nm/cycle) |
| --------------------------------- | ------------------- |
| Al₂O₃ (TMA+H₂O)                | 0.08~0.12           |
| TiO₂ (TiCl₄+H₂O)               | 0.04~0.06           |
| HfO₂ (TEMAHf+H₂O)               | 0.08~0.12           |
| ZnO (DEZ+H₂O)                    | 0.15~0.22           |
| SiO₂ (BDEAS+O₂ plasma)          | 0.10~0.15           |
| Mo (MoCl₅+H₂, or MoO₂Cl₂+H₂) | 0.03~0.06           |

### 4.4 K_max (최대 표면 흡착량)

**자동 추정 방법** (Film DB):

```
K_max = GPC [m] × ρ_film [kg/m³] × N_A / MW_film [g/mol]
```

 **오차 범위** : ±20~50%. 실제 K_max는 표면 -OH 밀도, 리간드 차폐(steric hindrance), 증착 온도에 따라 변합니다.

 **보다 정확한 방법** : In-situ QCM으로 한 cycle 동안의 질량 증가(Δm)를 측정하고, K_max = Δm × N_A / MW_film으로 직접 계산.

### 4.5 공정 온도 (T)와 기준 압력 (P)

 **온도** :

* ALD window 내의 증착 온도를 입력합니다
* 모델에서 온도는 분자 열속도(∝√T)에만 영향을 미치며, 반응 활성화 에너지는 고려하지 않습니다
* 일반적으로 150~400°C 범위

 **기준 압력** :

* "일정 압력 P×t" 모드에서 사용되는 precursor의 도징 중 챔버 압력
* 공정 레시피의 "dose pressure" 또는 Chamber Manometer 판독값 사용
* 단위: mTorr (1 mTorr = 0.133 Pa)
* 일반적으로 10~1000 mTorr 범위

### 4.6 Fill Tank 파라미터

Fill Tank 모드를 사용할 때 필요한 파라미터와 각각의 결정 방법:

#### Fill Tank 부피 (V_fill)

* **어디서 확인** : 장비 매뉴얼의 "canister volume", "precursor cylinder headspace", 또는 "manifold volume"
* **주의** : 밸브~챔버 사이의 배관 사체적(dead volume)도 포함해야 정확
* **일반 범위** : 10~200 cc
* **측정법** : 알려진 압력의 기체를 fill tank에 채운 뒤 밸브를 열어 팽창시키고 PV=nRT로 역산

#### Fill 전/후 압력 (P_before, P_after)

* **P_before** : Fill tank에 precursor 증기를 채운 뒤의 압력. Canister 온도에서의 precursor 증기압으로 결정됨
* **P_after** : Precursor를 챔버로 방출한 뒤 fill tank에 남은 압력
* **측정** : Fill tank에 부착된 게이지(Baratron 등)로 직접 판독
* **실전 예시** : TMA at 25°C → 증기압 ~11 Torr. P_before=10 Torr, P_after=5 Torr이면 ΔP=5 Torr

#### Pump Speed (S)

* **주의** : Pump spec sheet의 값은 **무부하(no-load) 펌핑 속도**입니다. 실제 챔버에서의 유효 S는 배관 conductance로 인해 **50~80%로 감소**합니다
* **측정법** : 챔버를 일정 P에서 안정화 → 가스 공급 차단 → 압력 감소 곡선 기록 → S_eff = V_c × (dP/dt) / P
* **일반 범위** : 50~500 L/s (Turbo pump), 10~50 L/s (Dry pump)

#### Valve Conductance (C)

* **입력 안 해도 되는 경우 (Fast-fill 조건)** : ALD 밸브가 완전 개방(fully open)이고 C >> S일 때. 이 경우 fill tank의 precursor가 거의 즉시 챔버와 평형을 이룸
* **입력이 필요한 경우** : 밸브가 partially open이거나, flow restrictor가 있거나, 도징 중 fill tank 압력이 서서히 변하는 것이 관측될 때
* **판단 기준** : 도징 중 fill tank 압력이 빠르게(수십 ms 이내) P_after에 도달하면 fast-fill → C 입력 불필요
* **대표값** : Swagelok ALD valve ~20-50 L/s (N₂ 기준)

---

## 5. 실전 사용 예시

### 예시 1: DRAM 커패시터 HfO₂ 코팅

 **상황** : DRAM 커패시터에 HfO₂ dielectric을 ALD로 코팅. 구조는 실린더 홀, 깊이 4 μm, 폭 80 nm.

 **Step 1 — 프리셋 사용** :

* 사이드바에서 "DRAM Capacitor (typical)" 프리셋 선택
* 값이 자동 입력됨: L=4μm, w=80nm, TEMAHf, HfO₂, T=250°C, P=100mTorr

 **Step 2 — 실제 구조에 맞게 수정** :

* TEM 분석 결과 실제 bottom CD = 75 nm → w를 75로 변경
* 실제 에칭 깊이 = 4.2 μm → L을 4.2로 변경

 **Step 3 — 결과 확인** :

* AR = 56:1, EAR = 56:1 (Cylindrical Hole이므로 EAR = AR)
* 필요 노출량 = 약 350 L
* 100 mTorr에서 필요 펄스 시간 = 약 0.46 s
* Kn = 8.2 → Molecular Flow → Gordon 모델 유효

 **Step 4 — Fill Tank으로 실제 노출량 확인** :

* Fill Tank 모드 선택
* V_fill=50cc, P_before=10 Torr, P_after=5 Torr, V_chamber=10 L
* t_dose=0.5 s, S_pump=100 L/s
* 포화도 = 85% → 부족! → t_dose를 0.7 s로 늘림 → 포화도 = 98%

### 예시 2: 3D NAND Al₂O₃ Blocking Oxide

 **상황** : 128L 3D NAND의 channel hole에 Al₂O₃ blocking oxide 코팅. 구조는 실린더 홀, 깊이 8 μm, 폭 110 nm.

 **Step 1** : "3D NAND Channel Hole" 프리셋 선택 후 값 수정

 **Step 2** :

* L=8μm, w=110nm, TMA, Al₂O₃
* s₀=0.01 (Cremers 2019), GPC=0.10 nm/cycle
* T=300°C, P=200 mTorr

 **Step 3 — 결과 판독** :

* EAR = 72.7:1
* 필요 노출량 = 약 1,200 L → 🔴 극도로 어려움
* Kn = 3.8 → 전환 영역 → 근사 적용 가능하나 보정 필요

 **Step 4 — Fill Tank 최적화** :

* 기존 조건으로 포화도 = 55% → 크게 부족
* 해결 전략: ΔP를 높이기 위해 fill tank 온도 상승(P_before↑) 또는 V_fill↑
* P_before=20 Torr, V_fill=100 cc로 변경 → 포화도 = 110% → 충분

### 예시 3: Flash Memory Mo Word Line

 **상황** : 3D NAND Flash의 word line에 Mo 금속을 MoO₂Cl₂ + H₂ ALD로 코팅.

 **Step 1** : "Flash WL Mo Fill" 프리셋 선택

 **Step 2** :

* L=6μm, w=50nm → AR=120:1
* MoO₂Cl₂ (MW=198.85), s₀=0.04, GPC=0.04
* T=500°C, P=300 mTorr

 **Step 3** :

* EAR = 120:1 → 🔴 극한 구조
* 필요 노출량 = 약 15,000 L
* 이 수준에서는 multi-pulse dosing 전략이 필수

---

## 6. 결과 해석 방법

### 6.1 결과 카드 읽는 법

계산 후 화면 상단에 5개의 메트릭 카드가 표시됩니다:

| 카드              | 의미                                    | 좋은 값 방향          |
| ----------------- | --------------------------------------- | --------------------- |
| AR                | 기하학적 종횡비                         | 작을수록 코팅 용이    |
| EAR               | 일반화 종횡비 (모델 입력)               | 작을수록 코팅 용이    |
| Required Exposure | 100% step coverage에 필요한 최소 노출량 | 작을수록 공정 용이    |
| Required Pulse    | 기준 압력에서의 최소 도징 시간          | 짧을수록 throughput↑ |
| Knudsen Kn        | 유동 영역 판정                          | 10 이상이면 모델 신뢰 |

### 6.2 신호등 색상 해석

각 결과 아래에 🟢🟡🔴 신호등과 함께 1줄 해석이 표시됩니다:

 **EAR 해석** :

* 🟢 EAR < 10: 대부분의 ALD로 쉽게 코팅
* 🟡 EAR 10~50: 충분한 노출량 필요 (DRAM급)
* 🔴 EAR 50~100: 매우 도전적 (최신 3D NAND급)
* 🔴 EAR > 100: 극한 구조, precursor 재설계 검토 필요

 **필요 노출량 해석** :

* 🟢 < 100 L: 표준 ALD 공정에서 일반적
* 🟡 100~1,000 L: 도전적, 도징 시간 연장 또는 압력 증대 필요
* 🔴 > 1,000 L: 극도로 어려움, multi-pulse 또는 공정 전략 변경 필요

 **Knudsen 수 해석** :

* 🟢 Kn > 10: 분자-벽 충돌 지배. 분자가 벽을 따라 탁구공처럼 튕기며 이동. Gordon 모델 정확
* 🟡 Kn 0.1~10: 전환 영역. 근사적으로 적용 가능
* 🔴 Kn < 0.1: 점성 유동. 분자가 유체처럼 거동. Navier-Stokes 기반 모델 필요

### 6.3 Dose Multiple (평탄 기판 대비 배수) — v5.1 신규

**이것은 무엇인가?**

ALD 공정에서 "평탄 기판 대비 100~1000배의 노출량이 필요하다"는 말을 자주 듣습니다. 이 계산기의 **Dose Multiple**은 정확히 이 비율을 정량화한 것입니다.

```
Dose Multiple = E_required(HAR 구조) / E_flat(평탄 기판 포화)
```

여기서:

* `E_flat = K_max × √(2πmkBT) / s₀` — 평탄 기판에서 한 cycle의 표면 반응을 완료하는 데 필요한 최소 노출량
* `E_required` — Gordon 모델이 계산한 HAR 구조 conformal 코팅 노출량

 **문헌 근거 (Cremers et al. 2019)** :

Gordon et al.은 HfO₂ ALD(TEMAHf, 200°C)에서 평탄 기판 포화에 3~43 L이 필요하다고 추정했습니다.
반면, 직경 0.17 μm, 깊이 7.3 μm의 홀(AR≈43:1)에 conformal 코팅을 하려면 9,000 L이 필요했습니다.
이는 평탄 기판 대비 약 **200~3,000배**에 해당합니다.

**이 배수가 크다고 해서 공정에 문제가 있는 것이 아닙니다.** 이것은 고종횡비 구조에서의 확산 제한(diffusion limitation)에 의한 물리적 필연입니다.

| Dose Multiple | 의미                     | 대표 조건             |
| ------------- | ------------------------ | --------------------- |
| < 10×        | 평탄 기판과 큰 차이 없음 | 낮은 AR (< 5:1)       |
| 10~100×      | 도징 시간 연장 필요      | 중간 AR (5~20:1)      |
| 100~1000×    | 고종횡비 표준 범위       | DRAM cap, 3D NAND     |
| > 1000×      | 극한 조건                | 200L+ NAND, ultra-HAR |

 **핵심 구분 — Dose Multiple ≠ Saturation Ratio** :

이 두 개념이 혼동되는 것이 매우 흔합니다:

| 비율                       | 정의                  | 분모                  | ≥100% 의미                     |
| -------------------------- | --------------------- | --------------------- | ------------------------------- |
| **Dose Multiple**    | E_req(HAR) / E_flat   | 평탄 기판 포화 노출량 | HAR 코팅이 평탄 대비 N배 어렵다 |
| **Saturation Ratio** | E_actual / E_req(HAR) | Gordon Required       | **conformal 코팅 달성**   |

* **Dose Multiple** 200×: "이 구조는 평탄 기판보다 200배 많은 노출이 필요합니다" → Gordon Required에 이미 반영됨
* **Saturation Ratio** 100%: "현재 공급 노출량이 Gordon이 예측한 필요량의 100%입니다" → **이때 conformal 코팅 달성**

 **결론** : "평탄 대비 100~1000배"는 Dose Multiple에 해당하며, 이는 Gordon Required Exposure에 이미 포함되어 있습니다. Saturation Ratio가 100% 이상이면 추가 배수를 곱할 필요 없이 conformal 코팅이 달성됩니다.

### 6.4 포화도 (Saturation Ratio) 해석

실제 노출량 / 필요 노출량의 비율로, 코팅 완성도를 나타냅니다:

| 포화도  | 의미      | 코팅 프로파일         | 조치                                             |
| ------- | --------- | --------------------- | ------------------------------------------------ |
| ≥ 100% | 충분      | 바닥까지 균일한 GPC   | 현 조건 유지 (과잉이면 시간 낭비)                |
| 90~99%  | 거의 달성 | 하단 5~10% 두께 감소  | 도징 시간 약간 증가                              |
| 70~89%  | 부족      | 하단 상당 부분 미코팅 | 도징 시간 또는 압력을 1.2~1.5배 증가             |
| 50~69%  | 크게 부족 | 절반 이상 미코팅      | 근본적 조건 변경 필요                            |
| < 50%   | 심각      | 유의미한 코팅 불가    | Fill tank 용량, 펌프, precursor 전략 전면 재검토 |

### 6.5 Fill Tank Model A vs B 차이 해석

| Model             | 정의                                             | 의미               |
| ----------------- | ------------------------------------------------ | ------------------ |
| A: P_eq × t_dose | 노출량 상한 (챔버 압력이 P_eq로 일정하다고 가정) | 낙관적 추정        |
| B: ODE 수치적분   | 실제 챔버 압력 변화를 적분한 노출량              | 현실적 추정 (권장) |

 **두 값의 차이가 클 때 (>20%)** :

* t_dose < τ인 경우: 도징 시간이 짧아 fill tank 내 precursor가 완전히 소진되지 않음 → **도징 시간을 τ의 2~3배로 늘리면** 노출량 효과적 증가
* t_dose > τ인 경우: 펌핑에 의해 precursor가 빠르게 제거됨 → **ΔP를 키우거나 fill tank 부피를 늘려** E_max 자체를 높여야 함

 **두 값의 차이가 작을 때 (<20%)** :

* Fast-fill 근사가 잘 적용되는 조건. 현재 밸브/배관 설계가 적절

### 6.6 침투 깊이 시간 해석

 **l(t) 그래프** :

* 파란 실선: 시간에 따른 침투 깊이
* 빨간 점선: 구조 전체 깊이 L
* 빨간 점: 현재 t_dose에서의 침투 깊이
* 초록 점선: 완전 코팅 달성 시각 t_full

 **핵심 지표** :

* **침투 깊이 @ t_dose** : 현재 도징 시간에서 precursor가 도달한 깊이. L과 비교하여 코팅 완성도 판단
* **t_full** : 구조 바닥까지 코팅이 완료되는 최소 시간. t_dose > t_full이면 OK
* **l_max (t→∞)** : Fill tank의 유한한 precursor 양으로 달성 가능한 최대 침투 깊이. l_max < L이면 어떤 도징 시간으로도 완전 코팅 불가

---

## 7. 그래프 탭 활용법

### Tab 1: Exposure vs EAR

 **질문** : "구조 shrink 후 EAR이 증가하면 노출량이 얼마나 필요한가?"

* 빨간 점선(현재 EAR)에서 필요 노출량 확인
* 고 EAR 영역에서 노출량 ∝ EAR² (quadratic 증가) 확인
* 차세대 구조 설계 시 필요 노출량 예측에 활용

### Tab 2: Pulse Time vs Pressure

 **질문** : "어떤 압력에서 도징해야 목표 시간을 달성하는가?"

* P × t = const 관계를 시각화
* 주황 영역(Kn < 1)에서는 Gordon 모델 부적합 → 해당 압력 범위 회피
* 장비의 펌핑 한계와 조합하여 최적 운전점 결정

### Tab 3: Pulse Time vs Temperature

 **질문** : "공정 온도를 바꾸면 도징 시간이 어떻게 변하는가?"

* T↑ → 열속도↑ → 같은 노출량에 약간 더 긴 펄스 필요 (약한 의존성)
* 하단 그래프: λ vs T → Kn 개선 여부 확인
* ALD window와 조합하여 최적 온도 결정

### Tab 4: EAR Structure Comparison

 **질문** : "같은 AR의 Trench와 Hole 중 어느 쪽이 코팅하기 쉬운가?"

* Hole: EAR = AR (가장 어려움)
* Trench: EAR = AR/2 (양면 열림)
* Pillar: EAR = AR/(2√2) (가장 쉬움)
* 공정 설계 시 구조 형상 변경 효과를 정량 비교

---

## 8. Fill Tank 모델 상세

### 8.1 기본 원리

Fill tank 방식 ALD에서 precursor는 고압의 fill tank에서 저압의 반응 챔버로 확장됩니다:

```
dP_f/dt = −C(P_f − P_c) / V_f        (fill tank 압력 감소)
dP_c/dt =  C(P_f − P_c) / V_c − S·P_c / V_c   (챔버 압력 변화)
```

### 8.2 핵심 파라미터

* **P_eq = ΔP × V_fill / V_chamber** : 평형 챔버 압력. Fill tank과 챔버가 즉시 평형을 이루었을 때의 압력
* **τ = V_chamber / S_pump** : 시상수. 챔버에서 precursor가 펌핑되는 특성 시간
* **E_max = P_eq × τ** : 최대 노출량. Fill tank 한 번의 방출로 달성 가능한 이론적 최대 노출량

### 8.3 실무 판단 기준

| 조건               | 의미                  | 대응                           |
| ------------------ | --------------------- | ------------------------------ |
| E_max > E_required | 완전 코팅 가능        | t_dose를 t_full 이상으로 설정  |
| E_max < E_required | 완전 코팅 불가        | ΔP↑, V_fill↑, 또는 S_pump↓ |
| t_dose/τ < 0.5    | 도징 시간이 너무 짧음 | t_dose 증가가 효과적           |
| t_dose/τ > 3      | 추가 시간 효과 미미   | ΔP 또는 V_fill 증대 필요      |

---

## 9. 단위 테스트 실행

```bash
# Standalone runner (streamlit/plotly 설치 불필요)
python run_tests.py

# 또는 pytest 사용
pytest test_gordon_calculator.py -v
```

테스트 항목 (총 44개):

* 단위 변환 정확도
* Gordon EAR 계산 검증
* K_max 추정 공식 검증
* 평균 자유 경로 계산 (경계 조건 포함)
* Fill Tank ODE vs 해석해 일치 (5개 t_dose에서 상대 오차 < 2%)
* Model A ≥ Model B 부등식
* 누적 노출량 단조 증가 및 E_max 수렴
* 침투 깊이 √E 점근 거동
* find_t_full 정방향/역방향 검증
* 입력 검증 (물리적 불가능 입력 방어)

---

## 10. 제한사항 및 주의사항

### 모델 한계

1. **s₀ = 1 가정** : Gordon 모델은 sticking coefficient가 1인 경우(완전 diffusion-limited)의 결과를 제공합니다. 실제 s₀ < 1이면 침투 프로파일이 더 완만해지므로, 모델이 제공하는 노출량은 **보수적 상한(하한 추정치)**입니다
2. **1차 Langmuir 흡착** : θ(coverage) 의존 반응 확률 s(θ) = s₀(1−θ)만 고려. 다층 흡착이나 복잡한 반응 메커니즘은 미반영
3. **등온 가정** : 구조 내부의 온도 구배를 고려하지 않음
4. **이상 기체 가정** : 저압(< 수 Torr)에서 유효. 고압에서는 실제 기체 보정 필요

### 실무 주의사항

1. **Kn < 1 영역** : 점성 유동에서는 Gordon 모델이 부정확합니다. 이 영역에서는 CFD(전산유체역학) 시뮬레이션 권장
2. **Fill Tank 압력 측정** : P_before와 P_after는 fill tank의 게이지 값이며, 배관 사이의 압력 손실은 별도 고려 필요
3. **Precursor 분해** : 일부 precursor(특히 금속 유기물)는 고온에서 열분해될 수 있으며, 이 경우 ALD가 아닌 CVD 성분이 포함됩니다. Gordon 모델은 순수 ALD만 고려합니다
4. **표면 화학 변화** : 구조 내부의 표면 화학이 깊이에 따라 다를 수 있습니다 (예: 상부는 SiO₂, 하부는 Si₃N₄). 이 경우 s₀가 위치에 따라 변하므로 모델 정확도가 저하됩니다

---

## 11. 참고문헌

1. R. G. Gordon, D. Hausmann, E. Kim, J. Shepard, "A Kinetic Model for Step Coverage by Atomic Layer Deposition in Narrow Holes or Trenches," *Chem. Vap. Deposition*  **9** , 73–78 (2003)
2. V. Cremers, R. L. Puurunen, J. Dendooven, "Conformality in Atomic Layer Deposition: Current Status Overview of Analysis and Modelling," *Appl. Phys. Rev.*  **6** , 021302 (2019)
3. S. M. Sze, M. K. Lee,  *Semiconductor Devices: Physics and Technology* , 3rd Ed., Wiley (2012)
4. B. A. Neaman,  *Semiconductor Physics and Devices: Basic Principles* , 4th Ed., McGraw-Hill (2012)

---

*ALD Gordon Model Calculator v5 — 2025*
