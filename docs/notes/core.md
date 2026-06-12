# Core (Θ–Ψ)

Единый рабочий документ: физическое ядро + операционный план + наблюдательные ориентиры.
смысл.md — авторский манифест и не редактируется.

## 1. Физическое ядро

- В качестве двух фундаментальных сущностей в гипотезе принимаются тензорное поле времени Θ_{μν} и энергетический субстрат Ψ.
- Θ_{μν} и Ψ обе первичны, но функционально несимметричны:
  Θ_{μν} задает структуру причинности и динамики, а в reduced слое
  динамизируется амплитуда R, являющаяся частью комплексного Ψ = R e^{iφ}.
- В полной версии канона Ψ рассматривается как комплексный order parameter
  Ψ = R e^{iφ}; в FRW-редукции остаётся амплитуда R и возможный фазовый
  заряд Q = a^3 R^2 \dot φ. В этом смысле |Ψ|^2 естественно читается как
  плотность амплитуды, а не как отдельная новая сущность.
- Скаляр θ в расчетах — это редуцированная проекция θ_eff
  (а не самостоятельная онтология).
- Операционное определение редукции:
  θ_eff^2 = Θ_{μν} u^μ u^ν, u^μ u_μ = -1.
- Каноническое определение Φ_eff:
  Φ_eff^can = ln(θ_eff/θ_0).
- Smooth Φ_eff, используемый в части solver'ов, считается только численной
  регуляризацией положительной ветви и не вводит новую физику.
- Канонический режим: FRW/weak-field воспроизводит проверенные локальные
  наблюдаемые тесты; отличие Θ–Ψ в преломляющем секторе включается через
  пространственные неоднородности, то есть через I_grad.
  Weak-field reference checks: GPS, Hafele-Keating, Pound-Rebka, Shapiro,
  solar redshift, solar deflection, Mercury, Lense-Thirring, binary pulsar
  timing, LLR, atomic limit, twin-paradox, c_T=1.
- Strong-field reference check: public event-level report on M87*, Sgr A* and
  GWTC-3 echo-scale. Conservative compact-bound scan remains diagnostic and
  currently has no allowed points.
- На однородном фоне I_grad=0, поэтому FRW сам по себе не включает
  преломляющий сектор.
- Время в гипотезе понимается не как шкала, а как положительный физический объект; допустимая
  динамика не должна переводить временной режим в нуль или отрицательную ветвь.
- Рабочая ветка редуцированного контура реализует этот принцип, в частности,
  через конечный индекс 1 <= n(θ_eff) < ∞.
- Weak-field matching зафиксирован как санитарное условие на наблюдаемые:
  n ≈ 1 + 2 Φ_eff в области слабых неоднородностей.
- Strong-field включается напрямую через составной инвариант режима
  J_refr = Φ_eff(θ_eff) I_grad(Θ,Ψ).
- Минимальный стационарный strong-field сектор задается функционалом
  L_stat = 1/2 I_grad - U(θ,Ψ) + αF(J_refr) и ведет к квазилинейной
  матричной системе для θ'' и Ψ''.
- Единый covariant completion, который покрывает и FRW, и strong-field:


S_cov = S_red + S_m[g_hat, matter]

- g и g_hat в этом completion не являются первичными сущностями
  гипотезы. Они используются как техническая ковариантная запись единого
  закона распространения процессов в Θ–Ψ-среде. Онтологический базис
  остается только Θ_{μν} и Ψ.
- Горизонт не вводится как физический объект Θ–Ψ. В сильном режиме
  каноническая переменная времени остается на положительной ветви, а
  наблюдаемый эффект читается как большое, но конечное замедление процессов:
  конечные n, redshift и delay.
- Derived microphysics layer, выводимый из базы Θ_{μν}, Ψ, задаёт локальный
  action S_micro с response-scalar R_micro, выведенным как минимальный
  локальный quadratic deviation scalar из locality / positivity / weak-field
  normalization, and with the preferred-slicing measure measure_dt_d3x
  fixed as derived volume normalization. Matching-правило:


n^2 = Z_s / Z_t

- Рабочий мост читается как S_micro -> S_red -> S_cov; S_red остаётся
  reduced/covariant completion, а S_micro - derived layer from the
  Θ_{μν}, Ψ basis.
- Algebraic reduction is fixed by R_micro -> 0 and the canonical sector
  weights, not by free coefficient functions or free potentials.
- Reduction-check: A_Θ, B_Θ, A_Ψ, B_Ψ, Z_t, Z_s -> 1, n -> 1,
  local quadratic blocks восстанавливаются, heavy modes интегрируются out.
- In code the derived micro layer is pinned to the canonical ledgers:
  S_red_from_micro = S_red, S_cov_from_micro = S_cov.
- Ordered derivation ladder is fixed as
  S_micro -> weak_field_projection -> integrate_heavy_modes -> assemble_S_red -> add_S_matter_completion -> S_cov.
- Micro perturbation check lives in theory/micro_perturbation.py
  and requires positive kinetic/response coefficients, positive sound-speed
  squares, and n² = Z_s / Z_t > 0.
- Single strong-field consequence of the derived micro layer is the compact branch
  summary with finite redshift, finite delay, and ok_no_horizon == True.
  The explicit derived strong-field ledger is micro_strongfield_primary_ledger()
  and micro_to_strongfield_primary_bridge().
- Stationary strong-field branch is closed by
  strongfield_branch_closure_ledger() and
  test_strongfield_branch_closure_pins_horizonless_criteria. Closure means:
  positive θ_eff, finite n, finite redshift/delay, outgoing characteristic,
  finite energy, regular center, elliptic radial operator, forward-observable
  map, and explicit falsification criteria.
- Dynamic collapse is closed as a canonical evolution contract by
  dynamic_collapse_closure_ledger() and
  test_dynamic_collapse_closure_contract_pins_evolution_criteria. Closure
  means finite-energy initial data, evolution of Θ–Ψ, tracked min θ_eff
  and max n, preserved outgoing characteristic, no blow-up before branch
  decision, hyperbolic evolution operator, and only two allowed outcomes:
  settle to stationary branch or disperse to weak-field. A numerical PDE
  solver is an implementation of this contract, not a new canon.
- Symbol audit classifies the micro layer as:
  derived = projection/matching/reduction outputs, including measure_dt_d3x
  and R_micro;
  assumed = R, chi; diagnostic = the compact-branch consequence.
- Замыкание пунктов 1-3 фиксируется так:


S_micro
  -> universal Z_t(R_micro), Z_s(R_micro)
  -> n^2 = Z_s/Z_t
  -> g_hat как technical propagation record
  -> weak-field observable coefficients


- Пункт 1 (g/g_hat): g и g_hat не добавляют третью сущность, если все
  matter-секторы получают их только как запись одного и того же локального
  закона распространения, выводимого из Z_t,Z_s.
- Пункт 2 (универсальность): в текущем derived micro layer универсальность
  обеспечивается тем, что для matter-сектора вводится один общий набор
  response coefficients Z_t(R_micro), Z_s(R_micro), а не отдельные
  коэффициенты для разных веществ или полей.
- Пункт 3 (weak-field): weak-field тест считается закрытым на уровне канона,
  если при R_micro -> 0 выполняется Z_t -> 1, Z_s -> 1, n -> 1, а
  первая ненулевая поправка к n нормируется как n ≈ 1 + 2 Φ_eff и тем
  самым даёт правильные наблюдаемые коэффициенты красного смещения, задержки,
  линзирования и орбитальной динамики.
- Статус закрытия: пункты 1-3 закрыты в каноне через
  propagation_law_closure_ledger() и тест
  test_propagation_law_closure_pins_points_1_2_3. Этот ledger запрещает
  независимую онтологию g_hat, фиксирует один общий Z_t,Z_s для
  matter-секторов и закрепляет weak-field нормировку до fit.
- Vacuum-offset sanity: абсолютная аддитивная константа в effective
  matter action не считается прямым Θ–Ψ источником. Source map зависит от
  response/boundary-разностей (R_micro, Z_t, Z_s, n, boundary data), а
  не от произвольного Λ_vac_offset. Казимир при этом сохраняется как
  boundary-разность спектров. Кодовая фиксация:
  vacuum_offset_source_ledger() и test_vacuum_offset_source_ledger_*.
- EPR/CHSH sanity: алгеброидный слой читает Θ как causal anchor/foliation,
  а Ψ_AB как неразложимое singlet-state. Первый тест фиксирует без
  fit-параметров: E(a,b) = -a·b, |S_CHSH| = 2√2, no-signalling marginals
  равны 1/2, а correlation generator лежит в ker ρ_Θ и не является
  сигналом. Слой также явно отвергает PR-box/super-Tsirelson корреляции,
  фиксирует GHZ/Mermin contextuality, foliation-order independence и no-fit
  angle profile. Кодовая фиксация: theory/epr_algebroid.py,
  tests/test_epr_algebroid.py, docs/13_epr_algebroid.md.
- Particle sanity: частица определяется как derived object, а не третья
  онтология: particle = stable localized Θ-projection of a Ψ-sector.
  Минимальные проверки: finite localized footprint, conserved labels
  (mass, charge, spin), detector event map, различение product и
  entangled multiparticle sectors, classical trajectory limit. Кодовая
  фиксация: theory/particle.py, tests/test_particle.py,
  docs/14_particle.md.
- Anti-fit guardrails: слабополевый профиль n(θ_eff) закрыт как
  канонический asinh с фиксированным наклоном 2; альтернативы
  linear/exp/tanh2 остаются только diagnostic-only. Saturating family
  F_v1(J) with removable J_* normalization, ADM eps*/alpha* mixings, sigma-расширение, bridge
  calibration coefficients и inverse-BAO required multipliers тоже считаются
  diagnostic-only, пока не выведены из действия и не зафиксированы до fit.
  В каноне сейчас: canonical_free_profile_count = 0,
  canonical_bridge_fit_parameters = 0. Кодовая фиксация:
  theory/anti_fit.py, tests/test_anti_fit.py.
- Action-derivation status: из текущего минимального действия уже следует
  диагональный ADM high-k блок
  K = G = diag(3/θ² + γ², 1, R²) и
  M = diag(0, V0 + 3Q²/(a⁶R⁴), 0). Поэтому канонические
  epsK*, epsG*, alpha_k2* равны нулю; ненулевые значения требуют
  новых операторов в действии. Форма F_v1(J) и strong-field bridge
  scales пока не выведены: действие задает аргумент J_refr и требования к
  F, но не выбирает саму функцию F(J), а bridge требует решения
  sourced static Euler-Lagrange equations. Кодовая фиксация:
  theory/action_derivation.py, tests/test_action_derivation.py.
- Refractive semigroup principle: если положить S(J)=1-F(J) и потребовать
  S(J1+J2)=S(J1)S(J2), то при непрерывности получается
  S(J)=exp(-J/J_*), а значит F(J)=1-exp(-J/J_*). J_* здесь не физика,
  а нормировка остатка отклика. Это уже не выбор формы, а law for the
  remaining response. Кодовая фиксация: theory/action_derivation.py,
  tests/test_action_derivation.py.
- Базовый однородный лагранжиан в reduced amplitude form:


L_eff = 1.5 * (θ̇/θ)^2 + 0.5 * Ṙ^2 - 0.5 * V0 * R^2 + 0.5 * γ^2 * θ̇^2,
где θ = θ_eff — FRW-редукция тензора Θ_{μν}, а R — reduced amplitude
representative комплексного Ψ = R e^{iφ}.


- Энергия и давление:


ρ = Σ_i (∂L/∂q̇_i) q̇_i - L
p = L


- Математический критерий «без истинного горизонта» (при конечной энергии):


L_eff ⊃ (3/2) (θ̇/θ)^2  =>  ρ ⊃ (3/2) (θ̇/θ)^2.


Если вдоль физической траектории ρ(t) <= ρ_max < ∞, то


|d ln θ / dt| = |θ̇/θ| <= sqrt(2ρ_max/3) = K.


Интегрируя:


ln θ(t) >= ln θ(t0) - K|t-t0|,
θ(t) >= θ(t0) exp(-K|t-t0|) > 0  для любого конечного t.


Следствие: состояние θ=0 не достигается за конечное время при конечной
энергии. Это и есть содержательная форма тезиса о положительности времени как
объекта. В рамках интерпретации n=n(θ) это дает динамический барьер для
формирования истинного горизонта: получается асимптотический режим сильного
замедления, а не физически достижимая поверхность «невозврата».

- SSOT по формулам редуцированного контура: background/frw_symbolic.py.
- Канонический физический слой (тензорное время): docs/11_derivation_path.md.

## 2. Численная схема

- background/frw_integrator.py — адаптивный RK45 в переменной x = artanh(γθ/√3).
- background/frw_background.py — FRW-контур с масштабным фактором и Ω-компонентами; невалидная ветка убрана конструктивно через гладкую сатурацию A.
- scripts/run_grid_scan.py — базовый скан + health-фильтры.

## 3. Состояние и приоритеты

- Текущий статус:
  - математика и pipeline рабочие;
  - явных конфликтов на наблюдательных проверках нет;
  - статистически подтвержденного отличающегося выигрыша на реальных
    event-level данных нет.

- Ближайшие шаги идут строго по порядку:
  1. Добить строгий ковариантный вывод: одно действие, одна вариация, один тензор энергии-импульса, одно правило связи материи с ĝ_{μν}.
  2. Проверить устойчивость: no-ghost, no-gradient, no-tachyon, отсутствие лишних мод и корректную эллиптичность strong-field ветки.
  3. Свести наблюдательные различия к 1-2 действительно отличающим предсказаниям, а всё остальное оставить как диагностику.
  4. Только после этого возвращаться к fit и сравнению с observational baselines; раньше это превращает гипотезу в подгонку.

- Практический смысл этого порядка:
  - если пункт 1 не закрыт, то fit подгоняет модель, а не теорию;
  - если пункт 2 не закрыт, то хорошая подгонка может скрывать патологию;
  - если пункт 3 не закрыт, то модель расползается в десятки слабых прокси вместо одного сильного теста;
  - если пункт 4 нарушен, то статистика начинает заменять физику.

## 4. Наблюдательные ориентиры для gamma

- Рабочая рамка по масштабу:
  - лаборатория: чувствительность грубо до γ ~ 10^-6;
  - астрофизика: γ ~ 10^-10 ... 10^-8;
  - космология/ранняя Вселенная: γ <= 10^-18 (ключевой лимит);
  - сверхкомпактные объекты: проверка конечного n(θ) и конечного redshift/delay.

- Для строгой калибровки нужны реальные таблицы и явная единичная нормировка в fit-пайплайне.
