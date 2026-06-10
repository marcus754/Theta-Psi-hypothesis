# 06. Наблюдаемые и сопоставление с данными

## Канон времени

- Фундаментальный объект времени в гипотезе: тензор `Θ_{μν}` (см. `docs/11_derivation_path.md`).
- В формулах этого раздела используется редуцированная FRW-переменная
  `θ ≡ θ_eff`, где `θ_eff^2 = Θ_{μν} u^μ u^ν`, `u^μ u_μ = -1`.
- Все записи вида `n(θ)` в этом файле следует читать как `n(θ_eff)`.
- Этот раздел описывает вычислительный редуцированный контур, а не полную
  тензорную динамику `Θ_{μν}`.


## Каноническая постановка (обязательно для всех разделов)

В текущем редуцированном вычислительном слое strong-field отклик строится
непосредственно из базы `Θ_{μν}, Ψ`:

```
J_refr(Θ,Ψ) = Φ_eff(θ_eff) · I_grad(Θ,Ψ),
L_refr = M_*^4 F(J_refr),
n = n(θ_eff, Θ, Ψ).
```

Для данных в репозитории важно различать:
- `empirical` observables: прямые data-facing величины, которые идут в
  `falsifiability`;
- `diagnostic` observables: внутренние диагностики и scans для triage и разработки.

### Синапсис
- Фоновые наблюдаемые: H(z), расстояния (SN, BAO). Единицы: 8πG=1, c=1; H0 переводить через omegas_from_km_s_Mpc.
- Рост структур: D(a) из D''+(2+d ln H/d ln a)D'−(3/2)Ω_m μ D=0; в базовом варианте μ=1 (тензорный сектор немодифицирован).
- Каузальность: c_T=1; скаляры люминальны в high‑k; сверхсвет запрещён.
- Сопоставление: сначала срез по «здоровью» (K≻0, c_s²≤1, ω²>0, PPN/BBN), затем фит CC/BAO/SN.
- Практика: построение H(a) и D(a) из Omegas; отчёт о единицах и конверсии H0 обязательны.


Цели сопоставления:
- Фон: H(z), BAO, SN — проверка расширения и эффективной энергии.
- Локальные тесты: PPN‑ограничения, гравитационные волны (c_T=1), тесты Солнечной системы.
- Спутниковая синхронизация: GPS clock correction как reference check; знак и порядок коррекции совпадают со стандартной релятивистской оценкой.
  В data-facing отчёте используем опубликованный эталон порядка 38 μs/day для GPS.
- Лабораторная и транспортная синхронизация: Pound-Rebka, Hafele-Keating и twin-paradox.
  Все три оформлены как reference checks.
- Орбитальная weak-field динамика: Lense-Thirring, binary pulsar timing, LLR.
- Атомный / химический предел: `Θ–Ψ` не должен менять стандартную релятивистскую
  химию; в этом режиме требуется `n(θ_eff) -> 1`, `γ_PPN = 1`, `β_PPN = 1`.
- Космологические пертурбации: спектры, c_s², стабильность.

Статус чтения этого блока:
- reference: GPS, Hafele-Keating, Pound-Rebka, Shapiro, solar redshift,
  solar limb deflection, Mercury, Lense-Thirring, binary pulsar timing, LLR;
- reference: atomic limit, twin-paradox, c_T=1;
- weak-field consistency only.

Статус strong-field:
- reference: public event-level report for M87*, Sgr A* and GWTC-3 echo-scale;
- diagnostic: conservative compact-bound scan from the EHT/LVK envelope;
- compact-bound scan currently has no allowed points.

Практика:
- Использовать `scripts/run_grid_scan.py` для предварительного среза параметров (γ, V0, Q) с автоматическими проверками.
- Дальше — выбор MCMC/сканера для подгонки H(z), BAO, SN, с фиксированным H0 (через omegas_from_km_s_Mpc).

### Демонстрационные пайплайны (в кодовых единицах)

- CC (H(z)): `fitting/model.py::Hz_dimless` и `fitting/data_io.py::load_cc_dimless`.
- SN (μ(z)): `fitting/model.py::mu_dimless` и `fitting/data_io.py::load_sn_dimless`.
- BAO (D_V(z)): `fitting/model.py::DV_dimless` и `fitting/data_io.py::load_bao_dimless`.
- Комбинированный логлайк и MCMC: `python -m fitting.run_mcmc --steps 200 --output results/fit_demo.txt`.
- Для реального H0 передайте `--H0_km_s_Mpc 70.0` (внутри будет использована `omegas_from_km_s_Mpc`), либо вручную `--H0_dimless`.
- Для сохранения трейса укажите `--save_csv results/chain.csv` (формат: step,gamma,V0,Q,logpost,accept).
- Проверка перигелия Меркурия:
  - аналитический PPN-блок: `python scripts/run_mercury_precession_check.py`;
  - численное интегрирование орбиты: `python scripts/run_mercury_numeric_check.py`.
- Совместная карта области параметров:
  `python scripts/run_joint_region_scan.py` (CC+SN+BAO+compact_star + diagnostic compact scan).

### Как подготовить реальные данные (офлайн)

1) Скачайте публичные таблицы (пример):
   - CC (H(z)): компиляция cosmic chronometers (z, H, σ_H) в км/с/Мпк.
   - SN: Pantheon/Pantheon+ (z, μ, σ_μ).
   - BAO: таблица D_V(z) (или приведите к D_V) с (z, DV, σ_DV).
2) Преобразуйте форматы в простые CSV (без сетевого доступа):
   - `python -m fitting.convert_public --cc CC_public.csv --H0_km_s_Mpc 70.0 --outdir fitting/data/real_dimless`
   - `python -m fitting.convert_public --sn Pantheon.csv --outdir fitting/data/real_dimless`
   - `python -m fitting.convert_public --bao BAO_DV.csv --outdir fitting/data/real_dimless`
   - strong-field события empirical-уровня (EHT ring catalog -> forward-модель на уровне событий):
     `python -m fitting.convert_public --forward_events strong_field_events.csv --outdir fitting/data/real_dimless`
     где `strong_field_events.csv` имеет колонки:
     `name,observable,target,sigma,mass_msun,distance_m,source`.
     В empirical-контуре сейчас поддерживается только `observable=ring_uas`.
3) Запустите комбинированный фит на конвертированных данных (замените пути `--data/--sn/--bao`).

### Одношаговый strong-field pipeline

Если нужен прямой прогон от таблицы событий до отчётов:

`python scripts/run_forward_public_pipeline.py --forward_events_csv strong_field_events.csv --outdir results/public_forward`

Будут автоматически собраны:
- `forward_event_catalog.json`
- `forward_event_scan.csv/.md`
- `forward_event_report.json/.md`
- `empirical_checks.md`
- `pipeline_summary.md`

Примечание:
- `echo`-подобные величины и conservative compact-bound scans не входят в
  weak-field empirical-контур; для strong-field они остаются diagnostic-артефактами.

Примечание. Если BAO DV даны в Мпк, умножьте на безразмерный масштаб (например, H0/c) через `--bao_scale_dimless` при конвертации.

## Линейный рост структур (квазистационарный субгоризонтальный предел)

В переменной N=ln a уравнение роста для фактора D(a) имеет вид

```
D'' + (2 + d ln H / d ln a) D' − (3/2) Ω_m(a) μ(a,k) D = 0,
```

где штрихи — производные по N, а μ(a,k) — модификация гравитации (в GR μ≡1).
Для фоновой космологии с `E(a)^2 = Ω_r a^{-4} + Ω_m a^{-3} + Ω_Λ` получаем

```
Ω_m(a) = Ω_m a^{-3} / E^2(a),
d ln H / d ln a = (1/2) d ln E^2 / d ln a = (a / (2 E^2)) dE^2/da,
dE^2/da = −4 Ω_r a^{-5} − 3 Ω_m a^{-4}.
```

В текущей реализации базовый вариант — GR (тензорный сектор не модифицирован), поэтому μ(a,k)=1.
Реализация: см. `src/growth.py`, функция `growth_factor(a_min, a_max, omegas, nsteps, mu_of_a=None)`
возвращает `(a_grid, D_norm)` с нормировкой `D(a=1)=1`.

Связь с P(k)

В первом приближении спектр сегодня может быть собран как

```
P(k, z=0) ∝ D^2(a=1) · T^2(k) · P_prim(k),
```

где `T(k)` — transfer‑функция (может быть заимствована из ΛCDM на ранних стадиях), а `P_prim(k) ∝ k^{n_s}`.
Отличия от ΛCDM на уровне роста отражаются в `D(a)` (и, при желании, в μ(a,k)).

Примечание: для строгого расчёта P(k) и CMB‑обсерваблов нужен полный Больцман‑код; текущая цель — обеспечить корректный фон и устойчивые моды для первых феноменологических сравнений.
