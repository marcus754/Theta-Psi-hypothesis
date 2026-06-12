# 06. Наблюдаемые и сопоставление с данными

## Канон времени

- Фундаментальный объект времени в гипотезе: тензор Theta_munu (см. docs/11_derivation_path.md).
- В формулах этого раздела используется редуцированная FRW-переменная theta ≡ theta_eff, где theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- Все записи вида n(theta) в этом файле следует читать как n(theta_eff, Theta, Psi).
- Этот раздел описывает вычислительный редуцированный контур, а не полную тензорную динамику Theta_munu.

## Каноническая постановка (обязательно для всех разделов)

В текущем редуцированном вычислительном слое сильное поле отклик строится
непосредственно из базы Theta_munu, Psi:
J_refr(Theta,Psi) = Phi_eff(theta_eff) · I_grad(Theta,Psi),
L_refr = M_*^4 F(J_refr),
n = n(theta_eff, Theta, Psi).

Для данных в репозитории важно различать:
- empirical наблюдаемые: прямые data-facing величины, которые идут в falsifiability;
- diagnostic наблюдаемые: внутренние диагностики и сканs для tria >= и разработки.

### Краткий обзор
- Фоновые наблюдаемые: H(z), расстояния (SN, BAO). Единицы: 8πG=1, c=1; H0 переводить через omegas_from_km_s_Mpc.
- Рост структур: D(a) из D''+(2+d ln H/d ln a)D'−(3/2)Ω_m mu D=0; в базовом варианте mu=1 (тензорный сектор немодифицирован).
- Каузальность: c_T=1; скаляры люминальны в high‑k; сверхсвет запрещён.
- Сопоставление: сначала срез по «здоровью» (K≻0, c_s²≤1, ω²>0, PPN/BBN), затем фит CC/BAO/SN.
- Практика: построение H(a) и D(a) из Omegas; отчёт о единицах и конверсии H0 обязательны.

Цели сопоставления:
- Фон: H(z), BAO, SN — проверка расширения и эффективной энергии.
- Локальные тесты: PPN‑ограничения, гравитационные волны (c_T=1), тесты Солнечной системы.
- Спутниковая синхронизация: GPS clock correction как reference check; знак и порядок коррекции совпадают со стандартной релятивистской оценкой.
  В data-facing отчёте используем опубликованный эталон порядка 38 mus/day для GPS.
- Лабораторная и транспортная синхронизация: Pound-Rebka, Hafele-Keating и twin-paradox. Все три оформлены как reference checks.
- Орбитальная слабое поле динамика: Lense-Thirring, binary pulsar timing, LLR.
- Атомный / химический предел: Theta–Psi не должен менять стандартную релятивистскую химию; в этом режиме требуется n(theta_eff) -> 1, gamma_PPN = 1, beta_PPN = 1.
- Космологические пертурбации: спектры, c_s², стабильность.

Статус чтения этого блока:
- reference: GPS, Hafele-Keating, Pound-Rebka, Shapiro, solar redshift, solar limb deflection, Mercury, Lense-Thirring, binary pulsar timing, LLR;
- reference: atomic limit, twin-paradox, c_T=1;
- слабое поле consistency only.

Статус сильное поле:
- reference: public event-level отчет for M87*, Sgr A* and GWTC-3 echo-scale;
- diagnostic: conservative compact-bound скан from the EHT/LVK envelope;
- compact-bound скан currently has no allowed points.

Практика:
- Использовать scripts/run_grid_скан.py для предварительного среза параметров (gamma, V0, Q) с автоматическими проверками.
- Дальше — выбор MCMC/сканера для подгонки H(z), BAO, SN, с фиксированным H0 (через omegas_from_km_s_Mpc).

### Демонстрационные пайплайны (в кодовых единицах)

- CC (H(z)): fitting/model.py::Hz_dimless и fitting/data_io.py::load_cc_dimless.
- SN (mu(z)): fitting/model.py::mu_dimless и fitting/data_io.py::load_sn_dimless.
- BAO (D_V(z)): fitting/model.py::DV_dimless и fitting/data_io.py::load_bao_dimless.
- Комбинированный логлайк и MCMC: python -m fitting.run_mcmc --steps 200 --output results/fit_demo.txt.
- Для реального H0 передайте --H0_km_s_Mpc 70.0 (внутри будет использована omegas_from_km_s_Mpc), либо вручную --H0_dimless.
- Для сохранения трейса укажите --save_csv results/chain.csv (формат: step,gamma,V0,Q,logpost,accept).
- Проверка перигелия Меркурия:
  - аналитический PPN-блок: python scripts/run_mercury_precession_check.py;
  - численное интегрирование орбиты: python scripts/run_mercury_numeric_check.py.
- Совместная карта области параметров:
  python scripts/run_joint_region_скан.py (CC+SN+BAO+compact_star + diagnostic compact скан).

### Как подготовить реальные данные (офлайн)

1) Скачайте публичные таблицы (пример):
   - CC (H(z)): компиляция cosmic chronometers (z, H, σ_H) в км/с/Мпк.
   - SN: Pantheon/Pantheon+ (z, mu, σ_mu).
   - BAO: таблица D_V(z) (или приведите к D_V) с (z, DV, σ_DV).
2) Преобразуйте форматы в простые CSV (без сетевого доступа):
   - python -m fitting.convert_public --cc CC_public.csv --H0_km_s_Mpc 70.0 --outdir fitting/data/real_dimless
   - python -m fitting.convert_public --sn Pantheon.csv --outdir fitting/data/real_dimless
   - python -m fitting.convert_public --bao BAO_DV.csv --outdir fitting/data/real_dimless
   - сильное поле события empirical-уровня (EHT ring catalog -> прямой-модель на уровне событий):
     python -m fitting.convert_public --прямой_events strong_field_events.csv --outdir fitting/data/real_dimless
     где strong_field_events.csv имеет колонки:
     name,observable,цель,sigma,mass_msun,distance_m,source. В empirical-контуре сейчас поддерживается только observable=ring_uas.
3) Запустите комбинированный фит на конвертированных данных (замените пути --data/--sn/--bao).

### Одношаговый сильное поле пайплайн

Если нужен прямой прогон от таблицы событий до отчётов:
python scripts/run_прямой_public_пайплайн.py --прямой_events_csv strong_field_events.csv --outdir results/public_прямой

Будут автоматически собраны:
- прямой_event_catalog.json
- прямой_event_скан.csv/.md
- прямой_event_отчет.json/.md
- empirical_checks.md
- пайплайн_summary.md

Примечание:
- echo-подобные величины и conservative compact-bound сканs не входят в слабое поле empirical-контур; для сильное поле они остаются diagnostic-артефактами.

Примечание. Если BAO DV даны в Мпк, умножьте на безразмерный масштаб (например, H0/c) через --bao_scale_dimless при конвертации.

## Линейный рост структур (квазистационарный субгоризонтальный предел)

В переменной N=ln a уравнение роста для фактора D(a) имеет вид
D'' + (2 + d ln H / d ln a) D' − (3/2) Ω_m(a) mu(a,k) D = 0,

где штрихи — производные по N, а mu(a,k) — модификация гравитации (в GR mu≡1).
Для фоновой космологии с E(a)^2 = Ω_r a^-4 + Ω_m a^-3 + Ω_Λ получаем Ω_m(a) = Ω_m a^-3 / E^2(a), 
d ln H / d ln a = (1/2) d ln E^2 / d ln a = (a / (2 E^2)) dE^2/da,
dE^2/da = −4 Ω_r a^-5 − 3 Ω_m a^-4.

В текущей реализации базовый вариант — GR (тензорный сектор не модифицирован), поэтому mu(a,k)=1.
Реализация: см. src/growth.py, функция growth_factor(a_min, a_max, omegas, nsteps, mu_of_a=None)
возвращает (a_grid, D_norm) с нормировкой D(a=1)=1.

## Связь с P(k)

В первом приближении спектр сегодня может быть собран как P(k, z=0) ∝ D^2(a=1) · T^2(k) · P_prim(k),

где T(k) — transfer‑функция (может быть заимствована из ΛCDM на ранних стадиях), а P_prim(k) ∝ k^n_s.
Отличия от ΛCDM на уровне роста отражаются в D(a) (и, при желании, в mu(a,k)).

Примечание: для строгого расчёта P(k) и CMB‑обсерваблов нужен полный Больцман‑код; текущая цель — обеспечить корректный фон и устойчивые моды для первых феноменологических сравнений.
