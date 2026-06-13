# Θ-Ψ Hypothesis Repository / Репозиторий гипотезы Тета-Пси

## English Version

### Θ-Ψ Hypothesis: A Minimal Two-Field Model of Time and Energy

The Θ-Ψ (Theta-Psi) framework provides a rigorous alternative to the geometric paradigm of General Relativity. It replaces the abstract concept of spacetime curvature with the physical dynamics of a temporal tensor field, Θ, and its energetic substrate, Ψ.

#### Key Concepts:
* **Refractive Origin of Gravity:** Gravity is not geometry; it is the refractive response of matter to a non-homogeneous temporal medium. The temporal tensor Θ defines the vacuum causal structure. Matter interacts with the spatial gradients of this field through a local refractive index, n. The gravitational acceleration is governed by the relation: a = grad ln n. This refractive law replaces the metric connection, providing a physical mechanism where matter is displaced toward regions of lower temporal frequency.
* **Non-Singular Strong Fields:** The theory prevents the formation of mathematical singularities. Due to the saturation of the refractive response, it predicts horizonless compact objects with finite redshift and a regular center. These objects remain causally open, offering falsifiable signatures for gravitational wave echo analysis.
* **Interdisciplinary Bridges:** The framework establishes a causal anchor for non-local correlations (EPR/CHSH), satisfies the Tsirelson bound (2√2) without free parameters, and links substrate phase dynamics directly to cosmological expansion.

---

### Prediction Comparison

#### Table 1. Compliance (Weak Field)
| Subject                   | Theta-Psi               | General Relativity       |
| :------------------------ | :---------------------- | :----------------------- |
| Mercury Precession        | 42.98 arcsec/century    | 42.98 arcsec/century     |
| Light Deflection          | 1.75 arcsec             | 1.75 arcsec              |
| Gravitational Redshift    | 2.12 x 10^-6            | 2.12 x 10^-6             |
| Shapiro Delay             | 250 microseconds        | 250 microseconds         |
| Grav. Wave Speed (c_T)    | Exactly 1.000...        | Exactly 1.000...         |

#### Table 2. Differences (Strong Field and Cosmology)
| Subject                   | Theta-Psi                       | GR and Lambda-CDM          |
| :------------------------ | :------------------------------ | :------------------------- |
| Event Horizon             | No (Physical Surface)           | Yes (Boundary of no return) |
| Redshift z                | Finite (approx 10^10)           | Infinite                   |
| GW Echoes                 | ~0.3 ms per 1 M_sun (2.5-4.0 ms)| 0 (Signal absorbed)        |
| Center Density            | Finite (No singularity)         | Infinite (Singularity)     |
| Planck's Constant (hbar)  | Derived (Charge quantum q0)     | Fundamental Postulate      |

---

### Main Directories

- `background/`, `checks/`, `fitting/`, `src/`, `theory/`: Computational core and verification.
- `observational/`: External layer for data comparison, falsifiability, and validation.
- `docs/`: Core hypothesis documentation (`00...16`).
- `docs/notes/`: Working notes, roadmap, status, and plans.
- `tests/`: Test suite.

---

### Architectural Layers

- `theory/`: Symbolics, covariant and reduced action, refractive-sector, and stability analysis.
- `background/`: FRW module and reduced cosmological dynamics.
- `checks/`: Physical tests (Mercury, PPN, stability) and profile-level diagnostics.
- `fitting/`: Inference layer: observables, likelihoods, MCMC, and bridge via `fitting/core_api.py`.
- `observational/`: External policy layer: validation, falsifiability, and model comparison.
- `scripts/`: Entrypoints and thin wrappers.

---

### Usage Commands

Environment Setup:

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

Full Validation Suite:

```bash
python scripts/entry_validate.py
```

Outputs:
- `results/validation_suite.md`
- Summary reports in `results/*.md`

Check Classes:
- `empirical`: Direct data-facing checks (`prediction_suite`, `forward_event_scan`, `forward_event_report`, `joint_region_scan`, `falsifiability_check`).
- `diagnostic`: Proxy/scan/calibration steps for internal debugging; these are not considered hard falsification criteria.

Mercury Check (Precession and numerical integrator):

```bash
python scripts/run_mercury_precession_check.py
python scripts/run_mercury_numeric_check.py
```

Parametric Fit (MCMC):

```bash
python scripts/entry_fit.py
```

Output:
- `results/fit_summary.txt`

Fit on current `real_dimless` data (Unified working mode):

```bash
python scripts/entry_fit.py \
  --steps 1200 --jobs 2 \
  --data fitting/data/real_dimless/cc_Ez_zlt1.csv \
  --sn fitting/data/real_dimless/sn_dimless_subsample.csv \
  --bao fitting/data/real_dimless/bao_DV_h0_over_c.csv \
  --sn_fit_mu0 \
  --output results/fit_real_Ez.txt
```

Important Notes:
- `cc_Ez_zlt1.csv` is already formatted as `E(z)=H/H0`, so `--H0_km_s_Mpc` is not required.
- The original `cc_dimless_zlt1.csv` uses Planck normalization (`~1e-60`) and requires a separate consistent mode.
- `bao_DV_h0_over_c.csv` is normalized to `DV_dimless = H0*DV/c`.
- By default, `entry_fit.py` runs the **joint strong-field** loop (`--joint_strong_field`) including `compact_star + forward + forward_events`.
- In the current canon, only `Θ` and `Ψ` are fundamental; legacy `sigma*` arguments are for CLI compatibility only.

Fast Fit for large datasets:

```bash
python scripts/entry_fit.py --jobs 4 --steps 800 --cc_stride 2 --sn_stride 20 --bao_stride 1 --output results/fit_summary_fast.txt
```

Strong-field Quick Mode:

```bash
python scripts/entry_strong_field.py --mode quick
```

Outputs:
- `results/forward_event_scan.md`
- `results/forward_event_report.md`
- `results/compact_constraint_scan.md`

Strong-field Public Mode:

```bash
python scripts/entry_strong_field.py --mode public --forward_events_csv fitting/data/forward_events_public_proxy.csv --outdir results/public_forward
```

---

## Русская версия

### Θ-Ψ Гипотеза — минимальная двухпольная модель времени и энергии

Это исследование того, можно ли построить согласованную физическую рамку, начиная только с двух фундаментальных сущностей и их асимметрии.

#### Базовые идеи:
* **Рефракционное происхождение гравитации:** Гравитация — это не геометрия, а рефракционный отклик материи на неоднородность временного режима. Тензор времени Θ определяет структуру причинности вакуума. Материя взаимодействует с градиентами этого поля через локальный индекс преломления n. Ускорение задается законом: a = grad ln n. Этот закон заменяет метрическую связность, объясняя движение тел в области с более низким темпом времени.
* **Безгоризонтные сильные поля:** Теория исключает формирование математических сингулярностей. Благодаря насыщению отклика, она предсказывает существование компактных объектов с конечным красным смещением и регулярным центром. Такие объекты остаются причинно открытыми, что дает проверяемые сигналы в виде эхо-колебаний гравитационных волн.
* **Междисциплинарные мосты:** Фреймворк устанавливает связь между причинным расслоением и квантовой контекстуальностью (EPR/CHSH), воспроизводит границу Цирельсона (2√2) без подгоночных параметров и связывает динамику фазы субстрата с космологическим расширением.

---

### Сравнение предсказаний

#### Таблица 1. Соответствие (слабое поле)
| Предмет                   | Тета-Пси               | ОТО                      |
| :------------------------ | :--------------------- | :----------------------- |
| Прецессия Меркурия        | 42.98 угл. сек за век  | 42.98 угл. сек за век    |
| Отклонение света          | 1.75 угл. сек          | 1.75 угл. сек            |
| Красное смещение          | 2.12 миллионных        | 2.12 миллионных          |
| Задержка Шапиро           | 250 микросекунд        | 250 микросекунд          |
| Скорость Грав. Волн (c_T) | Строго 1.000...        | Строго 1.000...          |

#### Таблица 2. Отличия (сильное поле и космология)
| Предмет                   | Тета-Пси                        | ОТО и Лямбда-CDM          |
| :------------------------ | :------------------------------ | :------------------------ |
| Горизонт событий          | Нет (физическая поверхность)    | Есть (граница невозврата) |
| Красное смещение z        | Конечное (около 10^10)          | Бесконечное               |
| Эхо слияния (GW echoes)   | ~0.3 мс на 1 M_sun (2.5-4.0 мс) | 0 (сигнал поглощен)       |
| Плотность в центре        | Конечная (без сингулярности)    | Бесконечная (сингулярность)|
| Постоянная Планка (hbar)  | Выводится (квант заряда q0)     | Фундаментальный постулат  |

---

### Основные каталоги

- `background/`, `checks/`, `fitting/`, `src/`, `theory/`: вычислительное ядро и проверки.
- `observational/`: внешний слой сравнения с данными, falsifiability и validation.
- `docs/`: основная документация гипотезы (`00...16`).
- `docs/notes/`: рабочие заметки, roadmap, статус, планы.
- `tests/`: тесты.

---

### Архитектурные слои

- `theory/`: символика, ковариантный и reduced action, refractive-sector, устойчивость.
- `background/`: FRW-модуль и редуцированная космологическая динамика.
- `checks/`: физические тесты (Меркурий, PPN, стабильность) и profile-level диагностика.
- `fitting/`: inference-слой: observables, likelihoods, MCMC, bridge через `fitting/core_api.py`.
- `observational/`: внешний policy-слой: validation, falsifiability, compare-models.
- `scripts/`: только entrypoints и thin wrappers.

---

### Готовые команды запуска

Подготовка окружения:

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

Полный контур проверок:

```bash
python scripts/entry_validate.py
```

Что на выходе:
- `results/validation_suite.md`
- сводные отчеты по шагам в `results/*.md`

Классы шагов:
- `empirical`: прямые data-facing проверки (`prediction_suite`, `forward_event_scan`, `forward_event_report`, `joint_region_scan`, `falsifiability_check`).
- `diagnostic`: proxy/scan/calibration шаги для внутренней отладки; они не считаются жёсткими критериями фальсификации.

Выборочная проверка Меркурия (прецессия и численный интегратор):

```bash
python scripts/run_mercury_precession_check.py
python scripts/run_mercury_numeric_check.py
```

Параметрический фит (MCMC):

```bash
python scripts/entry_fit.py
```

Фит на текущих `real_dimless` данных (единый рабочий режим):

```bash
python scripts/entry_fit.py \
  --steps 1200 --jobs 2 \
  --data fitting/data/real_dimless/cc_Ez_zlt1.csv \
  --sn fitting/data/real_dimless/sn_dimless_subsample.csv \
  --bao fitting/data/real_dimless/bao_DV_h0_over_c.csv \
  --sn_fit_mu0 \
  --output results/fit_real_Ez.txt
```

Важно:
- `cc_Ez_zlt1.csv` уже приведен к виду `E(z)=H/H0`, поэтому `--H0_km_s_Mpc` не нужен.
- `forward_event_*` сейчас является empirical-контуром только для `ring_uas`.
- В текущем каноне первичны только `Θ` и `Ψ`; старые `sigma*` аргументы не относятся к физическому ядру.
