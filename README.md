# Θ-Ψ Hypothesis Repository

## Основные каталоги

- `background/`, `checks/`, `fitting/`, `src/`, `theory/`: вычислительное ядро и проверки.
- `observational/`: внешний слой сравнения с данными, falsifiability и validation.
- `docs/`: основная документация гипотезы (`00...11`).
- `docs/notes/`: рабочие заметки, roadmap, статус, планы.
- `tests/`: тесты.

## Архитектурные слои

- `theory/`
  символика, ковариантный и reduced action, refractive-sector, устойчивость
- `background/`
  FRW-модуль и редуцированная космологическая динамика
- `checks/`
  физические тесты и profile-level диагностика
- `fitting/`
  inference-слой: observables, likelihoods, MCMC, bridge через `fitting/core_api.py`
- `observational/`
  внешний policy-слой: validation, falsifiability, compare-models
- `scripts/`
  только entrypoints и thin wrappers

## Готовые команды запуска

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
- `empirical`: прямые data-facing проверки, которые входят в `falsifiability`.
- `diagnostic`: proxy/scan/calibration шаги для внутренней отладки; они не считаются жёсткими критериями фальсификации.

Параметрический фит (MCMC):

```bash
python scripts/entry_fit.py
```

Что на выходе:
- `results/fit_summary.txt`

Фит на текущих `real_dimless` данных из репы (единый рабочий режим):

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
- `cc_Ez_zlt1.csv` уже приведен к виду `E(z)=H/H0`, поэтому `--H0_km_s_Mpc` не нужен;
- исходный `cc_dimless_zlt1.csv` в планковской нормировке (`~1e-60`) и требует отдельного согласованного режима.
- `bao_DV_h0_over_c.csv` приведен к `DV_dimless = H0*DV/c` (из `DV/r_d` через множитель `r_d*H0/c`, где взято `r_d=147.09` Мпк и `H0=70` км/с/Мпк).
- По умолчанию `entry_fit.py` запускает **joint strong-field** контур
  (`--joint_strong_field`) и включает `compact_star + forward + forward_events`.
- В текущем каноне первичны только `Θ` и `Ψ`; старые `sigma*` аргументы
  не относятся к физическому ядру и должны рассматриваться только как legacy-совместимость CLI.

Ускоренный фит на тяжелых наборах (параллель + прореживание):

```bash
python scripts/entry_fit.py --jobs 4 --steps 800 --cc_stride 2 --sn_stride 20 --bao_stride 1 --output results/fit_summary_fast.txt
```

Strong-field быстрый режим (встроенные данные):

```bash
python scripts/entry_strong_field.py --mode quick
```

Что на выходе:
- `results/forward_event_scan.md`
- `results/forward_event_report.md`
- `results/compact_constraint_scan.md`

Важно:
- `forward_event_*` сейчас является empirical-контуром только для `ring_uas`.
- `compact_constraint_scan.md` остаётся diagnostic-отчётом по inferred proxy-ограничениям и не входит в `falsifiability`.

Strong-field public-режим (готовый CSV из репы, только empirical ring-события):

```bash
python scripts/entry_strong_field.py --mode public --forward_events_csv fitting/data/forward_events_public_proxy.csv --outdir results/public_forward
```

Что на выходе:
- `results/public_forward/pipeline_summary.md`
- `results/public_forward/forward_event_scan.md`
- `results/public_forward/forward_event_report.md`

Одноразовая калибровка фиксированных bridge-коэффициентов (`background -> compact`) на event-каталоге:

```bash
python -m scripts.run_bridge_global_calibration \
  --forward_event_catalog_json fitting/data/forward_event_catalog.json \
  --gamma 0.0020339622 --V0 0.012364227 \
  --rounds 40 --reg_strength 0.2 \
  --output_json fitting/data/bridge_coeffs_global.json \
  --output_md results/bridge_global_calibration.md
```

Мульти-точечная калибровка bridge по сетке фоновых режимов:

```bash
python -m scripts.run_bridge_global_calibration_multi \
  --forward_event_catalog_json fitting/data/forward_event_catalog.json \
  --gamma 0.0020339622 --V0 0.012364227 \
  --grid_gamma_factors 0.7,1.0,1.4 \
  --grid_v0_factors 0.7,1.0,1.4 \
  --weight_center 1.0 \
  --rounds 28 --reg_strength 0.02 \
  --output_json fitting/data/bridge_coeffs_global_multi_best.json \
  --output_md results/bridge_global_calibration_multi_best.md
```

Head-to-head сравнение с внешним baseline на одном likelihood (`CC + SN + BAO`) в режиме fixed-bridge:

```bash
python -m scripts.run_head2head_comparison \
  --cc fitting/data/real_dimless/cc_Ez_zlt1.csv \
  --sn fitting/data/real_dimless/sn_dimless.csv \
  --bao_mean fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt \
  --bao_cov fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt \
  --include_strong_field --strong_field_weight 1.0 \
  --bridge_coeffs_json fitting/data/bridge_coeffs_global.json \
  --theta_samples 500 \
  --output_md results/head2head_comparison.md \
  --output_json results/head2head_comparison.json
```

Head-to-head в режиме **одинакового состава** (`same_content`):

```bash
python -m scripts.run_head2head_comparison \
  --comparison_mode same_content \
  --omega_r 8.4e-5 --omega_b 0.049 --omega_c 0.251 --omega_l 0.7 \
  --cc fitting/data/real_dimless/cc_Ez_zlt1.csv \
  --sn fitting/data/real_dimless/sn_dimless.csv \
  --bao_mean fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt \
  --bao_cov fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt \
  --include_strong_field --strong_field_weight 1.0 \
  --bridge_coeffs_json fitting/data/bridge_coeffs_global.json \
  --theta_samples 180 \
  --output_md results/head2head_mode_same_content.md \
  --output_json results/head2head_mode_same_content.json
```

Фальсификационный `train/test` протокол (калибровка bridge только на train, оценка на test без ретюна):

```bash
python -m scripts.run_falsification_protocol \
  --comparison_mode same_content \
  --split_mod 3 --test_fold 1 \
  --theta_samples_train 180 \
  --bridge_rounds 32 --bridge_reg_strength 0.2 \
  --output_md results/falsification_protocol/report.md \
  --output_json results/falsification_protocol/report.json
```

K-fold фальсификация (прогон всех `test_fold` и сводка по `Delta BIC`):

```bash
python -m scripts.run_falsification_kfold \
  --comparison_mode same_content \
  --bridge_model constant \
  --split_mod 3 \
  --theta_samples_train 180 \
  --bridge_rounds 32 --bridge_reg_strength 0.2 \
  --output_md results/falsification_kfold/report_same_content.md \
  --output_json results/falsification_kfold/report_same_content.json
```

Примечание по `bridge_recenter_on_train`:
- по текущим k-fold прогонам более стабильный режим: `--no-bridge_recenter_on_train` (он теперь дефолт).

Что на выходе:
- `results/bridge_global_calibration.md`
- `results/head2head_comparison.md`
- `results/head2head_comparison.json`

Примечание по SN covariance:
- если есть полный SN covariance matrix, передай `--sn_cov <path>`; без него используется диагональная схема по `sigma_mu`.
