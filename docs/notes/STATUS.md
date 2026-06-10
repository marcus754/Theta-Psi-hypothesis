# STATUS (v2.4)

Файл отвечает только на вопрос: какие артефакты считать источником истины.

Легенда:
- `Canonical` — можно использовать в выводах.
- `Provisional` — рабочий слой, не канон.
- `Exploratory` — исследовательский/legacy слой.

| Статус | Путь | Назначение |
|---|---|---|
| Canonical | `docs/notes/смысл.md` | Авторский манифест гипотезы (не редактировать автоматически) |
| Canonical | `docs/notes/core.md` | Единый конспект: ядро + план + наблюдательные ориентиры |
| Canonical | `docs/11_derivation_path.md` | Канонический физический слой: тензорное время `Θ_{μν}` |
| Canonical | `background/frw_symbolic.py` | SSOT по формулам L/EOM/ρ/p |
| Canonical | `docs/notes/frw_reduced_eom.md` | Читаемая форма фоновых уравнений |
| Canonical | `background/frw_integrator.py` | Адаптивная интеграция в переменной x |
| Canonical | `background/frw_background.py` | FRW-интегратор с Ω-компонентами |
| Canonical | `scripts/run_grid_scan.py` | Базовый параметрический скан |
| Canonical | `theory/adm_symbolic.py` | ADM-скалярный сектор редуцированного контура `θ_eff` |
| Canonical | `theory/covariant_action.py` | Символический ledger полного ковариантного `Θ–Ψ` действия |
| Canonical | `theory/perturbation_proof.py` | Полный ковариантный completion и perturbation ledger |
| Canonical | `theory/micro_action.py` | Derived microphysics ledger и локальное matching-правило |
| Canonical | `theory/micro_strongfield.py` | Single strong-field consequence of the derived micro layer |
| Canonical | `theory/micro_perturbation.py` | Derived microphysics perturbation ledger |
| Provisional | `theory/twofield_stability.py` | Минимальные FRW/stationary stability diagnostics |
| Exploratory | `theory/tensor_symbolic.py` | Символическая проверка редуцированных формул |
| Exploratory | `theory/sigma_gate.py` | Legacy gate для старого sigma-слоя |
| Canonical | `theory/optical_metric.py` | конечный `n(θ_eff)`, strong-field ветка |
| Provisional | `docs/00_overview.md` ... `docs/10_falsifiability.md` | Редуцированный вычислительный слой `θ_eff` (канон `Θ_{μν}` задан в `docs/11_derivation_path.md`) |
| Canonical | `scripts/run_validation_suite.py` | End-to-end валидация |
| Canonical | `scripts/run_falsifiability_check.py` | PASS/FAIL-критерии версии |
| Provisional | `checks/` | Veto/health-проверки (часть феноменологична) |
| Provisional | `fitting/forward_transfer.py` | explicit/event-level forward |
| Provisional | `fitting/strong_field_link.py` | background->compact link-модель |
| Provisional | `fitting/data/strong_field_bounds_registry.json` | strong-field реестр (часть записей exploratory) |
| Provisional | `docs/13_publication_draft.md` | Внешний paper-level черновик для публикации |
| Exploratory | — | сейчас нет выделенных exploratory-модулей |

## Правила использования

1. Формулы брать из `background/frw_symbolic.py`, не из пересказов.
2. Если `Canonical` и `Provisional` расходятся, приоритет у `Canonical`.
3. Документы `docs/notes` держать в минимальном наборе: `смысл.md`, `core.md`, `frw_reduced_eom.md`, `STATUS.md`.
4. Онтология времени: `Θ_{μν}` и `Ψ` — первичные сущности; `θ` в скриптах и FRW-доках трактуется как `θ_eff`.
