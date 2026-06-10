# 15. Частица как derived object

## Определение

В текущей Θ–Ψ рамке частица не вводится как третья первичная сущность.

Каноническая формулировка:

```text
particle = stable localized Θ-projection of a Ψ-sector.
```

То есть частица — это устойчивый способ, которым сектор `Ψ` проявляется через
`Θ` как локализуемое событие с повторяемыми квантовыми метками.

Первичными остаются только:

```text
Θ_{μν}, Ψ.
```

## Acceptance criteria

Derived particle sector проходит минимальную проверку, если одновременно:

1. есть локализованный конечный footprint:

```text
finite norm,
finite width,
finite energy.
```

2. сохраняются labels, идентифицирующие повторяемый сектор:

```text
mass, charge, spin.
```

3. есть detector event map:

```text
local Ψ-density + detector coupling -> Θ-local detector event.
```

4. различаются product и entangled multiparticle sectors:

```text
Ψ_AB = Ψ_A ⊗ Ψ_B      -> independent complete states,
Ψ_AB != Ψ_A ⊗ Ψ_B     -> one nonseparable Ψ-sector with multiple Θ-events.
```

5. существует classical trajectory limit:

```text
S_action / hbar >> 1,
decoherence_strength > 0.
```

## Смысл

Частица не обязана быть маленьким объектом, который всегда имеет все свойства
в пространстве. В Θ–Ψ она является устойчивым классом проявлений:

```text
Ψ-sector
  -> Θ-localization anchor
  -> detector event
  -> conserved labels
  -> repeatable particle identity.
```

Для запутанной пары:

```text
Ψ_AB != Ψ_A ⊗ Ψ_B,
local Θ-events = 2,
independent complete states = 0.
```

Это согласовано с EPR-layer: две регистрации могут быть локальными событиями,
но не обязаны быть двумя полными независимыми состояниями.

## Кодовая фиксация

```text
theory.particle.particle_ledger()
tests/test_particle.py
```
