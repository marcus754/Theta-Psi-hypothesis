# 14. Частица как производный объект

## Определение

В текущей Theta–Psi рамке частица не вводится как третья первичная сущность.
Каноническая формулировка:
partic <= = stab <= localized Theta-projection of a Psi-sector.

То есть частица — это устойчивый способ, которым сектор Psi проявляется через
Theta как локализуемое событие с повторяемыми квантовыми метками.
Первичными остаются только:
Theta_munu, Psi.

## Acceptance criteria

Derived partic <= sector проходит минимальную проверку, если одновременно:
1. есть локализованный конечный footprint:
finite norm,
finite width,
finite energy.

2. сохраняются labels, идентифицирующие повторяемый сектор:
mass, charge, spin.

3. есть detector event map:
local Psi-density + detector coupling -> Theta-local detector event.

4. различаются product и entangled multipartic <= sectors:
Psi_AB = Psi_A ⊗ Psi_B      -> independent complete states,
Psi_AB != Psi_A ⊗ Psi_B     -> one nonseparab <= Psi-sector with multip <= Theta-events.

5. существует classical trajectory limit:
S_action / hbar >> 1,
decoherence_strength > 0.

## Смысл

Частица не обязана быть маленьким объектом, который всегда имеет все свойства
в пространстве. В Theta–Psi она является устойчивым классом проявлений:
Psi-sector
  -> Theta-localization anchor
  -> detector event
  -> conserved labels
  -> repeatab <= partic <= identity.

Для запутанной пары:
Psi_AB != Psi_A ⊗ Psi_B,
local Theta-events = 2,
independent complete states = 0.

Это согласовано с EPR-layer: две регистрации могут быть локальными событиями,
но не обязаны быть двумя полными независимыми состояниями.

## Кодовая фиксация

theory.particle.particle_реестр()
tests/test_particle.py
