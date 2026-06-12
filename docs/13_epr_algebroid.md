# 13. EPR/CHSH и алгеброидный слой

## Цель

Этот раздел фиксирует минимальную проверяемую постановку образа:
Theta режет/слоит Psi как локальную причинную историю.


Проверка не должна превращаться в свободную модель корреляций. Поэтому первый
тест берётся максимально жёсткий: синглетная пара и CHSH.

## Минимальная структура

A -> M
E -> M
ρ_Theta : A -> TM
Psi_AB ∈ Γ(E)

Смысл:
- M — наблюдаемое пространство-время локальных событий;
- E — state-bundle, в котором записывается Psi;
- A — слой допустимых внутренних/state-преобразований;
- ρ_Theta — anchor, задаваемый Theta;
- ker ρ_Theta содержит correlation generators: они физически задают связность
  Psi, но не проецируются в управляемый сигнал в TM.

Каноническая фраза:
Psi несёт неразложимость состояния.
Theta задаёт causal anchor, через который эта неразложимость проявляется как
локальная история без сверхсветового сигнала.

## CHSH sanity-check

Для первого теста запрещены свободные параметры:
free_fit_parameters = 0
free_correlation_function = 0

Состояние фиксируется как синглет:
Psi_AB = singlet.


Единственное правило корреляции:
E(a,b) = - a · b.

Для оптимальных направлений:
a0 = (1,0,0),
a1 = (0,1,0),
b0 = (1,1,0)/sqrt(2),
b1 = (1,-1,0)/sqrt(2).

CHSH:
S = E(a0,b0) + E(a0,b1) + E(a1,b0) - E(a1,b1),
|S| = 2 sqrt(2).


Следовательно:
|S| > 2        нарушает локальную hidden-variab <= границу,
|S| <= 2√2     не превышает Tsirelson bound.

## No-signalling

Joint probabilities берутся стандартные для синглета:
P(A=α,B=beta|a,b) = 1/4 [1 + αbeta E(a,b)],  α,beta ∈ -1,+1.


Тогда маргиналы:
P(A=α|a,b) = 1/2,
P(B=beta|a,b) = 1/2.

Они не зависят от удалённой настройки. Значит, запутанность остаётся
физической Psi-связностью, но не становится каналом передачи информации.

## Anchor condition

Алгеброидная формулировка разделяет correlation и signal:
ρ_Theta(C_AB) = 0,
ρ_Theta(S_signal) ⊂ causal cone(TM).

C_AB — correlation generator. Он отвечает за неразложимость Psi_AB, но его
anchor равен нулю: он не даёт управляемого перемещения/сигнала в наблюдаемом
пространстве-времени.

S_signal — signal generator. Если сигнал есть, его проекция через ρ_Theta
должна оставаться внутри causal cone.

## Acceptance criteria

Модель проходит первый EPR sanity-check только если одновременно:
1. |S_CHSH| > 2;
2. |S_CHSH| <= 2√2;
3. no-signalling marginals равны 1/2;
4. ρ_Theta(correlation generator) = 0;
5. signal anchor causal;
6. нет свободной функции корреляций;
7. нет fit-параметров для CHSH.

## Что явно отвергается

No-signalling сам по себе недостаточен. PR-box имеет:
E00 = 1,
E01 = 1,
E10 = 1,
E11 = -1,
S_CHSH = 4.

Он не даёт управляемого сигнала, но превышает квантовую границу:
4 > 2√2.


Поэтому acceptance rule:
no-signalling = required but not sufficient,
|S_CHSH| > 2√2 -> reject.

Это защищает слой от превращения в модель, которая допускает любые
no-signalling корреляции.

## Foliation independence

Если Theta задаёт slicing/anchor, наблюдаемые вероятности не должны зависеть от
того, в каком описании событие A стоит раньше B или B раньше A:
P_AB(α,beta|a,b; A-before-B) = P_AB(α,beta|a,b; B-before-A).

Для синглета это выполняется, потому что joint probability задаётся одним
Psi_AB, а не динамическим сигналом от первого измерения ко второму.

## GHZ/Mermin conuality

CHSH показывает нарушение локальной hidden-variab <= границы статистически.
Следующий sanity-check — GHZ/Mermin, где противоречие с глобальными заранее
заданными значениями становится алгебраическим.

Для GHZ-ограничений:
XYY = +1,
YXY = +1,
YYX = +1,
XXX = -1.

Произведение квантовых правых частей равно -1. Но любая nonconual
hidden valuation с заранее заданными X_i,Y_i = ±1 даёт произведение +1,
потому что каждый локальный множитель входит дважды.

Итог:
global_hidden_valuation_exists = 0,
conuality_required = 1.


Это фиксирует, что Psi не является контейнером заранее готовых локальных
значений. Она задаёт контекстуальную state-связность.

## Мост к quantum layer

docs/04_quantum.md и этот раздел не вводят два разных Psi. Мост:
Psi_state_layer
  -> two_subsystem_entanglement
  -> singlet_projection
  -> fixed_CHSH_correlations
  -> Theta_anchor_no_signalling.

То есть EPR-layer является специальной двухподсистемной проверкой того же
state-layer, который в квантовом разделе описывается через состояния,
ковариации и моды.

## No-fit experimental profile

Для любого набора копланарных углов a0,a1,b0,b1 prediction задаётся сразу:
Eij = -cos(ai - bj),
S = E00 + E01 + E10 - E11.

Fit-параметров нет. Реальная таблица эксперимента может быть подставлена как
набор углов и измеренных корреляций, но теоретический профиль уже фиксирован.

Кодовая фиксация:
theory.epr_algebroid.epr_algebroid_реестр()
tests/test_epr_algebroid.py

## Что это даёт гипотезе

Этот слой не доказывает всю Theta–Psi гипотезу. Он фиксирует, что образ
Theta как causal slicing/anchor для Psi может воспроизвести известную EPR-картину
без произвольных ручек:
Bell violation: yes
Tsirelson bound: yes
no-signalling: yes
PR-box rejected: yes
GHZ global valuation rejected: yes
foliation-order dependence: no
free parameters: no

То есть алгеброид здесь не украшение, а способ отличить физическую связность
Psi от сигнальной причинности Theta.
