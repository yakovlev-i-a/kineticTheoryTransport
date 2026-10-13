# kineticTheoryTransport

Библиотека OpenFOAM-v14, реализующая **mixture-averaged** транспортную модель на
основе стандартной кинетической теории газов (Chapman–Enskog /
Lennard–Jones + поправка Стокмайера) для реагирующих многокомпонентных
течений.

Порт форка `standardKinetic/of-11` на OpenFOAM-v14 в виде **отдельной
библиотеки** (без правки нативных исходников OF).

## Атрибуция и лицензия

**Физика (формулы) — классика кинетической теории газов**, не наше изобретение:

- Chapman–Enskog для $\mu_k$ и $D_{jk}$: Hirschfelder, Curtiss & Bird (1954);
- интегралы столкновений $\Omega^{(1,1)*},\Omega^{(2,2)*}$ и их аппроксимации:
  Monchick & Mason (1961), Neufeld et al. (1972);
- правило Вильке для вязкости смеси: Wilke (1950);
- комбинационное усреднение теплопроводности: Mathur, Tondon & Saxena (1967);
- сводное изложение всех формул: **Kee, Coltrin & Glarborg, "Chemically
  Reacting Flow" (2003)**.

**Конкретная реализация (структура классов, числовые коэффициенты интегралов
столкновений, поправка Стокмайера) является портом** GPL-форка
`standardKinetic/of-11` из работы:

> Fei Wang, Liang An, Tat Leung Chan. *An OpenFOAM solver incorporating
> detailed transport model for reacting flow simulations*, Computer Physics
> Communications **309** (2025) 109472.
> https://doi.org/10.1016/j.cpc.2024.109472

Форк распространяется под **GPLv3** (полный текст в `LICENSE`). Настоящая
библиотека — производная работа и также лицензируется под GPLv3. Если вы
используете код в научной работе, пожалуйста, цитируйте Wang et al. (2025)
и Kee et al. (2003).

## Что считает

Из шести молекулярных параметров CHEMKIN (`geometry`, `welldepth`, `diameter`,
`dipole`, `polar`, `relaxation`) вычисляются:

- чисто-видовые $\mu_k$, $\kappa_k$, $D_{kk}$;
- бинарные $D_{jk}$ (Chapman–Enskog, интегралы столкновений
  $\Omega^{(1,1)*}$, $\Omega^{(2,2)*}$ с полярной поправкой Стокмайера);
- вязкость смеси $\mu$ — правило Вильке;
- теплопроводность смеси $\kappa$ — комбинационное усреднение (Mathur et al.);
- mixture-averaged $D_{k,mix} = (1-Y_k)/\sum_{j\ne k} X_j/D_{jk}$.

## Состав

```
kineticTheoryTransport/
├── Make/{files,options}
├── specie/transport/standardKineticTransport/   # species transport
├── mixtures/kineticTheoryMixture/               # runtime-интерфейс Djk
├── mixtures/kineticTheoryMulticomponentMixture/ # mixture (mu/kappa/Djk)
├── thermophysicalTransport/kineticTheoryFickianFourier/  # transport model
├── thermos.C                    # инстанцирование psi/rho thermo
├── thermophysicalTransport.C    # регистрация kineticTheoryFickianFourier
└── validation/                  # тест-кейс + валидация против Cantera
```

## Сборка

```bash
source /opt/openfoam14/etc/bashrc
cd kineticTheoryTransport
./Allwmake
```

`Allwmake` собирает `libkineticTheoryTransport.so` и устанавливает
user-level override шаблона `chemistryModel` (см. ниже).

## Использование в кейсе

`constant/physicalProperties`:

```cpp
libs ("libkineticTheoryTransport.so");

thermoType
{
    type            hePsiThermo;
    mixture         kineticTheoryMulticomponentMixture;
    transport       standardKinetic;
    thermo          janaf;
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}
```

В `constant/thermo` у каждого вида вместо `transport { As ...; Ts ...; }`
записываются кинетические параметры:

```cpp
transport
{
    geometry        1;      // 0 atom / 1 linear / 2 nonlinear
    welldepth       97.53;  // epsilon/kB [K]
    diameter        3.621;  // sigma [Angstrom]
    dipole          0.0;    // [Debye]
    polar           1.76;   // polarizability [A^3]
    relaxation      4.0;    // Zrot(298) [-]
}
```

`constant/thermophysicalTransport`:

```cpp
laminar
{
    model           kineticTheoryFickianFourier;
}
```

## Важная деталь: chemistryModel через dynamicCode

В OF-v14 класс `chemistryModel` (стандартная stiff-ODE химия) шаблонизируется
на thermoPhysics и компилируется **на лету** (dynamicCode). Чтобы эта
динамическая компиляция нашла `standardKineticTransport.H`, `Allwmake`
устанавливает user-level override шаблона в

```
~/.OpenFOAM/14/codeTemplates/dynamicCode/chemistryModel
```

добавляющий `-I<lib>/lnInclude` в `EXE_INC`. При первом запуске кейса
OpenFOAM один раз соберёт `standard<standardKinetic<...>>` и закеширует его.

## Валидация

`validation/validate_vs_cantera.py` сверяет формулы кинетической теории с
Cantera (gri30, mixture-averaged). Результат:

| величина | T=300 K | T=1800 K |
|---|---|---|
| $\mu$ | 1.24% | 0.40% |
| $\kappa$ | 1.28% | 0.22% |
| $D_{k,mix}$ | 0.14–1.35% | 0.01–0.28% |

Тест-кейс `validation/` — 2D counterflow flame (GRI 3.0, 53 вида,
325 реакций), запускается `foamRun -libs ("libkineticTheoryTransport.so")`.
