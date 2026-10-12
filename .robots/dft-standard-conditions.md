# Reações publicadas (DFT) nas quatro condições padrão: tentativa revertida e requisitos

Revisão feita pelo Claude em 12/10/2026, a pedido do autor. A implementação vai ser refeita em outra sessão; este arquivo é o ponto de partida dela.

## 1. Objetivo

Rodar a microcinética com Monte Carlo das reações publicadas, usando as energias DFT (`imports/calc-dft.xlsx`, aba `1M STATE FINAL`), nas quatro condições padrão da biblioteca (A1, B1, A2 e B2; `library-conditions.md`). A saída alimenta as abas `validation A1 dft` … `validation B2 dft` da planilha de agregados. Essas abas geram o `DFT_for-labelling.csv`, que é rotulado com as mesmas regras da biblioteca (G1/G2/Z3).

As rodadas antigas equivalentes estão em `exports_persist/STANDARD_STATE_WRONG_SUBTRACT_INSTEAD_OF_SUM/validation_dft-opt_A1` … `_B2`. São 16 reações, com `im2-bitet0`, e **σ = 2** (desvio de `A` nas amostras: 2.06). Elas usam a conversão de estado padrão no sentido errado e servem só para comparação.

## 2. O que foi revertido

Dois commits de outro agente (Copilot), de 12/10, 16:19:
- `c587f54` "Run DFT reactions across standard conditions" (`main.ipynb`);
- `1d08840` "Update standard kinetics configuration" (`src/config.py`).

Os dois foram revertidos por `3421c92` e `9a8378e`. O código voltou a ser idêntico ao de `2da32bf`. Também foi descartada uma mudança não commitada do `main.ipynb`, que transformava a célula 5 (`DFT_DATASET_SETTINGS`) em célula *raw*. A versão do agente pode ser recuperada com `git show c587f54` e `git show 1d08840`.

O que o agente fez:
- na célula 4, deixou ativo o conjunto DFT e acrescentou `DFT_REACTIONS` (16 reações) e `DFT_RUN_ALL_STANDARD_CONDITIONS = True`;
- na célula 10, passou a selecionar todas as reações de `DFT_REACTIONS`, e não mais uma por vez;
- na célula 14, criou um laço sobre as condições: chama `configure_standard_condition(condition)`, relê temperatura, tempo, concentrações e σ do `config` e exporta para `exports/<condição>/`;
- na célula 15, desligou a exportação com `if False and RUN_MODE == "mode2"`;
- no `config.py`, mudou `PERTURBATION_SIGMA_DEFAULT` de 3.0 para 2.0; o resto é formatação (aspas simples para duplas).

## 3. Problemas encontrados

1. **σ errado nas reações DFT.** O laço chama `configure_standard_condition`, que põe em `PERTURBATION_SIGMA_DEFAULT` o `sigma` da condição (3.0 nas quatro). Depois o laço lê `SIGMA` desse valor, e as reações DFT rodariam com σ = 3. Por isso a mudança do padrão para 2.0 no `config.py` não tinha efeito. Ela só mudava o valor padrão dos argumentos de `src/perturbation.py`, que é fixado na importação e que o notebook não usa, porque passa `sigma` explicitamente.
2. **`*_mode2_samples.csv` e `*_timed_out.txt` deixavam de ser gravados.** O pipeline (`run_reaction_pipeline_mode2_parallel2`) grava sozinho só as trajetórias e as energias perturbadas. As amostras (`A`, `x_ts1`, `x_ts2`, `timed_out`) vêm da célula 15, que ficou desligada. A coluna `timed_out` é a contagem confiável de *timeouts*: o `.txt` só é criado quando há *timeout*, e em 12/10 havia `.txt` velhos de uma rodada anterior nas pastas `lib_b1` e `lib_b2`.
3. **Indentação na célula 15.** Na versão do agente, só a linha do `samples.csv` ficou dentro do laço das reações. Os blocos do `timed_out.txt`, das energias e dos agregados ficaram no laço das condições e seriam gravados apenas para a última reação de cada condição.
4. **As rodadas nas condições experimentais deixavam de funcionar.** Mesmo com `DFT_RUN_ALL_STANDARD_CONDITIONS = False`, o laço chamava `configure_standard_condition(STANDARD_CONDITION)` e sobrescrevia temperatura, tempo, concentrações e σ definidos pela célula 5 (`DFT_DATASET_SETTINGS`).

## 4. Problema anterior, que continua no código (não foi o agente)

**A célula 5 força σ = 2 para qualquer conjunto de dados.** Quando `DFT_SELECTED_REACTION` está definido, ela aplica `options.update({'SIGMA': 2.0})` ("DFT Only"), mesmo com o conjunto GFN1-xTB e `FORCE_DFT_CONDITIONS = True`. A rodada `exports_persist/1M_STANDARD-STATE_thesis/validation_xtb` (12/10, 16:05), das reações publicadas com energias xTB, saiu assim com σ = 2 (desvio de `A` de 2.06, amostras idênticas às da rodada DFT). A regra é σ = 3 para xTB (`dft-validation.md` §4a). O autor vai refazer essa rodada.

## 5. Requisitos para a nova implementação

1. **Dois modos, os dois funcionando:**
   - **condições experimentais:** uma reação por vez ou todas, com as condições da célula 5 (`DFT_DATASET_SETTINGS`), para o conjunto DFT e para o xTB (`FORCE_DFT_CONDITIONS`);
   - **quatro condições padrão:** as reações publicadas em A1, B1, A2 e B2, com exportação para uma subpasta por condição.

   Num modo, os parâmetros do outro não podem sobrescrever os seus.
2. **σ definido pelo conjunto de dados, e não pela condição:** 2.0 com `calc-dft.xlsx` e 3.0 com as energias GFN1-xTB. Vale para os dois modos e corrige também o problema da §4. Não mexer no `PERTURBATION_SIGMA_DEFAULT` nem no `sigma` de `STANDARD_CONDITIONS`, que são da biblioteca.
3. **Exportação completa por reação** (`samples.csv` com `timed_out`, `perturbed_energies.csv`, `trajectories.pkl.gz` e `timed_out.txt` quando houver), com a indentação certa.
4. **Lista de reações:** as 16 do conjunto DFT, com `im2-bitet0`, como nas rodadas antigas `validation_dft-opt_A1…B2`. As chaves com sufixo (`im4-oc1-1eq`, `im3-j29-2eq`) só existem no modo experimental.
5. Manter o `input()` que pergunta antes de apagar `exports/`, para não misturar arquivos de rodadas diferentes.

## 6. Como conferir a rodada nova

Scripts de referência: os usados em 12/10 (descritos em `phd-thesis/.robots/machine-analysis/ch2-microkinetics-standard-state-A1.md` §6 e §10.1).
- **σ:** desvio de `A` nas amostras de cerca de 2.06 (σ = 2, DFT) ou 3.09 (σ = 3, xTB), com semente 42.
- **Condição de cada pasta:** razão ROH:Mo e tempo final na trajetória sem perturbação. A temperatura e a concentração saem da queda inicial do bispirrolido comparada com k1(T)·c₀ (razão ≈ 1).
- **Energias:** em relação a `validation_dft-opt_A1…B2` antigas, TS1 e TS2 ficam 3.788 kcal/mol abaixo, e MAP e BIS não mudam.
- **Amostras:** as mesmas das rodadas antigas (mesmo `A` e `x`, se σ e semente forem iguais).
