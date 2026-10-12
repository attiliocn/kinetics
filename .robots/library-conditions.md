# Biblioteca: as quatro condições padrão da microcinética

Contexto para rodar a microcinética da biblioteca (conjunto GFN1-xTB, `imports/calculated-exchange.xlsx`) nas quatro condições usadas na tese (11/10/2026). Ver também `standard-state.md` e `dft-validation.md` nesta pasta.

---

## 1. As condições

Definição da tese: capítulo de Métodos, `sec:methodsMonteCarlo`, em `phd-thesis/chapters/methods.tex`. São duas condições, branda e forçada, cada uma com 1 e 2 equivalentes de fenol (ROH):

| Rodada | Condição | T (°C) | Tempo | bispirrolido (mol/L) | ROH (mol/L) | Equivalentes |
|---|---|---|---|---|---|---|
| A1 | branda | 25 | 6 h | 0.01 | 0.01 | 1 |
| B1 | branda | 25 | 6 h | 0.01 | 0.02 | 2 |
| A2 | forçada | 70 | 12 h | 0.1 | 0.1 | 1 |
| B2 | forçada | 70 | 12 h | 0.1 | 0.2 | 2 |

Nas quatro: `me2pyr`, `map` e `bis` começam em zero, σ = 3 kcal/mol (energias do protocolo composto GFN1-xTB//TPSS0), 1000 amostras e semente 42.

Na classificação da tese, A1 e B1 são as condições brandas e A2 e B2 as forçadas. Reações com conversão mediana de até 10 % numa condição branda são classificadas só pelas forçadas (`sec:methodsKineticClasses`).

Os valores acima conferem com os `config.py` salvos nas rodadas anteriores (`exports_persist/STANDARD_STATE_WRONG_SUBTRACT_INSTEAD_OF_SUM/library_A1_soft-1eq`, `library_A2_harsh_1eq`, `library_B1_soft-2eq`, `library_B2_harsh_2eq`). Essas rodadas usaram energias com a conversão de estado padrão no sentido errado (`standard-state.md`) e servem só para comparação.

## 2. Onde mudar a condição

Na biblioteca, as condições vêm dos padrões de `src/config.py`. O bloco `DFT_DATASET_SETTINGS` da célula 5 do `main.ipynb` vale só para a validação DFT.

| Parâmetro em `src/config.py` | A1 | B1 | A2 | B2 |
|---|---|---|---|---|
| `INITIAL_CONCENTRATIONS['bispyr']` | 0.01 | 0.01 | 0.1 | 0.1 |
| `INITIAL_CONCENTRATIONS['roh']` | 0.01 | 0.02 | 0.1 | 0.2 |
| `TEMPERATURE_DEFAULT` | 25 | 25 | 70 | 70 |
| `SIMULATION_TIME_DEFAULT` | `3600 * 6` | `3600 * 6` | `3600 * 12` | `3600 * 12` |
| `PERTURBATION_SIGMA_DEFAULT` | 3.0 | 3.0 | 3.0 | 3.0 |

Para trocar a condição sem editar vários valores manualmente, chame a função `configure_standard_condition()` no `src/config.py`:

```python
from src import config

config.configure_standard_condition("A1")
```

Os identificadores aceitos são `"A1"`, `"B1"`, `"A2"` e `"B2"` (a função também aceita minúsculas). A chamada atualiza `INITIAL_CONCENTRATIONS`, `TEMPERATURE_DEFAULT`, `SIMULATION_TIME_DEFAULT` e `PERTURBATION_SIGMA_DEFAULT`; faça-a antes de copiar esses padrões para as variáveis da célula 4 do `main.ipynb`. Depois de alterar o módulo, reinicie o *kernel* para que ele seja carregado de novo.

**Estado em 11/10/2026:** o `config.py` está na B2. O `main.ipynb` lê esses padrões na célula 4 (`TEMPERATURE_C`, `SIMULATION_TIME`, `INITIAL_CONCENTRATIONS`, `SIGMA`).

## 3. Como rodar cada condição

1. Célula 4: deixe ativo o bloco da biblioteca (`DATASET_PATH = "imports/calculated-exchange.xlsx"`, `DATASET_SHEET = "raw energies"`) e mantenha `FORCE_DFT_CONDITIONS = False`. Com o arquivo do DFT ou com `FORCE_DFT_CONDITIONS = True`, a célula 5 sobrescreve as condições e o σ.
2. Edite o `src/config.py` com a coluna da condição (§2) e reinicie o *kernel*.
3. Confira a lista "Configured options" impressa pela célula 4: temperatura, tempo, concentrações e σ.
4. Rode com `RUN_MODE = "mode2"`. O número de amostras está fixo em `n_samples=1000` na célula 14 (não usa `N_SAMPLES`).
5. Exporte cada condição numa pasta de rodada nova em `exports_persist/1M_STANDARD-STATE_thesis/`, seguindo os nomes antigos (`library_A1_soft-1eq`, `library_A2_harsh_1eq`, `library_B1_soft-2eq`, `library_B2_harsh_2eq`), com uma cópia do `config.py` usado, como nas rodadas anteriores. Não sobrescreva as pastas de `exports_persist/`.

**Planilha da biblioteca:** `imports/calculated-exchange.xlsx` (aba `raw energies`) foi refeita em 11/10/2026, 12:30. A versão anterior, com a conversão para 1 mol/L no sentido errado, ficou em `imports/calculated-exchange WRONG 1 MOLAR.xlsx` e não deve ser usada. A conversão certa soma RT ln 24.46 = 0.003 019 E_h a cada espécie (`standard-state.md`). **Conferido em 11/10/2026:** as 1235 células numéricas da aba `raw energies` diferem da versão antiga por +0.006 037 E_h (= 2 × 0.003 019, porque a antiga subtraía o valor), sem exceção. A planilha nova está no sentido certo.
