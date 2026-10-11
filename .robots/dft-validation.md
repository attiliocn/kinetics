# Validação DFT: reações publicadas

Contexto para rodar ou alterar a microcinética das reações publicadas calculadas em DFT (11/10/2026). Ver também `standard-state.md` nesta pasta.

---

## 1. O que é

Reações de troca de ligantes publicadas: bispirrolido Mo(NR)(CHCMe₂Ph)(Me₂Pyr)₂ + ROH → monoarilóxido (MAP) e bisarilóxido (BIS), com liberação de 2,5-dimetilpirrol. Quinze combinações de imido e ROH foram calculadas com o mecanismo *trans* imido em DFT completo: geometrias e frequências `ωB97X-D4/def2-SVP CPCM (benzeno)`, energia eletrônica `ωB97X-D4/def2-TZVP CPCM (benzeno)`, *G* a 298.15 K. A microcinética com Monte Carlo compara o produto previsto (MAP × BIS) com o relatado em cada artigo.

**Rótulos:**
- imidos: `im1` = 2,6-*i*Pr₂C₆H₃; `im2` = 2,6-Me₂C₆H₃; `im3` = C₆F₅; `im4` = 2-CF₃C₆H₄; `im5` = adamantila;
- ROH: `oc0`–`oc4` = tetra-hidrobinaftol (5,5′,6,6′,7,7′,8,8′-octa-hidro-1,1′-binaftol) monoprotegido com TBS, com H, F, Cl, Br e I nas posições 3 e 3′; `ad1` = 2-fenilfenol; `15` = 2,6-difenilfenol; `29` = 2,3,5,6-tetrafenilfenol; `31` = 4-bromo-2,3,5,6-tetrafenilfenol; `26` = 2,6-dimesitilfenol (HMTOH); `bitet0` = binaftol monoprotegido sem substituinte em 3,3′.

## 2. Arquivos

| Arquivo | O que é |
|---|---|
| `imports/calc-dft.xlsx`, aba `1M STATE FINAL` | *G* absolutos (E_h) por reação e espécie, em 1 mol/L. **Refeita em 11/10/2026 com a conversão no sentido certo**: +RT ln 24.46 = +0.003 019 E_h por espécie em relação a 1 atm (conferido em `roh_oc1` e no bispirrolido `im5`). |
| `imports/calc-dft.xlsx`, aba `1M STATE FINAL ERRADO` | versão anterior, com a conversão subtraída (barreiras bimoleculares 3.79 kcal/mol altas demais). Não usar. |
| `imports/calc-dft.xlsx`, abas `Sheet1` (1 atm) e `1M state` | intermediárias. |
| `exports_persist/1M_STANDARD-STATE_thesis/validation_dft-opt/` | **rodada corrigida (11/10/2026)**: aba `1M STATE FINAL` corrigida, σ = 2, condições da §5a. Acordo no produto majoritário em 11 de 14 reações com dado (erros nas 6, 8 e 12, com ΔG‡₂ − ΔG‡₁ dentro de ~2 kcal/mol); a 14 dá MAP:BIS ≈ 93:7 (relatado 9:1). |
| `exports_persist/STANDARD_STATE_WRONG_SUBTRACT_INSTEAD_OF_SUM/` | todas as rodadas **anteriores à correção** (validação DFT e xTB, biblioteca nas quatro condições). Só para comparação. |

Origem da planilha: parser `~/projects/ufscar/ufscar.schrock-exchange/notebooks/2026_parsed-notebooks/exchange-parser/dft-optimized.ipynb` (célula 8: a conversão para 1 mol/L deve ser `+=`; ver `standard-state.md`).

## 3. Problema conhecido (11/10/2026): rótulos com `_`

**Resolvido no `main.ipynb` em 11/10/2026:** a célula 7 normaliza o índice logo depois de `load_dataset` (`dataset.index = dataset.index.astype(str).str.replace('_', '-')`).

A aba `1M STATE FINAL` nova usa rótulos com sublinhado (`im5_oc1`, `im1_ad1`…). O `main.ipynb` usa hífen em tudo: `dft_labels` (célula 5), as chaves de `DFT_DATASET_SETTINGS` e a seleção `selected_rates.loc[[base_label]]` (célula 10). `load_dataset` (`src/data.py`) lê a primeira coluna como índice, sem conversão. **Do jeito que está, `dataset.loc[dft_labels]` (célula 7) falha com `KeyError`.** Corrigir uma das pontas: trocar `_` por `-` na coluna A da aba, ou normalizar o índice ao carregar (`dataset.index = dataset.index.str.replace('_', '-')`).

## 4. Como rodar (células 4 e 5 do `main.ipynb`)

**Estado em 11/10/2026:** célula 4 já aponta para `calc-dft.xlsx` / `1M STATE FINAL`; célula 5 com o bloco de condições da §5a e **σ = 2.0** ativo (o bloco de 3.0 ficou comentado). Para voltar ao conjunto GFN1-xTB, inverter as duas coisas.

1. Célula 4: comentar o bloco do conjunto GFN1-xTB e descomentar o do DFT (`DATASET_PATH = "imports/calc-dft.xlsx"`, `DATASET_SHEET = "1M STATE FINAL"`).
2. Célula 5: **σ = 2 kcal/mol para o DFT.** Hoje a célula força `SIGMA = 3.0` ("XTB Only") e deixa o bloco de 2.0 comentado: inverter os comentários. Regra: 2 kcal/mol para otimização DFT, 3 kcal/mol para GFN1-xTB.
3. `DFT_SELECTED_REACTION`: uma reação por vez, com as chaves da tabela da §5. Chaves com sufixo (`-1eq`, `-2eq`) e com `j` (`im3-j29`) são reduzidas ao rótulo do conjunto pela célula 10 (`im3-j29-1eq` → `im3-29`).
4. Exportar cada reação para uma pasta de rodada nova (não sobrescrever as de `exports_persist/`).

## 4a. Incerteza σ (decisão de 11/10/2026: σ = 2 kcal/mol)

**Decisão:** manter σ = 2 kcal/mol para as reações otimizadas em DFT. Critério: σ igual ao MAD da calibração de cada protocolo (DFT ≈ 2; protocolo composto GFN1-xTB = 3), o mesmo usado para o conjunto xTB.

**Base (calibração do sistema modelo contra DLPNO-CCSD(T1)/CBS, geometria ωB97X-D4/def2-SVP CPCM e energia def2-TZVP; energia eletrônica em vácuo, independente do estado padrão):**

| | N | MD | MAD | RMSD | parte aleatória √(RMSD² − MD²) |
|---|---|---|---|---|---|
| intermediários | 9 | −2.33 | 2.33 | 2.98 | 1.86 |
| barreiras | 6 | +1.79 | 1.86 | 2.39 | 1.58 |

**Como σ chega às barreiras** (`src/perturbation.py`): um único A ~ N(0, σ) desloca MAP e BIS; cada TS recebe A·x, com x ~ U(0, 1); os reagentes ficam fixos. Desvio resultante: intermediários σ; cada barreira (b1 = TS1; b2 = TS2 − MAP) σ/√3 ≈ 0.58σ; b2 − b1 ≈ 1.08σ. Com σ = 2: 1.15 kcal/mol por barreira, abaixo da dispersão aleatória do *benchmark* nas barreiras (1.58); σ ≈ 2.7–3 reproduziria essa dispersão.

**Ressalvas:** o viés sistemático (ωB97X-D4 superestima as barreiras em ~1.8 kcal/mol em média) não é representado pelo sorteio em torno de zero; amostra pequena (N = 9 e 6). Teste opcional: rodar as 17 reações também com σ = 3 e conferir se o produto majoritário muda.

## 5. Reações e condições (conferidas em 11/10/2026)

"#" = numeração usada nesta validação (17 comparações; 8/9 e 14/15 são a mesma combinação em duas condições). Conferência feita contra os procedimentos dos artigos e seus materiais suplementares (massas, volumes e quantidades de matéria).

| # | Chave no notebook | Comentário atual no notebook | Fonte (composto) | Notebook: BP / ROH (mol/L), t, T | Relatado | Situação |
|---|---|---|---|---|---|---|
| 1 | `im5-oc1` | Entry 1 (46) | Ibrahem et al. 2009, *J. Am. Chem. Soc.* 131, 3844 (5d) | 0.02 / 0.02, 1 h, 22 °C | C₆H₆, 0.02 M, 1:1, 22 °C, 1 h | ok |
| 2 | `im5-oc2` | Entry 2 (43) | idem (5a) | idem | idem | ok |
| 3 | `im5-oc3` | Entry 3 (44) | idem (5b) | idem | idem | ok |
| 4 | `im5-oc4` | Entry 4 (45) | idem (5c) | idem | idem | ok |
| 5 | `im2-oc0` | Entry 5 (49) | Wang et al. 2013, *Angew. Chem. Int. Ed.* 52, 1939 (7c) | 0.01 / 0.01, 1 h, 22 °C | "prepared in analogy to 7b": 1:1, 0.01 M em C₆D₆, 1 h (temperatura não escrita; caixa de luvas) | ok |
| 6 | `im2-oc1` | Entry 6 (48) | idem (7b) | 0.01 / 0.01, 1 h, 22 °C | SI-1 2.68 mg (5.00 µmol) + álcool 10 2.2 mg (5.00 µmol) em 0.500 mL de C₆D₆ (0.01 M), 1 h (temperatura não escrita; caixa de luvas) | ok |
| 7 | `im2-oc3` | Entry 7 (72) | Yu et al. 2012, *J. Am. Chem. Soc.* 134, 2788 (8a) | 0.02 / 0.02, 1 h, 22 °C | **não relatado** (preparado *in situ* "according to published procedures") | suposição: condição padrão |
| 8 | `im4-oc1-1eq` | Entry 9A (50) | Wang et al. 2013 (8) | 0.01 / 0.01, 1 h, 22 °C | suplementar: "prepared in analogy to 7b" (1:1, 0.01 M, 1 h); texto do artigo (eq. 1): **1.1 equiv** de ROH; o RMN do suplementar (38/16/46) é o mesmo resultado do texto | ambíguo: 0.01 (suplementar) ou 0.011 (texto); diferença pequena |
| 9 | `im4-oc1-2eq` | Entry 9B (50) | idem (11) | 0.02 / 0.04, 1 h, 22 °C | 10.0 + 20.0 µmol em 0.500 mL (0.02 / 0.04 M), 22 °C, 1 h | ok |
| 10 | `im1-oc1` | Entry 10 (40) | Sattely et al. 2009, *J. Am. Chem. Soc.* 131, 943 (39) | 0.02 / 0.02, 1 h, 22 °C | benzeno, 1:1, 22 °C, 1 h; concentração não confirmada | provável |
| 11 | `im1-oc3` | Entry 11 (21) | Malcolmson et al. 2008, *Nature* 456, 933 (13b) | 0.02 / 0.02, 1 h, 22 °C | 16.9 µmol em 845 µL (0.02 M), 1:1, 22 °C, 1 h | ok |
| 12 | `im1-ad1` | Entry 12 (24) | Lee et al. 2009, *J. Am. Chem. Soc.* 131, 10652 (Tabela 5, entrada 2) | 0.02 / 0.02, 1 h, 22 °C | 10.0 µmol em 500 µL (0.02 M), 1:1, 22 °C, 1 h | ok |
| 13 | `im3-j15` | Entry 13 (71) | Qin et al. 2023, *J. Am. Chem. Soc.* 145, 3748 (Mo(bisAr)-3) | 0.1 / 0.2, 6 h, 22 °C | 0.100 + 0.200 mmol em ~1 mL de tolueno, 22 °C, 6 h | ok |
| 14 | `im3-j29-1eq` | Entry 14A (68) | Nguyen et al. 2017, *Nature* 552, 347 (Mo-5), proporção citada por Qin 2023 | 0.03 / 0.03, 2 h, 22 °C | suplementar de Nguyen: preparo para RMN, 0.0280 + 0.0280 mmol em 1 mL de C₆D₆ (**0.028 M**), 22 °C, 2 h; preparo para reações, 0.100 + 0.100 mmol em 1 mL de C₆H₆ (**0.1 M**), 22 °C, 2 h | 0.03 M = preparo para RMN; **usar 0.1 M** (preparo *in situ* para reações; §5a) |
| 15 | `im3-j29-2eq` | Entry 14B (68) | Qin et al. 2023 (Mo(bisAr)-1) | 0.1 / 0.2, 12 h, 70 °C | 0.300 + 0.600 mmol em 3.00 mL de benzeno (0.1 / 0.2 M), 70 °C, 12 h | ok |
| 16 | `im3-j31` | Entry 15 (70) | idem (Mo(bisAr)-2) | 0.1 / 0.2, 12 h, 70 °C | 0.200 + 0.400 mmol em 2 mL de tolueno, 70 °C, 12 h | ok |
| 17 | `im3-j26` | Entry 16 (64) | Zhang et al. 2014, *J. Am. Chem. Soc.* 136, 16493 (9a) | 0.1 / 0.1, 1 h, 22 °C | solução estoque 0.1 M em benzeno, 22 °C, 1 h | ok |
| — | (sem chave) | Entry 8 (38) — removed | Sattely et al. 2009 (binaftol sem 3,3′) | — | — | fora da validação: o cálculo `im2-bitet0` usa o imido 2,6-Me₂, e o artigo usa 2,6-*i*Pr₂ |

Todas começam com `me2pyr`, `map` e `bis` nulos (correto). Os cálculos usam CPCM (benzeno); as reações 13 e 16 foram feitas em tolueno.

**Reação 14 (Mo-5):** Nguyen 2017 dá dois preparos *in situ* (massas conferem com o bispirrolido C₆F₅, 597.5 g/mol, e o tetrafenilfenol, 398.5 g/mol): para RMN a 0.028 M e para reações a 0.1 M, ambos 1:1, 22 °C, 2 h. **Nguyen não quantifica o produto** (só o sinal do alquilideno *syn*, δ 12.38). A proporção MAP:BIS ≈ 9:1 vem de Qin 2023, que repetiu "a reported procedure" sem dizer qual dos dois. Os dois preparos são *in situ* (sem isolamento). **Decisão (11/10):** simular o preparo para reações (0.1 M, 22 °C, 2 h), que é o procedimento anterior ao de Qin (a seletividade de Qin vem da mistura a −30 °C, que o modelo isotérmico não representa) e tem a mesma concentração do procedimento modificado de Qin. Inferência: Qin não diz qual preparo repetiu.

**Correções sugeridas no notebook** (decisão do autor): na `im4-oc1-1eq`, manter 0.01 (suplementar, "in analogy to 7b") ou usar 0.011 (texto do artigo, 1.1 equiv); para a `im3-j29-1eq`, usar 0.1 / 0.1 M (§5a). Os números entre parênteses nos comentários do notebook são entradas de um levantamento bibliográfico à parte; os das chaves `im4-oc1-2eq`, `im3-j29-1eq` e `im3-j29-2eq` não correspondem a essas reações (deveriam ser 51, 67 e 69).

## 5a. Condições de referência para a rodada corrigida (11/10/2026)

**Aplicadas no `main.ipynb` (célula 5) em 11/10/2026.**

Critério: a condição simulada é a do **preparo** do complexo, onde a distribuição relatada foi medida (e não a do uso em catálise). Bloco pronto para substituir `DFT_DATASET_SETTINGS` na célula 5 (as chaves são as atuais; mudam a `im3-j29-1eq` e os comentários).

| # | Chave | BP (mol/L) | ROH (mol/L) | t | T (°C) | Mudança |
|---|---|---|---|---|---|---|
| 1 | `im5-oc1` | 0.02 | 0.02 | 1 h | 22 | — |
| 2 | `im5-oc2` | 0.02 | 0.02 | 1 h | 22 | — |
| 3 | `im5-oc3` | 0.02 | 0.02 | 1 h | 22 | — |
| 4 | `im5-oc4` | 0.02 | 0.02 | 1 h | 22 | — |
| 5 | `im2-oc0` | 0.01 | 0.01 | 1 h | 22 | — |
| 6 | `im2-oc1` | 0.01 | 0.01 | 1 h | 22 | — |
| 7 | `im2-oc3` | 0.02 | 0.02 | 1 h | 22 | — (suposição: sem procedimento no artigo) |
| 8 | `im4-oc1-1eq` | 0.01 | 0.01 | 1 h | 22 | — (1:1 pelo suplementar; o texto do artigo diz 1.1 equiv → 0.011) |
| 9 | `im4-oc1-2eq` | 0.02 | 0.04 | 1 h | 22 | — |
| 10 | `im1-oc1` | 0.02 | 0.02 | 1 h | 22 | — (concentração não confirmada) |
| 11 | `im1-oc3` | 0.02 | 0.02 | 1 h | 22 | — |
| 12 | `im1-ad1` | 0.02 | 0.02 | 1 h | 22 | — |
| 13 | `im3-j15` | 0.1 | 0.2 | 6 h | 22 | — |
| 14 | `im3-j29-1eq` | **0.1** | **0.1** | 2 h | 22 | **0.03 → 0.1**: preparo *in situ* de Nguyen para reações (o 9:1 de Qin é do procedimento anterior ao dele, sem mistura a frio) |
| 15 | `im3-j29-2eq` | 0.1 | 0.2 | 12 h | 70 | — |
| 16 | `im3-j31` | 0.1 | 0.2 | 12 h | 70 | — |
| 17 | `im3-j26` | 0.1 | 0.1 | 1 h | 22 | — |

```python
    def _settings(bispyr, roh, hours, temp_c):
        return {
            'INITIAL_CONCENTRATIONS': {'bispyr': bispyr, 'roh': roh, 'me2pyr': 0.0, 'map': 0.0, 'bis': 0.0},
            'SIMULATION_TIME': 3600 * hours,
            'TEMPERATURE_C': temp_c,
        }

    # Conditions of the in situ preparation of each complex (where the reported distribution was measured)
    DFT_DATASET_SETTINGS = {
        'im5-oc1': _settings(0.02, 0.02, 1, 22),      # 1  Ibrahem 2009, 5d (F)
        'im5-oc2': _settings(0.02, 0.02, 1, 22),      # 2  Ibrahem 2009, 5a (Cl)
        'im5-oc3': _settings(0.02, 0.02, 1, 22),      # 3  Ibrahem 2009, 5b (Br)
        'im5-oc4': _settings(0.02, 0.02, 1, 22),      # 4  Ibrahem 2009, 5c (I)
        'im2-oc0': _settings(0.01, 0.01, 1, 22),      # 5  Wang 2013, 7c
        'im2-oc1': _settings(0.01, 0.01, 1, 22),      # 6  Wang 2013, 7b
        'im2-oc3': _settings(0.02, 0.02, 1, 22),      # 7  Yu 2012, 8a (no procedure reported; default)
        'im4-oc1-1eq': _settings(0.01, 0.01, 1, 22),  # 8  Wang 2013, 8 (SI: 1:1; main text: 1.1 equiv -> roh 0.011)
        'im4-oc1-2eq': _settings(0.02, 0.04, 1, 22),  # 9  Wang 2013, 11
        'im1-oc1': _settings(0.02, 0.02, 1, 22),      # 10 Sattely 2009, 39 (concentration not confirmed)
        'im1-oc3': _settings(0.02, 0.02, 1, 22),      # 11 Malcolmson 2008, 13b
        'im1-ad1': _settings(0.02, 0.02, 1, 22),      # 12 Lee 2009
        'im3-j15': _settings(0.1, 0.2, 6, 22),        # 13 Qin 2023, Mo(bisAr)-3
        'im3-j29-1eq': _settings(0.1, 0.1, 2, 22),    # 14 Nguyen 2017, Mo-5, in situ for reactions (MAP:BIS ~9:1, Qin 2023)
        'im3-j29-2eq': _settings(0.1, 0.2, 12, 70),   # 15 Qin 2023, Mo(bisAr)-1
        'im3-j31': _settings(0.1, 0.2, 12, 70),       # 16 Qin 2023, Mo(bisAr)-2
        'im3-j26': _settings(0.1, 0.1, 1, 22),        # 17 Zhang 2014, 9a
    }
```

Antes de rodar, conferir também os itens da §4 (σ = 2.0, aba `1M STATE FINAL`, rótulos `_` × `-`).

## 6. Resultado relatado em cada reação

Proporções em % (BP / MAP / BIS) e como foram medidas.

| # | Relatado | Medida |
|---|---|---|
| 1 | 46 / 7 / 47 | RMN de ¹H da mistura *in situ* |
| 2 | 22 / 56 / 22 | idem |
| 3 | 30 / 62 / 8 | idem |
| 4 | 29 / 67 / 4 | idem |
| 5 | 24 / 50 / 26 | RMN de ¹H *in situ* |
| 6 | 10 / 81 / 9 | idem |
| 7 | não quantificado | — |
| 8 | 38 / 16 / 46 | RMN de ¹H *in situ* |
| 9 | 3 / < 2 / 97 | idem |
| 10 | não quantificado (catalisador usado em metátese) | — |
| 11 | MAP, 95 % de conversão, sem BIS | conversão *in situ* |
| 12 | 28 / 72 / 28 (soma 128 % no artigo) | RMN de ¹H *in situ* |
| 13 | BIS 81 | rendimento isolado |
| 14 | MAP : BIS ≈ 9 : 1 | proporção citada por Qin 2023 para o procedimento original |
| 15 | BIS > 98 (84 isolado) | RMN e rendimento isolado |
| 16 | BIS 99 | rendimento isolado |
| 17 | não quantificado (MAP formado e usado como catalisador) | — |

## 7. Barreiras esperadas depois da correção do estado padrão

kcal/mol, 1 mol/L; b1 = TS1 − (BP + ROH); b2 = TS2 − (MAP + ROH).

| # | Reação | b1 | b2 |
|---|---|---|---|
| 1 | im5-oc1 | 10.2 | 5.3 |
| 2 | im5-oc2 | 20.8 | 21.5 |
| 3 | im5-oc3 | 13.5 | 20.1 |
| 4 | im5-oc4 | 16.1 | 23.6 |
| 5 | im2-oc0 | 12.0 | 15.5 |
| 6 | im2-oc1 | 15.2 | 13.2 |
| 7 | im2-oc3 | 14.6 | 25.8 |
| 8, 9 | im4-oc1 | 14.3 | 14.6 |
| 10 | im1-oc1 | 11.6 | 16.1 |
| 11 | im1-oc3 | 14.6 | 23.5 |
| 12 | im1-ad1 | 12.4 | 11.5 |
| 13 | im3-15 | 16.9 | 16.6 |
| 14, 15 | im3-29 | 17.9 | 19.9 |
| 16 | im3-31 | 16.1 | 20.2 |
| 17 | im3-26 | 20.5 | 37.0 |

Conferência rápida depois de carregar a aba: `barriers_df` (célula 10) deve dar esses valores (± 0.1). Em relação às rodadas antigas de `exports_persist/`, b1 e b2 caem 3.79 kcal/mol, e b2 − b1 não muda.
