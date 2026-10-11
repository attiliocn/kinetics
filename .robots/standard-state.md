# Estado padrão em solução: de 1 atm para 1 mol/L

Material de estudo de 11/10/2026, montado em duas etapas:
1. fontes que o autor reuniu ao verificar o parser das reações publicadas do cap. 2 (PDFs e textos de `external-data/`, conferidos pelo Claude);
2. levantamento do NotebookLM com prompt neutro, que não dizia qual sinal era o certo (NB-28, `notebooklm/RESPOSTA-NB28.txt`).

Contexto da verificação: `inventory/CHECKS.md`, "Parser das reações publicadas em DFT e estado padrão de 1 mol/L". Situação de cada citação: **[PDF]** conferida no PDF ou texto; **[captura]** vista em captura enviada pelo autor; **[NB-28]** devolvida pelo NotebookLM, ainda não conferida.

**Em uma linha:** a energia livre de cada espécie sobe RT ln(24.46) = 1.89 kcal/mol a 298 K ao passar do gás ideal a 1 atm para 1 mol/L; numa etapa em que duas espécies viram uma (associação, barreira bimolecular), ΔG e ΔG‡ **caem** 1.89 kcal/mol.

---

## 1. O problema

- Programas de química quântica (ORCA, Gaussian) calculam *G* com a termoquímica de gás ideal, no estado padrão de **1 atm** (às vezes 1 bar) e 298.15 K. Só a parte translacional depende da concentração: q_trans = (2πmk_BT/h²)^(3/2) V.
- Em solução, o estado padrão do soluto é **1 mol/L**. As constantes bimoleculares da microcinética, em M⁻¹ s⁻¹, pressupõem energias nesse estado padrão.
- Gás ideal a 1 atm e 298.15 K: V° = RT/p° = 24.46 L/mol, isto é, 1/24.46 = 0.0409 mol/L (0.0403 mol/L a 1 bar). Passar para 1 mol/L é **comprimir** o gás 24.46 vezes.

## 2. A correção

Por espécie (Bryantsev 2008, eqs. 9 e 10 [PDF]):

$$\Delta G^{\circ\rightarrow *} = -T\Delta S^{\circ\rightarrow *} = RT\ln\frac{V^\circ}{V^*} = RT\ln(24.46) = 1.89~\mathrm{kcal/mol}\quad(298.15~\mathrm{K})$$

$$G^{*} = G^{\circ} + \Delta G^{\circ\rightarrow *}$$

(° = 1 atm; * = 1 mol/L.) É a `eq:standardState` dos Fundamentos da tese.

- **Por que soma:** comprimir reduz a entropia translacional, ΔS = −R ln(24.46) = −6.35 cal mol⁻¹ K⁻¹, e *G* sobe.
- Valores: 1.894 kcal/mol = 0.003 019 E_h a 298.15 K; 2.98 kcal/mol a 423 K (Garay-Ruiz 2020 [NB-28]).
- Num ΔG, a correção se cancela quando o número de espécies não muda e vale **Δn × 1.89 kcal/mol** quando muda (Δn = espécies nos produtos − espécies nos reagentes).

| Etapa | Δn | ΔG(1 mol/L) − ΔG(1 atm) |
|---|---|---|
| A + B → AB (associação) | −1 | −1.89 |
| A + B → TS (barreira bimolecular) | −1 | −1.89 |
| AB → A + B (dissociação) | +1 | +1.89 |
| A → TS (barreira unimolecular) | 0 | 0 |
| A + B → C + D | 0 | 0 |

## 3. Três formas de conferir o sinal

1. **Potencial químico.** μ = μ°(p°) + RT ln(p/p°); com p = cRT, μ = [μ°(p°) + RT ln(c°RT/p°)] + RT ln(c/c°). Com c° = 1 mol/L, c°RT/p° = 24.46, e μ°(1 mol/L) = μ°(1 atm) + RT ln(24.46).
2. **Constantes de equilíbrio.** K_c = K_p (RT/p°)^(−Δn). Para A + B → AB (Δn = −1), K_c = 24.46 K_p: a associação é mais favorável na referência de 1 mol/L, e ΔG°_c = ΔG°_p − RT ln(24.46). É a eq. 10.55 de Cramer 2004: "ΔG°′ = ΔG° − RT ln(24.5)" para A + B → C [NB-28].
3. **Eyring com estado padrão explícito.** k = (k_BT/h)(c°)^(1−n) exp(−ΔG‡°/RT), com n = molecularidade. A constante física não depende de c°; para n = 2, ΔG‡(1 mol/L) = ΔG‡(1 atm) − RT ln(24.46). Somar 1.89 à barreira em vez de subtrair dá k cerca de 600 vezes menor (exp(3.79/0.593) a 298 K).

## 4. O que dizem as fontes

**Somam 1.89 kcal/mol a cada espécie (a barreira bimolecular cai):**

| Fonte | O que diz | Situação |
|---|---|---|
| Bryantsev, Diallo e Goddard 2008, eqs. 9–10 | ΔG°→* = RT ln(24.46) = 1.89 kcal/mol; "The conversion from an ideal gas standard state of 1 atm to a standard state of 1 M yields ΔG*_g = ΔG°_g + ΔG°→*", aplicado "to each gas-phase reactant and product" | [PDF] fonte primária |
| Eisenstein, Ujaque e Lledós 2020, eq. 2 | G = E(BS2) + [G(BS1) − E(BS1)] + ΔG^(1atm→1M); "Gibbs energy change for compression of 1 mol of an ideal gas from 1 atm to the 1 M solution phase standard state. It amounts 1.89 kcal mol⁻¹ at 298 K for each species and only affects reactions where the number of moles change (Δn ≠ 0)"; cita Bryantsev (ref. 65) | [PDF] *Top. Organomet. Chem.* 67, 1–38, Seção 5 |
| Hatanaka, Yoshimura e Maeda 2020 (capítulo do AFIR no mesmo volume, *Top. Organomet. Chem.* 67, 57–80, p. 63) | "Eq. 4 is the expression for unimolecular reactions, and the factor (RT/p0)^(n−1) is multiplied in the nth-order reactions, where p0 corresponds to the standard atmosphere, i.e., 1 atm." Para n = 2, multiplicar k por RT/p° = 24.46 L/mol equivale a baixar a barreira RT ln(24.46) | [PDF]; confirma o sentido pela constante de velocidade |
| Bursch et al. 2022, Seção 2 | "the solvation process refers to a change of standard state (from the usually assumed gas phase at 1 atm to 1 mol L⁻¹ in solution). Hence, a concentration-induced free-energy shift of RT ln V_M = 1.89 kcal mol⁻¹ (V_M: molar volume of an ideal gas, T = 298 K) has to be added (in addition to the pV term in H298 discussed above) for each species when common solvent models (COSMO-RS, SMD, ALPB/GBSA) are used with default settings. This generally affects all non-unimolecular reactions modelled in solution."; cita Harvey 2019 (ref. 134) | [PDF] |
| Harvey et al. 2019 | "the output from many quantum-chemistry codes needs to be corrected, as they usually assume a gas-phase standard state when computing the translational entropy […] for solutions, a standard state corresponding to [X]₀ = 1 M yields a different value; the switch from one standard state to the other is easily performed by adding a correction term RT ln V_m⁰ = 1.89 kcal mol⁻¹ at 298 K"; recomenda "the usual ideal gas expressions, combined with continuum solvent treatment of solvation Gibbs energy, and careful attention to the use of the correct reference standard state" | [PDF] |
| Marenich, Cramer e Truhlar 2009 (SMD), eq. 1 | ΔG°_S = ΔG_ENP + G_CDS + ΔG°_conc: "Since here the same concentration (1 mol/L) is used in both the gaseous and solution phases, ΔG°conc is zero. If we instead used a gas-phase standard state of 1 atm, ΔG°conc would be 1.89 kcal/mol." É a energia de transferência (gás → solução), por isso o +1.89 aparece nela (ver o parágrafo sobre a energia de solvatação abaixo) | [PDF] |
| Jaraíz et al. 2019, *ACS Catal.* 9, 4804, ESI Seção S1 (COPASI) | "The G0 values calculated with Gaussian at standard conditions of 1.0 atm, are transformed within COPASI (to use concentrations in units of mol/L instead of pressure) according to: G0[mol/L] = G0[atm] − RT·ln(24.5) and therefore, for a forward reaction of order (a+b): G0a[mol/L] = G0a[atm] − (a+b−1)·RT·ln(24.5)." **A fórmula da barreira está no sentido certo** (bimolecular, a + b = 2: a barreira cai RT ln 24.5); a primeira fórmula, lida como *G* de cada espécie, daria o sentido oposto e não leva à segunda. Ver §5.5 | [PDF] (ESI) |
| Garay-Ruiz e Bo 2020, eq. 6 | ΔG⁰(1M) = ΔG⁰(1atm) + RT ln(cRT); "situations with a high degree of association (with few entities) will be less destabilized than very disassociated (separated reactants) states" | [NB-28] |
| Lledós 2021, eq. 1 | último termo "amounts 1.89 kcal·mol⁻¹ for each species. It is only relevant when the number of moles changes in the reaction (Δn ≠ 0)" | [NB-28] |
| Cramer 2004, eqs. 10.53–10.55 | ΔG°′ = ΔG° + RT ln(Q°′/Q°); para A + B → C, ΔG°′ = ΔG° − RT ln(24.5) | [NB-28] |
| Cramer, lista CCL (30/11/2011) | G de −3.010 00 vai a "about −3.006 99 a.u." a 1 mol/L: "I need to pay the entropy penalty to compress my concentration from 1/24.5 M to 1 M, which is about 1.9 kcal/mol" | [PDF] `external-data/cramer.txt`; não citável como artigo |
| Fórum do ORCA, dois tópicos (12/2025 e 05/2025), respostas de Miquel Garcia (mgarcia) | Sobre a energia de solvatação, não sobre *G* de uma espécie: "If we denote 1mol/L as "*", and 1 atm as "0" […] deltaG^*_{gas->liquid} + 1.89 = deltaG^0_{gas->liquid}"; o SMD é parametrizado com 1 mol/L nas duas fases, "However, if the reference state for the gas-phase is of 1 atm, then you have to add a concentration term […] equal to 1.89 kcal/mol". A pergunta da usuária (Tetiana) parte da mesma premissa: "the gas-phase energy of the solute increases by RT ln 24.5 when converting from a 1 atm to 1 mol L⁻¹ standard state" | [PDF] `external-data/orca-forum-atm.txt`; não citável como artigo |

**Sentido inverso (subtraem de cada espécie; a barreira bimolecular sobe):**

| Fonte | O que diz | Situação |
|---|---|---|
| Besora e Maseras 2018, Seção 4 | "this can be done by simply increasing the barrier by 1.89 kcal/mol if the temperature of the process is 298.15 K" | [PDF] |
| Jaraíz 2020, "DFT-Based Microkinetic Simulations: A Bridge Between Experiment and Theory in Synthetic Chemistry", *Top. Organomet. Chem.* 67, 81–106 (volume de Lledós e Ujaque), pp. 86–87 | "The G0 values, calculated with a QM software at standard conditions of 1.0 atm, are transformed to use concentrations in units of mol/L instead of pressure in atm according to G0[mol/L] = G0[atm] − RT·ln(24.5)". É a primeira fórmula de Jaraíz 2019, **sem a segunda** (a da barreira): a busca no capítulo inteiro não encontrou a forma com (a + b − 1) nem outra menção a molecularidade (ver §5.5) | [PDF] |

**Confirma pela relação K_c × K_p, sem falar de solução:** Filot 2026 (*Introduction to Microkinetic Modeling*, livro voltado a catálise heterogênea e fase gasosa) [PDF]. Eq. 2.159: K_c(T) = (p°/k_BT)^Δν exp(−ΔG°_r/k_BT), com ΔG°_r no estado padrão de pressão (eq. 2.151, "typically p = 1 atm"). Para uma associação (Δν = −1), K_c = (k_BT/p°) exp(−ΔG°_r/k_BT); em unidades molares, k_BT/p° vira RT/p° = 24.46 L/mol a 298 K, e ΔG°(1 mol/L) = ΔG°(1 atm) − RT ln(24.46). Na teoria do estado de transição (Seção 5, eqs. 5.103–5.107), a constante de equilíbrio K‡ de uma etapa bimolecular carrega o fator V/(k_BT), ou seja, depende do volume de referência, e a unimolecular não. O livro não menciona o valor de 1.89 kcal/mol nem a conversão para solução.

**Não desempatam:**
- Ribeiro, Marenich, Cramer e Truhlar 2011 define o estado padrão de concentração fixa (eqs. 8–9) e dá barreiras em 1 M, sem a fórmula da conversão [PDF].
- Ryu, Park, Kim, Park, Kim e Baik 2018 (tutorial de *Organometallics*) não trata da conversão de 1 atm para 1 mol/L. Diz que os cálculos supõem 1 M: "Quantum chemical calculations are performed assuming a standard state that is ideally normalized. That is, when two reactants are brought together to react, the concentration is implied to be identical and normal: i.e., 1 M" [PDF]. A frase não confere com a termoquímica dos programas, que usa o gás ideal a 1 atm (§1). O tutorial trata de outra correção, a das concentrações reais (§6).

**Cuidado com a energia de solvatação (fórum do ORCA):** ΔG_solv é uma diferença *líquido − gás*, então o sinal da correção nela parece invertido: ΔG°(gás a 1 atm → solução a 1 mol/L) = ΔG*(1 mol/L → 1 mol/L) **+ 1.89**, porque G°(gás, 1 atm) = G*(gás, 1 mol/L) − 1.89 entra com sinal negativo. Isso é coerente com somar 1.89 a cada espécie; quem lê "somar 1.89" sem ver a que grandeza se aplica pode inverter a correção da energia de cada espécie. Na prática do ORCA: a termoquímica impressa no cálculo de frequências usa o gás ideal a 1 atm, mesmo com CPCM ou SMD, e por isso a conversão por espécie é necessária para chegar a 1 mol/L (Bursch 2022: "has to be added […] for each species when common solvent models […] are used with default settings").

**Leitura:** as fontes de termodinâmica de solvatação e os guias de boas práticas somam; as duas fontes com o sinal inverso são textos de microcinética, e o mesmo volume (Lledós 2020) traz os dois sinais em capítulos diferentes. No levantamento neutro (NB-28), só Besora 2018 apareceu com o sentido inverso.

**Questão encerrada (11/10):** fonte primária, guias de boas práticas, artigo do modelo SMD e fórum do ORCA conferidos no texto original; todos somam 1.89 por espécie.

**Para citar na tese:** Bryantsev 2008 (primária), Bursch 2022 e, se quiser, Harvey 2019 ou Marenich 2009; Bursch 2022 tem a frase mais direta ("has to be added […] for each species"). Não usar Besora 2018 para o sinal.

## 5. O que Besora e Maseras 2018 dizem e por que está errado

### 5.1 O texto

Besora e Maseras 2018, "Microkinetic modeling in homogeneous catalysis", *WIREs Comput. Mol. Sci.* 2018, e1372 (DOI 10.1002/wcms.1372), Seção 4, "Practical aspects", p. 10 [PDF]:

> "For first-order elementary steps, the rate constant has units of reciprocal time, often s⁻¹. The rate constant comes directly from the application of the Eyring–Polanyi equation with the corresponding ΔG‡, and no special care is required when taking directly the values resulting from application of a quantum chemistry code. Second-order elementary steps are trickier as the rate constant units have an additional dependence on reciprocal concentration, usually M⁻¹. In this case, the units for the rate constant derived from the Eyring–Polanyi equation will depend on the reference state that has been used in the calculation of ΔG‡. Unfortunately, the free energy values supplied by most quantum chemistry programs do not use 1 M as reference state, but the concentration of an ideal gas at 298.15 K and 1 atm. If we want to have concentrations in M in our rate constant, we have to change the reference state. Luckily, this can be done by simply **increasing the barrier by 1.89 kcal/mol** if the temperature of the process is 298.15 K. There is strictly no requirement to use the 1 M reference state, and a consistent system could be deviced with a different reference state and different units for the rate constants. […] However, we advice placing all free energies with respect to the 1 M reference state."

### 5.2 O que está certo no trecho

- Etapas de primeira ordem não dependem do estado padrão: k em s⁻¹ sai direto da equação de Eyring.
- Em etapas de segunda ordem, as unidades de k (M⁻¹ s⁻¹) dependem do estado padrão usado em ΔG‡.
- Os programas entregam *G* no estado padrão de gás ideal a 1 atm, e não em 1 mol/L.
- A conversão vale 1.89 kcal/mol a 298.15 K; qualquer estado padrão serve, desde que usado de forma consistente; é recomendável pôr tudo em 1 mol/L.

### 5.3 O que está errado: o sentido

"Increasing the barrier" inverte a correção. Para uma etapa bimolecular A + B → TS:

1. **Por espécie.** Cada espécie sobe RT ln(24.46) ao passar para 1 mol/L (§2). ΔG‡(1 mol/L) = [G_TS + k] − [G_A + k] − [G_B + k] = ΔG‡(1 atm) − k. A barreira **cai** 1.89.
2. **Pela equação de Eyring com estado padrão explícito.** k_r = (k_BT/h)(c°)⁻¹ exp(−ΔG‡/RT). A constante física é a mesma em qualquer referência. Com c° = 1/24.46 mol/L (gás a 1 atm), k_r = (k_BT/h) × 24.46 M⁻¹ × exp(−ΔG‡(1 atm)/RT) = (k_BT/h) × 1 M⁻¹ × exp(−[ΔG‡(1 atm) − RT ln 24.46]/RT). Para obter k_r em M⁻¹ s⁻¹ com c° = 1 mol/L, a barreira tem de **diminuir** RT ln(24.46).
3. **Leitura física.** A referência de 1 mol/L é 24.46 vezes mais concentrada que a de 1 atm; o encontro bimolecular é mais provável, e a barreira na nova referência é menor.

**Exemplo numérico** (ΔG‡(1 atm) = 20.0 kcal/mol, 298.15 K, RT = 0.5925 kcal/mol):

| Barreira usada em Eyring | Valor | k_r (M⁻¹ s⁻¹) | em relação ao correto |
|---|---|---|---|
| ΔG‡ − 1.89 (conversão correta) | 18.11 | ≈ 3.3 × 10⁻¹ | 1 |
| ΔG‡ sem conversão (1 atm) | 20.00 | ≈ 1.4 × 10⁻² | 1/24.5 |
| ΔG‡ + 1.89 (Besora) | 21.89 | ≈ 5.6 × 10⁻⁴ | 1/600 |

Com a regra de Besora, uma etapa bimolecular fica 24.5 vezes **mais lenta** que se nada fosse corrigido, em vez de 24.5 vezes mais rápida.

### 5.4 As fontes que contradizem o trecho

Bryantsev 2008 (eqs. 9–10), Eisenstein, Ujaque e Lledós 2020 (eq. 2), Bursch 2022, Harvey 2019 e Marenich 2009 [PDF ou captura]; Filot 2026 pela relação K_c × K_p [PDF]; a `eq:standardState` da tese (§4). O mesmo volume em que está o capítulo de Eisenstein (Lledós 2020) traz, em outro capítulo, "G°[mol/L] = G°[atm] − RT ln(24.5)", com o mesmo erro de Besora.

### 5.5 De onde o erro pode ter vindo (hipóteses, não verificadas)

- **A notação de Jaraíz 2019 (ESI, Seção S1) [PDF].** O texto dá duas fórmulas em sequência: "G0[mol/L] = G0[atm] − RT·ln(24.5)" e, "therefore", para uma etapa de ordem (a + b), "G0a[mol/L] = G0a[atm] − (a+b−1)·RT·ln(24.5)". A segunda está certa: a barreira bimolecular cai RT ln 24.5, e a unimolecular não muda. A primeira só é coerente com a segunda se G0 for lido como o ΔG de uma etapa de associação (Δn = −1), e não como o *G* de cada espécie; lida por espécie, ela daria o sentido oposto (TS − c − (A − c) − (B − c) = ΔG‡ + c). O capítulo de Jaraíz 2020 no volume de Lledós e Ujaque reproduz só a primeira fórmula (conferido no PDF), e, isolada, ela induz a subtrair a correção de cada espécie, que é o que o código da tese faz.

- **A energia de solvatação.** Na energia de transferência gás → solução, o termo aparece somado: ΔG°(gás a 1 atm → solução a 1 mol/L) = ΔG* + 1.89 (Marenich 2009; fórum do ORCA). "Somar 1.89" é certo nesse contexto e vira o oposto quando transportado para a barreira.
- **"Somar 1.89 a cada espécie" × "somar 1.89 à barreira".** As duas frases soam parecidas e têm efeitos opostos numa etapa bimolecular.
- Não foi encontrada errata do artigo.

### 5.6 Como o erro entra no código

Seguir a frase de Besora em hartree, antes de montar o perfil, significa **subtrair** 1.89 de cada espécie: ts − (a − k) − (b − k) − … dá ΔG‡ + k numa etapa bimolecular (álgebra do autor). É o que fazem `exchange-parser/dft-optimized.ipynb` e `exchange-parser/main.ipynb` (`G -= std_state_correction`) (§8).

### 5.7 O que não é afetado

- Etapas unimoleculares e ΔG de etapas sem mudança no número de espécies.
- A diferença entre duas barreiras bimoleculares (b2 − b1), e portanto a competição entre elas.
- Afetados: as barreiras bimoleculares em valor absoluto, as constantes de velocidade bimoleculares (fator de cerca de 600 em relação ao correto) e, com elas, as conversões e escalas de tempo da microcinética.

## 6. Correções relacionadas, que não são a de concentração

- **Solvente como reagente.** O solvente puro tem concentração própria: água, 55.34 mol/L. Bryantsev 2008 dá RT ln(55.34) = 2.38 kcal/mol por molécula de água a partir de 1 mol/L (eq. 14 e esquemas) [PDF]; Harvey 2019 dá RT ln(55.5/0.0409) = 4.26 kcal/mol a partir de "1 bar = p/RT = 0.0409 M" [PDF] (0.0409 mol/L é o valor a 1 atm; a 1 bar seriam 0.0403, diferença desprezível aqui). É o mesmo efeito, com referências diferentes. Stahn 2023 generaliza para qualquer solvente pela densidade e pela massa molar [NB-28]. Não se aplica às reações da tese (benzeno, sem participação do solvente).
- **Entropia da "gaiola".** Correções empíricas que reduzem a entropia translacional e rotacional em solução (½ ou ⅔ de S_gás, volume livre, Wertz, Mammen). Harvey 2019 e Lledós 2020 e 2021 recomendam não aplicá-las com solvatação contínua, por dupla contagem; o recomendado é IGRRHO + modelo contínuo + conversão de estado padrão, com mRRHO para modos de baixa frequência (Bursch 2022; Grimme 2021) [NB-28]. É o que dizem os Fundamentos da tese (com Lledós 2020).
- **Correção empírica de +3.2 kcal/mol.** Sciortino e Maseras 2023 e Pérez-Soto et al. 2020 somaram 3.2 kcal/mol à energia livre de todas as espécies para reproduzir rendimentos experimentais ("increasing the free energy of all the species by 3.2 kcal mol⁻¹") [NB-28]. É um ajuste empírico, não a conversão de estado padrão. O parser da biblioteca tem essa linha comentada, sem uso.
- **Alternativas à conversão das energias.** Multiplicar a constante de ordem n por (RT/p°)^(n−1) (volume de Lledós 2020) [NB-28]; no modelo de *energy span*, introduzir concentrações nas taxas ou nas energias, G = G° + RT ln[R] (Solel 2019) [NB-28].
- **Concentrações reais diferentes de 1 mol/L.** Ryu 2018 corrige o perfil pelas concentrações experimentais com Δ_rG = Δ_rG° + RT ln Q (eq. 9) [PDF]. Exemplo do tutorial: epóxido 1500 vezes mais concentrado que o catalisador na formação de carbonato cíclico; o efeito equivale a 5.1 kcal/mol, e a maior barreira passa de 33.2 para 28.1 kcal/mol (5.1 kcal/mol corresponde a RT ln 1500 a cerca de 80 °C; a 298 K seriam 4.3). Numa microcinética, as concentrações reais entram nas equações de velocidade, e essa correção não se aplica às energias; ela serve para ler um perfil de energia livre sem simulação cinética.
- **Entropia translacional na associação.** Ryu 2018: a termoquímica de gás ideal penaliza A + B → C em cerca de 11 kcal/mol de −TS a 25 °C (cerca de 15 a 100 °C e 25 a 300 °C), quase só pela entropia translacional; associações de três e quatro fragmentos, cerca de 22 e 33 kcal/mol. Um erro comum é tratar um complexo de reagentes fracamente ligado como uma partícula só, o que elimina essa penalidade [PDF].
- **Limite da própria conversão.** Supõe solução ideal (atividade igual à concentração); falha em soluções concentradas ou iônicas (Fundamentos).

## 7. Aplicação ao mecanismo da troca (cap. 2)

Referência: bispirrolido + 2 ROH (três espécies). Cada ponto do perfil é somado às espécies livres que o acompanham:

| Ponto | Espécies | Δn em relação à referência | Correção na energia relativa |
|---|---|---|---|
| a1 + ROH | 2 | −1 | −1.89 |
| TS1 + ROH | 2 | −1 | −1.89 |
| MAP + ROH + pirrol | 3 | 0 | 0 |
| a2 + pirrol | 2 | −1 | −1.89 |
| TS2 + pirrol | 2 | −1 | −1.89 |
| BIS + 2 pirrol | 3 | 0 | 0 |

- b1 = TS1 − (BP + ROH) e b2 = TS2 − (MAP + ROH) são bimoleculares: **as duas caem 1.89 kcal/mol** em relação a 1 atm.
- b2 − b1 não muda: a competição entre a primeira e a segunda troca não depende do estado padrão. Muda a escala de tempo: as etapas bimoleculares ficam 24.46 vezes mais rápidas que com as energias de 1 atm.
- Os ΔG de cada troca (BP + ROH → MAP + pirrol) não mudam (Δn = 0).

**Exemplo (entrada 1, `im5-oc1`, energias da `tab:RealSystemBestConformers`, kcal/mol):**

| | 1 atm | 1 mol/L (conversão correta) | 1 mol/L (código de 11/10) |
|---|---|---|---|
| a1 | 12.3 | 10.4 | 14.2 |
| TS1 = b1 | 12.1 | 10.2 | 14.0 |
| MAP | −12.3 | −12.3 | −12.3 |
| b2 | 7.2 | 5.3 | 9.1 |

## 8. O código da tese (verificação de 11/10/2026)

- Os parsers `exchange-parser/dft-optimized.ipynb` (reações publicadas) e `exchange-parser/main.ipynb` (biblioteca) fazem `G -= std_state_correction` em cada espécie, seguindo Besora 2018. Subtrair k de cada espécie e montar o perfil dá ΔG‡(1 atm) **+ k** numa etapa bimolecular (álgebra do autor: ts1 − k − (a − k) − (b − k) = ts1 − a − b + k).
- A conversão das fontes da §4 (primeira tabela) corresponde a `G += std_state_correction`, que dá ΔG‡(1 atm) − k.
- Diferença: 2 × 1.89 = **3.79 kcal/mol** em b1 e em b2; b2 − b1 não muda. A aba `1M STATE FINAL` lida pela microcinética tem o mesmo tratamento.
- Decisão do autor (11/10): refazer os cálculos com a conversão corrigida antes da versão da banca (decisão do autor, 11/10, que substitui a anterior de corrigir só depois da defesa). Detalhes em `inventory/CHECKS.md` e `machine-analysis/ch2-dft-reactions-sheet.md`.

## 9. Perguntas prováveis na arguição

- **"Por que a correção soma, se as barreiras caem?"** Soma por espécie; a barreira bimolecular tem uma espécie a menos no TS que nos reagentes, então cai 1.89 (tabela da §2).
- **"Qual o efeito no produto previsto?"** Nenhum na competição entre as duas trocas (b2 − b1 não muda); muda a conversão no tempo simulado.
- **"Por que não usar 1 atm direto?"** As constantes bimoleculares em M⁻¹ s⁻¹ e as concentrações iniciais da microcinética (0.01–0.1 mol/L) pressupõem o estado padrão de 1 mol/L.
- **"E a entropia em solução?"** A correção de concentração é diferente da correção da "gaiola"; esta é controversa com solvatação contínua (dupla contagem) e não foi aplicada.
- **"Por que não usou a correção empírica de 3.2 kcal/mol de Maseras?"** É um ajuste para reproduzir rendimentos de um sistema específico, não uma correção de estado padrão; aplicada igualmente a todas as espécies, ela só desloca energias em etapas com mudança no número de espécies.

## Fontes

| Fonte | Chave do `.bib` | Situação |
|---|---|---|
| Bryantsev, Diallo e Goddard 2008, *J. Phys. Chem. B* 112, 9709 | `bryantsevCalculationSolvationFree2008` | [PDF] em `external-data/` |
| Eisenstein, Ujaque e Lledós 2020, "What Makes a Good (Computed) Energy Profile?", *Top. Organomet. Chem.* 67, 1–38 | `eisensteinWhatMakesGood2020` | [PDF] em `external-data/` |
| Hatanaka, Yoshimura e Maeda 2020, capítulo do AFIR, *Top. Organomet. Chem.* 67, 57–80 | não está no `.bib` (só o volume) | [PDF] (volume em `external-data/`) |
| Bursch et al. 2022, "Best-Practice DFT Protocols for Basic Molecular Computational Chemistry" | `burschBestPracticeDFTProtocols2022` | [PDF] em `external-data/` |
| Marenich, Cramer e Truhlar 2009, "Universal Solvation Model Based on Solute Electron Density…" (SMD) | `marenichUniversalSolvationModel2009` | [PDF] em `external-data/` |
| Harvey et al. 2019, "Scope and Challenge of Computational Methods…" | `harveyScopeChallengeComputational2019` | [PDF] em `external-data/` |
| Garay-Ruiz e Bo 2020, "Revisiting Catalytic Cycles…" | `garay-ruizRevisitingCatalyticCycles2020` | [NB-28] |
| Lledós 2021, "Computational Organometallic Catalysis…" | `lledosComputationalOrganometallicCatalysis2021` | [NB-28] |
| Cramer 2004, *Essentials of Computational Chemistry*, eqs. 10.53–10.55 | `cramerEssentialsComputationalChemistry2004` | [NB-28] |
| Cramer, mensagem na lista CCL, 30/11/2011 | — | [PDF] `external-data/cramer.txt` |
| Ribeiro, Marenich, Cramer e Truhlar 2011, *J. Phys. Chem. B* 115, 14556 | — (conferir se está no `.bib`) | [PDF] em `external-data/` |
| Besora e Maseras 2018, *WIREs Comput. Mol. Sci.* | `besoraMicrokineticModelingHomogeneous2018` | [PDF] em `external-data/` |
| Jaraíz 2020, capítulo de microcinética em *New Directions in the Modeling of Organometallic Reactions*, *Top. Organomet. Chem.* 67, 81–106 | não está no `.bib` (só o volume, `lledosNewDirectionsModeling2020`) | [PDF] (volume em `external-data/`) |
| Sciortino e Maseras 2023; Pérez-Soto et al. 2020 (correção empírica de 3.2 kcal/mol) | `sciortinoMicrokineticModellingComputational2023`; `perez-sotoChallengeReproducingCalculations2020` | [NB-28] |
| Ryu, Park, Kim, Park, Kim e Baik 2018, "Pitfalls in Computational Modeling of Chemical Reactions and How To Avoid Them", *Organometallics* | `ryuPitfallsComputationalModeling2018` | [PDF] em `external-data/` |
| Filot 2026, *Introduction to Microkinetic Modeling* (eqs. 2.151, 2.159; 5.103–5.107) | `filotIntroductionmicrokinetic2026` | [PDF] em `external-data/` |
| Stahn 2023 (CPCM-X); Solel 2019 (*energy span*); Grimme 2021 (mRRHO) | não procuradas no `.bib` | [NB-28] |

**Sobre a resposta do NB-28:** o Gemini dá "E. Besora" (é Maria Besora); uma das fontes do notebook é uma conversa com o Gemini ("Calculating Solvation Energy", 2024), a retirar; não citou o capítulo de microcinética do volume de Lledós 2020 (talvez não selecionado); equações vieram como `\\(…\\)`.
