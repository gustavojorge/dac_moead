# Operador Worst Turbine Mutation no MoWFLOP

Este documento descreve o operador de mutação **Worst Turbine Mutation**: o que Ju e Liu propuseram, como o operador foi adaptado ao nosso problema biobjetivo e como foi implementado neste repositório.

Implementação: [worst_turbine_mutation.cpp](../../src/moead/src/global_modules/genetic_operators/mutation/worst_turbine_mutation.cpp).

> **Base bibliográfica.** O operador é baseado na **"relocation of the worst turbines" do algoritmo AGA de Ju e Liu (2019) [1]**, seções 3.3 (passo 4) e 6. O que foi mantido do artigo e o que foi adaptado ao nosso problema está na seção 2, e a referência completa está na seção 8.

---

## 1. O que Ju e Liu (2019) propuseram [1]

### 1.1 Motivação

Ju e Liu, em *Wind farm layout optimization using self-informed genetic algorithm with information guided exploitation* (Applied Energy 248, 2019, [doi:10.1016/j.apenergy.2019.04.084](https://doi.org/10.1016/j.apenergy.2019.04.084)), argumentam que a mutação do algoritmo genético convencional (CGA) é cega. Ela move uma turbina qualquer para uma posição qualquer e, por isso, pode piorar o layout (Fig. 6 do artigo). A proposta é que cada indivíduo "examine a si mesmo" e corrija o seu gargalo:

> *"The individuals will check the main factors in themselves to determine the least fit part for 'survival', and adjust themselves for better fitness."* (seção 3.2)

### 1.2 O operador no AGA

O artigo propõe dois algoritmos:

- **AGA (Adaptive Genetic Algorithm):** acrescenta ao GA convencional o passo 4, *relocation of the worst turbines*:

  > *"Find k wind turbines with lowest power output in each parent's layout and denote it as tᵢ, (i = 1, …, k) and relocate wind turbine tᵢ to a randomly selected locations hosting no wind turbine."* (seção 3.3)

  A mutação aleatória convencional continua existindo no passo 5.

- **SIGA (Self-Informed GA):** faz a mesma realocação, mas escolhe o destino com um modelo substituto (MARS) que estima a potência esperada de cada célula, construído por simulação de Monte Carlo.

Sobre o número de turbinas movidas, os autores fixaram **k = 1**:

> *"Although it is possible to find several of the worst turbines and move them to a better location, in this paper, we limited the number of worst turbines to be 1."* (seção 3.3)

e deixaram k > 1 como trabalho futuro:

> *"The exploration on how to do the relocation of multiple wind turbines in a principled way requires more efforts of investigation and can be our future work."* (seção 6)

### 1.3 Cenário e resultados do artigo

| | Ju e Liu (2019) | Este repositório |
|---|---|---|
| Objetivos | Um: eficiência de conversão de potência | Dois: custo das fundações e potência |
| Posições | Grid discreto de células quadradas (ex.: 21 × 21), turbina no centro da célula | Grid discreto de fundações (~160 m) dentro de polígonos |
| Número de turbinas | Fixo (ex.: 60) | Fixo por zona (`turbines_per_zone`) |
| Áreas | Uma | Até 3 zonas, mais um parque fixo |
| Custo por posição | Não há | Cada fundação tem um custo (depende da profundidade) |
| Distância mínima | Não há (uma turbina por célula) | Não há (uma turbina por fundação) |
| Wake | Jensen | Jensen, com soma quadrática |
| Vento | Distribuições de direção e velocidade (casos 1 a 5) | Um único cenário (10 m/s, 30°) |
| Algoritmo | GA mono-objetivo | MOEA/D |

O AGA difere do GA convencional **só** pelo passo de realocação, então a comparação entre os dois isola o efeito do operador. No caso 1, parque médio, 30 execuções de 200 gerações, a eficiência média foi de 85,6% no CGA, 87,9% no AGA e 88,2% no SIGA. Os outros casos mostram a mesma tendência.

Limites dessa evidência para o nosso uso:
- o problema do artigo é mono-objetivo (só potência);
- a comparação é por número de gerações, e não por número de avaliações.

---

## 2. Adaptação ao nosso problema

### 2.1 O que mantemos do Ju

| Elemento do artigo | No nosso operador |
|---|---|
| Mover a pior turbina, e não uma qualquer | Mantido |
| **k = 1** | Mantido |
| Pior em potência = menor potência gerada, com o wake | Mantido |
| Destino sorteado entre posições livres (versão AGA) | Mantido |

### 2.2 O que adaptamos

| Ponto | Ju e Liu (2019) | Nossa adaptação | Motivo |
|---|---|---|---|
| Objetivo do gargalo | Só potência | Potência **ou** custo, um por chamada | O problema é biobjetivo. O gargalo de custo é a turbina na fundação mais cara. |
| Escolha do objetivo | Não se aplica | Probabilidade dada pelo peso λᵢ do subproblema do MOEA/D | Não cria hiperparâmetro (seção 3.2) |
| Escopo | O parque inteiro | Uma zona sorteada, como a mutação atual | Com o parque inteiro, a zona mais funda receberia sempre o foco em custo. Turbinas não trocam de zona. |
| Destino | Célula livre sorteada (AGA) ou guiada por MARS (SIGA) | Fundação livre sorteada **na mesma zona** | Versão AGA. O SIGA exige um modelo substituto treinado por instância e fica fora do escopo. |
| Onde entra no algoritmo | Passo extra sobre os pais, além da mutação aleatória | **Alternativa** à mutação atual, escolhida por argumento | Permite comparar `random` contra `worst` com o mesmo orçamento. |
| Turbinas fixas | Não há | Nunca são escolhidas | Não podem ser movidas. Entram só no cálculo do wake. |

### 2.3 Inconsistências corrigidas da versão anterior deste documento

A primeira versão deste documento, gerada por uma ferramenta de busca de literatura, tinha pontos que não valem para o nosso problema:

- **Distância mínima σ:** não existe no nosso problema (README, seção 4.4). Todas as checagens de σ foram retiradas.
- **Fórmula de potência `f₂ = ΣΣ Eᵢⱼ zᵢ zⱼ` e AEP:** não é o nosso modelo. A potência é a soma, sobre todas as turbinas, de `0,6127 · P_curva · v³`, com o déficit de velocidade combinado por soma quadrática e um único cenário de vento.
- **Opção 2 (k = 2, pior em potência e pior em custo no mesmo passo):** descartada. Não está "alinhada ao texto do Ju": os autores usaram k = 1 e deixaram k > 1 como trabalho futuro. Além disso, mover duas turbinas por objetivos diferentes mistura os efeitos e torna a perturbação maior com a mesma avaliação.
- **Chamar o operador duas vezes (custo e depois potência):** descartado, porque equivale à opção 2 e gasta duas avaliações.
- **α = 0,5 fixo:** substituído pelo peso λᵢ do subproblema.
- **Código de referência:** dependia de uma `E_matrix` e de uma `InstanceInfo` que não existem no projeto, e não avaliava a solução nem contava a avaliação. O código real seguirá a estrutura de `mutation.cpp`.

---

## 3. Algoritmo proposto

A interface é `void worst_turbine_mutation(Solution &solution, double cost_probability)`. No MOEA/D, ao processar o subproblema i, a mutação recebe `cost_probability = λᵢ₁`, o peso do custo.

1. **Zona:** sorteia uma zona z.
2. **Foco:** com probabilidade `cost_probability`, o foco é o custo. Caso contrário, é a potência.
3. **Pior turbina na zona z:**
   - **custo:** a turbina cuja fundação tem o maior custo;
   - **potência:** a turbina com a menor potência gerada, considerando o wake de todas as turbinas, móveis e fixas.
4. **Caso sem o que corrigir:** se o foco é potência e a pior turbina gera a potência de vento livre (não sofre wake), não há gargalo de potência. O foco passa para o custo.
5. **Destino:** sorteia uma fundação livre da zona z.
6. **Movimento:** move a turbina para o destino, atualizando `grid`, `index`, `x` e `y`. O id é mantido.
7. **Avaliação:** calcula custo e potência, insere no arquivo Pareto, incrementa `countRevalue` e salva os snapshots, como a mutação atual.

### 3.1 Exemplo

Subproblema i = 80, com λ₈₀ ≈ (0,81; 0,19):
- em cerca de 81% das chamadas, a turbina na fundação mais cara da zona sorteada vai para uma fundação livre;
- em cerca de 19%, a turbina que gera menos potência vai para uma fundação livre.

### 3.2 Por que o peso λ do subproblema

O vetor λᵢ = (i/99, 1 − i/99) diz o quanto o subproblema i valoriza o custo e a potência. O filho gerado na iteração i só substitui soluções dos vizinhos B(i), que têm pesos parecidos. Usar λᵢ₁ como probabilidade do foco faz a mutação tender a melhorar o objetivo que importa para quem vai receber o filho:

| Subproblema | λ (custo, potência) | Foco em custo |
|---|---|---|
| 0 | (0; 1) | 0% |
| 50 | (0,51; 0,49) | 51% |
| 99 | (1; 0) | 100% |

Uma probabilidade fixa α seria um hiperparâmetro a calibrar (por exemplo, com o IRACE). Com λᵢ, a probabilidade vem do próprio MOEA/D.

### 3.3 Garantias

- **Número de turbinas por zona exato:** uma turbina sai de uma fundação e entra em outra da mesma zona.
- **Uma turbina por fundação:** o destino é sempre uma fundação livre.
- **Ids únicos:** o id da turbina movida é mantido.
- **O layout sempre muda:** a turbina vai para uma fundação diferente, então não há filho idêntico e toda chamada conta 1 avaliação, como a mutação atual.

---

## 4. Hiperparâmetros

| Parâmetro | Valor | Observação |
|---|---|---|
| k (turbinas movidas por chamada) | 1 | Fixo, como no Ju |
| Probabilidade do foco em custo | λᵢ₁ do subproblema | Não é escolhido: vem do MOEA/D |
| `mutation_type` | `random` (padrão) ou `worst` | Escolha do operador (categórico) |

O operador não cria hiperparâmetros numéricos novos.

---

## 5. Implementação

### 5.1 Potência por turbina guardada na avaliação

Encontrar a turbina de menor potência exige a potência de cada turbina com o wake. O `calculate_power` já calcula esse valor (`power_produced(windResulted, turbine)`), mas o descarta depois de somar. Recalcular dentro da mutação custaria um cálculo de wake inteiro a mais por chamada.

Esse valor é guardado no campo `generated_power` da `Turbine` durante o `calculate_power`. Assim, encontrar a pior turbina não custa nada. A potência total da solução continua sendo calculada exatamente como antes.

Os valores estão sempre atualizados quando a mutação é chamada, porque toda solução que chega a ela foi avaliada antes:
- população inicial: avaliada em `generate_solution`;
- filho de crossover: avaliado no próprio crossover;
- cópia de um pai (crossover não aplicado ou Block Swap devolvendo a base): o pai já foi avaliado.

O campo `power` da `Turbine` continua sendo o valor da curva de potência e não muda.

### 5.2 Arquivos

| Arquivo | Mudança |
|---|---|
| `src/global_modules/genetic_operators/mutation/worst_turbine_mutation.cpp` | Operador (novo) |
| `headers/global_modules/genetic_operators/worst_turbine_mutation.h` | Header (novo) |
| `headers/wind_farm_entities/turbine.h` | Campo `generated_power` |
| `src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp` | `calculate_power` guarda a potência de cada turbina móvel |
| `Makefile` | Novo arquivo em `SRC_FILES` |
| `main.cpp` | Global `mutation_type`, leitura de `argv[6]`, validação e impressão no log |
| `src/metaheuristics/moead/moead.cpp` | Escolha da mutação; `worst` recebe `lambda_vector[i].first` |

A mutação atual (`mutation.cpp`) não muda.

### 5.3 Uso

```bash
cd scripts
./exe/moead A ../dataset/moead/A/1/ 30 10 onepoint worst    # Worst Turbine Mutation
./exe/moead A ../dataset/moead/A/1/                         # mutação aleatória (padrão)
```

A mutação é o 6º argumento, então é preciso passar também o ângulo, o vento e o crossover.

### 5.4 Validação

O operador foi testado com o código real nas instâncias A (2 zonas, 34 turbinas) e 240 (1 zona, 77 turbinas). Para cada probabilidade de foco em custo (0; 0,25; 0,5; 0,75; 1), foram 400 chamadas: metade sobre soluções aleatórias e metade sobre soluções aleatórias com 20 mutações aleatórias.

Em todas as chamadas, sem nenhuma falha:
- exatamente uma turbina foi movida, na mesma zona, para uma fundação que estava livre, e manteve o id;
- a turbina movida era sempre a pior em custo ou a pior em potência da zona sorteada;
- o número de turbinas por zona, a coerência entre `grid` e turbinas e a unicidade dos ids se mantiveram;
- `generated_power` bateu com um recálculo independente da potência de cada turbina;
- o custo e a potência total bateram com um recálculo independente (a potência total, feito com a fórmula original, confirma que guardar `generated_power` não alterou o cálculo);
- cada chamada contou exatamente 1 avaliação.

Frequência do foco em custo:

| Probabilidade de foco em custo | A | 240 |
|---|---|---|
| 0 | 2,8% | 1,0% |
| 0,25 | 28,2% | 23,2% |
| 0,5 | 52,2% | 52,0% |
| 0,75 | 75,8% | 76,0% |
| 1 | 100% | 100% |

As frequências acompanham a probabilidade. Com probabilidade 0, os 1% a 3% contados como "custo" são casos em que a mesma turbina era a pior nos dois critérios. O caso sem perda por wake (seção 3, passo 4) não ocorreu nesses testes: com o vento de 30° e o parque fixo, sempre havia alguma turbina sob wake.

O MOEA/D completo rodou até o fim com orçamento reduzido (50 mil avaliações) com as duas mutações, e um valor inválido no 6º argumento encerra o programa com erro. Uma rodada curta não permite comparar o desempenho das mutações.

---

## 6. Próximos passos

- **Comparação experimental:** `random` contra `worst`, 20 execuções por instância, em instâncias de 1 e 2 zonas, comparando o hipervolume. A comparação com o Block Swap fica separada, para isolar o efeito de cada operador.
- **SIGA:** guiar o destino por um mapa de potência esperada por fundação, como o MARS do Ju. Com um único cenário de vento, esse mapa poderia ser calculado de forma mais simples. Fica para depois.
- **Destino guiado pelo custo:** no foco em custo, sortear o destino só entre fundações mais baratas que a atual. É um desvio do AGA e precisaria de justificativa.

---

## 7. Texto sugerido para a metodologia

> *"Além da mutação aleatória, que move uma turbina sorteada para uma fundação livre, foi implementado um operador de mutação baseado na realocação da pior turbina, adaptado do Adaptive Genetic Algorithm de Ju e Liu (2019). Os autores propõem, em um problema mono-objetivo, realocar a turbina de menor potência para uma posição livre sorteada, com uma turbina por vez. Para o problema biobjetivo considerado, o operador escolhe, a cada aplicação, entre a turbina de menor potência e a turbina instalada na fundação mais cara de uma zona sorteada. O foco em custo é escolhido com probabilidade igual ao peso do custo no subproblema do MOEA/D que gera o descendente, o que direciona a mutação ao objetivo priorizado pela vizinhança, sem introduzir parâmetros adicionais."*

O operador deve ser apresentado como uma **adaptação biobjetivo** do Ju: a evidência do artigo é mono-objetivo, e o efeito no nosso problema precisa ser demonstrado experimentalmente.

---

## 8. Referências

[1] X. Ju, F. Liu. *Wind farm layout optimization using self-informed genetic algorithm with information guided exploitation*. Applied Energy 248, pp. 429–445, 2019. https://doi.org/10.1016/j.apenergy.2019.04.084
