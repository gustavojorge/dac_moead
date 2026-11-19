# dac_moead — MOEA/D para o Wind Farm Layout Optimization Problem (WFLOP) multiobjetivo

Este repositório contém a implementação em C++ do **MOEA/D** (Multi-Objective Evolutionary Algorithm based on Decomposition) aplicada ao problema de **posicionamento de turbinas em parques eólicos offshore** com dois objetivos: **minimizar o custo de fundação** e **maximizar a potência gerada**.

Também contém os scripts shell para compilar e executar os experimentos em lote (várias instâncias × 20 execuções).

O código vem da branch [`moead-runner`](https://github.com/gustavojorge/MO_WFLOP-experiment-runner/tree/moead-runner) do MO_WFLOP-experiment-runner: só o MOEA/D, sem busca local e sem as versões híbridas.

---

## Sumário

1. [Estrutura do repositório](#1-estrutura-do-repositório)
2. [Como o MOEA/D funciona (e como está implementado)](#2-como-o-moead-funciona-e-como-está-implementado)
3. [Como executar](#3-como-executar)
4. [Modelagem do parque eólico](#4-modelagem-do-parque-eólico)
5. [Operadores genéticos](#5-operadores-genéticos)
6. [Hiperparâmetros utilizados](#6-hiperparâmetros-utilizados)
7. [Problemas conhecidos](#7-problemas-conhecidos)
8. [Referência de comandos do supercomputador](#8-referência-de-comandos-do-supercomputador)

---

## 1. Estrutura do repositório

```
.
├── makeall.sh                     # Compila src/moead e move o executável para scripts/exe/moead
├── scripts/
│   ├── batch.sh                   # Dispara lotes de instâncias em background (nohup)
│   ├── main.sh                    # Para cada instância do lote, roda 20 execuções
│   └── smaller_scripts/
│       ├── create_folders.sh      # Cria dataset/moead/<instância>/<run>/
│       └── moead.sh               # Chama ./exe/moead para (instância, run) e grava log.txt
├── dataset/
│   ├── moead/<inst>/<run>/        # Saída das execuções (criada pelos scripts, não versionada)
│   └── reference/moead/<inst>/<run>/  # Resultados de referência (A, 172, 240, 456), ver seção 3.5
└── src/moead/
    ├── main.cpp                   # Ponto de entrada: lê a instância e chama moead()
    ├── instance_info.cpp          # Leitura da instância, curvas das turbinas, turbinas fixas
    ├── Makefile
    ├── headers/                   # .h espelhando a árvore de src/
    ├── metafeatures/              # Extração de metafeatures (não compilado)
    └── src/
        ├── metaheuristics/
        │   ├── moead/             # moead.cpp e modules/
        │   └── nsga2/             # compilado, mas não chamado pelo main
        ├── global_modules/
        │   ├── generate_initial_population/   # Solução aleatória + funções objetivo + wake effect
        │   ├── genetic_operators/             # crossover/ e mutation/
        │   ├── dominates.cpp, isEqual.cpp
        ├── adaptive_grid/         # Arquivo Pareto limitado com grid adaptativo
        ├── instance_generation/   # Gerador de instâncias sintéticas (Python)
        └── utils/                 # Validador, comparação e plot (não compilados)
```

O executável contém apenas o MOEA/D. O `nsga2.cpp` entra na compilação, mas o `main.cpp` não o chama.

---

## 2. Como o MOEA/D funciona (e como está implementado)

### 2.1 Ideia geral

O MOEA/D transforma o problema multiobjetivo em **N subproblemas escalares**. Cada subproblema *i* tem um vetor de pesos λᵢ = (λᵢ₁, λᵢ₂) e é resolvido pela **função de Tchebycheff**:

```
g_te(x | λᵢ, z*) = max_k { λᵢk · |f_k(x) − z*_k| }
```

onde `z*` é o ponto de referência (ideal). Os subproblemas com pesos próximos têm soluções ótimas próximas, então cada subproblema troca informação só com os seus **T vizinhos** (os vetores de peso mais próximos). Cada subproblema guarda uma solução; a população é o conjunto dessas N soluções.

Neste código **os dois objetivos são maximizados**: o custo é guardado negativo (`fitness.first = −custo`) e a potência positiva (`fitness.second = potência`). Por isso `z*` é o **máximo** de cada objetivo.

### 2.2 Fluxo do algoritmo

Implementação: [moead.cpp](src/moead/src/metaheuristics/moead/moead.cpp).

```
1. Inicialização
   1.1 População inicial: N = 100 soluções aleatórias (create_initial_population)
       → cada uma é avaliada (100 avaliações) e enviada ao arquivo Pareto global
   1.2 Vetores de peso λ uniformes: λᵢ = (i/(N−1), 1 − i/(N−1))
   1.3 Vizinhança B(i): os T = 10 vetores mais próximos de λᵢ (distância euclidiana, sem incluir o próprio i)
   1.4 z* = (max f1, max f2) da população inicial

2. Enquanto avaliações < 1.000.000 (checado a cada geração):
   Para cada subproblema i = 0..N−1:
     2.1 Seleciona dois pais distintos k, l ∈ B(i)
     2.2 Com prob. 0.5: filho1 = crossover(k, l), filho2 = crossover(l, k)   (2 avaliações)
         Com prob. 0.5: mutação em filho1 e filho2                         (2 avaliações)
         (se nenhum operador for aplicado, os filhos são cópias dos pais e não há avaliação)
     2.3 Atualiza z* com o máximo de cada objetivo no arquivo Pareto global
     2.4 Para cada j ∈ B(i): se g_te(filho | λⱼ, z*) ≤ g_te(xⱼ | λⱼ, z*), xⱼ ← filho
         (feito para filho1 e depois para filho2; sem limite de substituições)
```

B(i) guarda apenas os índices dos vizinhos. Quando o filho vence, ele é **copiado** para `population[j]`: a solução do subproblema j muda, e o peso λⱼ continua o mesmo.

### 2.3 Módulos principais

| Módulo | Arquivo | Função |
|---|---|---|
| Vetores de peso | [generate_weight_vectors.cpp](src/moead/src/metaheuristics/moead/modules/generate_weight_vectors.cpp) | Gera N pares (w1, w2) igualmente espaçados com w1 + w2 = 1 |
| Vizinhança | [generate_neighborhood.cpp](src/moead/src/metaheuristics/moead/modules/generate_neighborhood.cpp) | Para cada λᵢ, ordena os demais por distância euclidiana e pega os T primeiros |
| Ponto de referência | [get_best_z_point.cpp](src/moead/src/metaheuristics/moead/modules/get_best_z_point.cpp) | z* = máximo de cada objetivo na população |
| Tchebycheff | [tchebycheff.cpp](src/moead/src/metaheuristics/moead/modules/tchebycheff.cpp) | `calculate_gte`: max λk·\|fk − z*k\| (sem normalização dos objetivos) |
| População externa (EP) | [updateEP.cpp](src/moead/src/metaheuristics/moead/modules/updateEP.cpp) | Implementação clássica da EP — **não usada**; o papel de EP é do arquivo Pareto global |
| Arquivo Pareto global | [BoundedParetoSet.cpp](src/moead/src/adaptive_grid/BoundedParetoSet.cpp), [ParetoSet.cpp](src/moead/src/adaptive_grid/ParetoSet.cpp) | Guarda até 300 soluções não dominadas; quando cheio, remove uma solução da célula mais populosa de um grid adaptativo (32×32 células) |

### 2.4 Arquivo Pareto global (`pareto`)

Variável global `BoundedParetoSet * pareto` declarada em [main.cpp](src/moead/main.cpp). **Toda avaliação** (geração aleatória, crossover, mutação) chama `pareto->adicionarSol(...)`. É esse arquivo que é salvo como resultado do algoritmo, e não a população do MOEA/D.

---

## 3. Como executar

### 3.1 Compilação

```bash
./makeall.sh
```

Roda `make rebuild` em `src/moead` (g++ com `-std=c++17 -Werror -O2`), move o executável para `scripts/exe/moead` e dá permissão de execução aos scripts.

### 3.2 Estrutura de instâncias esperada

As instâncias ficam no repositório [MO_WFLOP-archive/instances](https://github.com/gustavojorge/MO_WFLOP-archive/tree/main/instances) e não são versionadas aqui.

O executável lê as instâncias de **`../instances`, relativo ao diretório de trabalho** (`pathFolders = ".."` em [instance_info.cpp](src/moead/instance_info.cpp)). Como os scripts executam `./exe/moead` a partir de `scripts/`, a pasta deve ficar na **raiz do repositório**.

> **Atenção:** no MO_WFLOP-archive a pasta se chama `instances/sites/`, mas o código procura `instances/site/` (sem "s"). Depois de copiar, crie um link: `ln -s sites instances/site`.

```
instances/                      # na raiz do repositório (não versionada)
├── site/                       # (= sites/ no MO_WFLOP-archive)
│   ├── A/ ... J/               # instâncias de referência
│   ├── 0/ ... 505/             # instâncias sintéticas
│   │   ├── availablePositions.txt  # lido: fundações candidatas das zonas
│   │   ├── fixed_wf.txt            # lido: turbinas fixas
│   │   ├── turbines_per_zone.txt   # lido: só para instâncias numéricas
│   │   ├── geometry.txt            # não lido: contorno dos polígonos das zonas
│   │   ├── info.json               # não lido: metadados (só nas instâncias A–J)
│   │   └── plot.png                # não lido: visualização da instância
├── wtg/
│   ├── NREL-15-240.txt         # curva da turbina móvel (15 MW, rotor 240 m)
│   └── NREL-10-179.txt         # curva da turbina fixa (10 MW, rotor 179 m)
└── wind/
    └── RVO_*.txt               # rosas dos ventos (não lidas pelo código atual)
```

Formato dos arquivos (separados por espaço ou tab, **um registro por linha, terminando com quebra de linha**):

| Arquivo | Colunas | Observação |
|---|---|---|
| `availablePositions.txt` | `x y profundidade custo zona` | Fundações candidatas. Zona começa em 1 e o arquivo deve estar ordenado por zona: o número de zonas é lido da última linha. Máximo de 3 zonas. |
| `fixed_wf.txt` | `x y 0.0 z z` | Turbinas já existentes (fixas). Só `x` e `y` são usados; `z` é o número da "zona" fixa (ex.: 3). |
| `turbines_per_zone.txt` | `n1 n2 n3` | Quantidade de turbinas a posicionar por zona. Para A–J os valores estão fixos no código (`turbines_instace`). |
| `geometry.txt` | `x y` por vértice, linha em branco entre polígonos | Contorno externo de cada zona. |
| `wtg/*.txt` | `vento(m/s) potência(MW) Ct` | Curva de potência e coeficiente de empuxo por velocidade de vento. |
| `wind/RVO_*.txt` | `ângulo velocidade probabilidade` | Distribuição do vento em 12 setores de 30° e velocidades de 0 a 30 m/s. |

A seção [4.1](#41-instâncias-zonas-e-turbinas-fixas) explica o que esses dados representam, com exemplos.

As instâncias sintéticas são geradas por [script.py](src/moead/src/instance_generation/script.py) (detalhes na seção 4.1). O script cria `../instances/<N>/`; mova para `instances/site/<N>/`.

### 3.3 Linha de comando do executável

```bash
cd scripts
./exe/moead <instância> <pasta_saida/> [<ângulo> <vento> [<crossover> [<mutação>]]]
```

| Arg | Exemplo | Descrição |
|---|---|---|
| `instância` | `A`, `240` | Nome da pasta em `instances/site/` |
| `pasta_saida/` | `../dataset/moead/A/1/` | Prefixo de saída — **precisa da barra final** (padrão `./`) |
| `ângulo` | `330` | Direção do vento em graus. Só é lido se `vento` também for passado. Padrão: **30** |
| `vento` | `13` | Velocidade do vento em m/s. Padrão: **10** |
| `crossover` | `blockswap` | Operador de crossover: `onepoint` (padrão), `blockswap` ou `blockcopy` (ver [seção 5](#5-operadores-genéticos)). Só é lido se `ângulo` e `vento` também forem passados; outro valor encerra com erro |
| `mutação` | `worst` | Operador de mutação: `random` (padrão) ou `worst` (ver [seção 5](#5-operadores-genéticos)). Só é lido se os quatro argumentos anteriores forem passados; outro valor encerra com erro |

Exemplos com o cenário de vento padrão: `./exe/moead A ../dataset/moead/A/1/ 30 10 blockswap` (Block Swap) e `./exe/moead A ../dataset/moead/A/1/ 30 10 onepoint worst` (Worst Turbine Mutation).

O `moead.sh` passa só os dois primeiros argumentos, então os experimentos em lote rodam sempre com **vento de 10 m/s, ângulo de 30°, crossover `onepoint` e mutação `random`** (ver [problemas conhecidos](#7-problemas-conhecidos)).

### 3.4 Execução em lote

```bash
cd scripts
./batch.sh
```

- `batch.sh` roda `create_folders.sh`, cria `scripts/logs/` e dispara **52 processos em background** com `nohup`: um para A–J + 0, cinquenta lotes de 10 instâncias (1–500) e um para 501–505. Os logs vão para `scripts/logs/batch_*.txt`.
- `main.sh` roda **20 execuções** de cada instância do lote, chamando `smaller_scripts/moead.sh <instância> <run>`, e acrescenta a saída a `logs/<instância>.txt`.
- `moead.sh` executa `./exe/moead <instância> ../dataset/moead/<instância>/<run>/` e grava a saída, com o tempo (`time`), em `log.txt` na mesma pasta.
- `create_folders.sh` cria `dataset/moead/{A..J,0..505}/{1..20}/`, mas não faz nada se `dataset/moead` já existir (por exemplo, numa segunda rodada). Como `moead.sh` faz `mkdir -p`, isso não impede a execução.

### 3.5 Saídas

Em `dataset/moead/<inst>/<run>/`:

| Arquivo | Conteúdo |
|---|---|
| `log.txt` | Escrito pelo `moead.sh`: instância, run, nº de subproblemas, vizinhos, turbinas fixas e móveis, vento, ângulo, crossover e tempo (`time`) |
| `infoRun.txt` | Por geração: avaliações acumuladas e tamanho do arquivo Pareto |
| `<inst>_moead_<N>.txt` | Arquivo Pareto a cada 100.000 avaliações (`N = 100000, …, 1000000`). Cada linha: `custo potência` (custo em valor absoluto) |
| `<inst>_moead_layout.txt` | Gerado na avaliação 1.000.000: coordenadas `x y` das turbinas móveis de cada solução do arquivo, com uma linha em branco entre soluções |

**Resultados de referência:** `dataset/reference/moead/` tem as 20 execuções das instâncias A, 172, 240 e 456, com a mesma estrutura de arquivos. Eles foram gerados por uma versão anterior, que **sorteava o vento da rosa dos ventos** (ex.: A/run 1 → 13 m/s, 330°), e não são comparáveis com execuções do código atual, que usa 10 m/s e 30°. Ficam fora de `dataset/moead/` para que o lote não os sobrescreva.

---

## 4. Modelagem do parque eólico

Código: [generate_rSolution.cpp](src/moead/src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp), [instance_info.cpp](src/moead/instance_info.cpp), [headers/wind_farm_entities/](src/moead/headers/wind_farm_entities/).

### 4.1 Instâncias: zonas e turbinas fixas

Uma instância descreve um sítio offshore com **até 3 zonas**. As coordenadas estão em metros, numa faixa de ~100–150 km. Há dois tipos de zona:

- **Zonas disponíveis** (1, 2, …): polígonos onde o algoritmo pode instalar turbinas novas.
- **Zona fixa** (o número seguinte, ex.: 3): um parque já existente, fora dos polígonos. Suas turbinas não se movem, mas entram no cálculo do wake e da potência.

![Instância 10](https://raw.githubusercontent.com/gustavojorge/MO_WFLOP-archive/main/instances/sites/10/plot.png)

*Instância 10: duas zonas disponíveis (contorno magenta, cor = profundidade em m) e o parque fixo (pontos laranja).*

**Como uma zona é construída** (gerador [script.py](src/moead/src/instance_generation/script.py) e [gen_layout.py](src/moead/src/instance_generation/gen_layout.py)):

1. **Polígono**: de 3 a 9 vértices aleatórios, escalado para ter entre 20 e 550 km² (a soma das zonas fica entre 100 e 1500 km²). Zonas que se sobrepõem são reposicionadas. O contorno vai para `geometry.txt`.
2. **Áreas proibidas**: de 2 a 5 **buracos** (polígonos convexos dentro da zona) e **estruturas** (faixas que cortam a zona a partir da borda, como cabos ou dutos), num total de 5 a 8 recortes. Na figura, são as reentrâncias e os buracos brancos.
3. **Fundações candidatas**: a zona recortada é coberta por um **grid regular de ~160 m**. Cada centro de célula dentro da zona vira uma fundação em `availablePositions.txt`.
4. **Profundidade**: gerada por Perlin noise, em torno de −25 m (−32 a −19 m na instância 10).
5. **Custo da fundação**: linear na profundidade, `custo = 659934 − 72606.6 · profundidade`. Quanto mais fundo, mais caro (~2,0 M a ~3,0 M por fundação na instância 10).
6. **Turbinas a instalar por zona**: `int(U(0.005, 0.02) · nº de fundações da zona)`, salvo em `turbines_per_zone.txt`.
7. **Parque fixo**: de 5 a 75 turbinas num grid retangular com **1440 m** entre turbinas, colado a uma aresta de uma zona sorteada e sem invadir nenhuma zona. Fica em `fixed_wf.txt`.

**Exemplos reais**

| Instância | Zona 1 | Zona 2 | Turbinas fixas |
|---|---|---|---|
| `10` | 11.182 fundações, 63 turbinas | 8.335 fundações, 82 turbinas | 64 |
| `A` | 50,4 km², 2 buracos, 2.187 fundações, 26 turbinas | 26,8 km², 4 buracos, 1.009 fundações, 8 turbinas | 42 |

As instâncias A–J têm um `info.json` com área, número de buracos e opções de quantidade de turbinas para densidades de potência de 4, 6 e 8 MW/km² (ex.: zona 1 de A: 13, 20 ou 26 turbinas de 15 MW). As quantidades usadas estão fixas em `turbines_instace` ([instance_info.cpp](src/moead/instance_info.cpp#L26)).

**Como o código usa esses dados**

- `availablePositions.txt` → `foundations[z]` (`Foundation {x, y, cost}`). A posição da fundação no arquivo, dentro da zona, é o índice `p` usado em `grid[z][p]`.
- `fixed_wf.txt` → `fixd` (turbinas NREL 10 MW).
- `geometry.txt`, `info.json` e `plot.png` **não** são lidos: os polígonos já estão representados pelas fundações que sobraram dentro deles.
- O grid de 160 m é menor que o diâmetro do rotor (240 m), então duas fundações vizinhas podem receber turbinas. Não há distância mínima, e a proximidade só é penalizada pelo wake.

### 4.2 Turbinas e vento (dados do wake effect)

- **Turbinas móveis**: as que o algoritmo posiciona. Modelo NREL 15 MW: diâmetro 240 m, altura do cubo 150 m, curva `wtg/NREL-15-240.txt`.
- **Turbinas fixas** (`fixd`): modelo NREL 10 MW: diâmetro 179 m, altura 119 m, curva `wtg/NREL-10-179.txt`.
- **Curvas** (`vento potência Ct`): para a velocidade do vento livre `U`, o código lê a potência (`P_curva`, em MW) e o coeficiente de empuxo `Ct` de cada modelo. Se `U` não estiver na tabela, usa a média entre as linhas vizinhas. Esses valores são os mesmos para todas as turbinas do mesmo modelo, pois são lidos com o vento livre, antes do wake.
- **Vento**: um único cenário por execução, com velocidade `wind` (m/s) e direção `angle` (graus; o vetor do vento é `(cos θ, sin θ)`). Vem da linha de comando ou dos valores padrão, **10 m/s e 30°**, que são os usados pelo `moead.sh`. Os arquivos `wind/RVO_*.txt` (rosas dos ventos de sítios offshore holandeses, ex.: HKN = Hollandse Kust Noord) têm a distribuição `ângulo velocidade probabilidade`, mas o sorteio a partir deles está comentado em `instance_info.cpp`.

### 4.3 Solução

```cpp
struct Solution {
    vector<vector<int>> grid;          // grid[z][p] = id da turbina na fundação p da zona z, 0 = vazia
    vector<vector<Turbine>> turbines;  // turbines[z] = turbinas da zona z (com x, y, índice da fundação)
    pair<double,double> fitness;       // (−custo, potência)
};
```

Os ids das turbinas móveis começam em `fixd.size()`, para não colidir com os ids das fixas (0..F−1) nem com o 0 usado como "vazio". Os ids precisam ser únicos na solução, porque o wake ignora pares com o mesmo id.

A solução inicial é aleatória: para cada zona, sorteia `turbines_per_zone[z]` fundações distintas.

### 4.4 Restrições do layout

Todas são garantidas **pela representação e pelos operadores** (não há penalização nem soluções inviáveis):

1. Turbinas só podem ficar em **fundações candidatas** da sua zona (posições discretas; áreas proibidas — buracos e estruturas — já foram removidas na geração da instância).
2. **Número exato de turbinas por zona** (`turbines_per_zone`) — crossover e mutação preservam essa cardinalidade.
3. **No máximo uma turbina por fundação**.
4. Turbinas **não mudam de zona**: os operadores agem zona a zona.
5. Turbinas fixas não podem ser alteradas.

Não existe restrição explícita de distância mínima entre turbinas; a proximidade é penalizada indiretamente pelo efeito esteira.

### 4.5 Funções objetivo

**f1 — custo** (`calculate_cost`): soma do custo das fundações ocupadas, guardada negativa para ser maximizada.

```
f1(x) = − Σ custo(fundação de cada turbina móvel)
```

**f2 — potência** (`calculate_power`): soma da potência de **todas** as turbinas (móveis + fixas), cada uma com o vento reduzido pelas esteiras das demais:

```
P_i = 0.6127 · P_curva · v_i³
```

### 4.6 Efeito esteira (wake effect)

Modelo do tipo Jensen/Park ([calculate_interference](src/moead/src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp#L46)). Para a turbina *a* gerando esteira sobre a turbina *b*:

1. **Está a jusante?** Projeta o vetor `b − a` na direção do vento `(cos θ, sin θ)`. Se a projeção `s ≤ 0`, *b* está atrás de *a* em relação ao vento → sem interferência.
2. **Coeficiente de decaimento**: `α = 0.5 / ln(h_a / z0)`, com `z0 = 0.0005` (rugosidade offshore).
3. **Diâmetro da esteira** a uma distância `s` a jusante: `D_w = D_a + α · s`.
4. **Está dentro da esteira?** Se `(D_w + D_b)/2 ≤ distância de b ao eixo da esteira` → sem interferência.
5. **Déficit de velocidade**:
   ```
   δ_ab = (1 − √(1 − Ct_a)) / (1 + α · d_ab / (D_b / 2))²
   ```
   onde `d_ab` é a distância euclidiana entre as turbinas.
6. **Combinação** (soma quadrática de todas as esteiras sobre *b*, vindas de móveis e fixas):
   ```
   v_b = U · (1 − min(√Σ δ²_ab, 1))
   ```

A potência total exige comparar todos os pares de turbinas: custo O((n + F)²) por avaliação.

### 4.7 Orçamento de avaliações

- **Máximo**: `stop_criteria = 1.000.000` avaliações ([main.cpp](src/moead/main.cpp#L11)). O contador global é `countRevalue`.
- **O que conta como 1 avaliação**: cada chamada de `calculate_cost` + `calculate_power` sobre uma nova solução, isto é, cada solução da população inicial, cada filho de crossover e cada mutação. No Block Swap e no BlockCopy, o filho que sai com o mesmo layout do pai-base não é avaliado nem contado.
- **Como é gasto**:
  - 100 avaliações na população inicial;
  - por subproblema: 0, 2 ou 4 avaliações (crossover com prob. 0.5 gera 2 filhos; mutação com prob. 0.5 muta os 2) → em média **~200 avaliações por geração**, ou seja, cerca de **5.000 gerações**.
- **Parada**: o laço só verifica o critério no início de cada geração, então pode ultrapassar 1.000.000 em até ~400 avaliações.
- **Snapshots**: a cada 100.000 avaliações o arquivo Pareto é salvo em disco; na avaliação 1.000.000 também é salvo o layout.

---

## 5. Operadores genéticos

Arquivos: [crossover.cpp](src/moead/src/global_modules/genetic_operators/crossover/crossover.cpp), [block_swap_crossover.cpp](src/moead/src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp), [block_copy_crossover.cpp](src/moead/src/global_modules/genetic_operators/crossover/block_copy_crossover.cpp), [mutation.cpp](src/moead/src/global_modules/genetic_operators/mutation/mutation.cpp) e [worst_turbine_mutation.cpp](src/moead/src/global_modules/genetic_operators/mutation/worst_turbine_mutation.cpp). Todos avaliam a solução gerada, inserem no arquivo Pareto global e incrementam o contador de avaliações.

O crossover usado é escolhido pelo 5º argumento da linha de comando (`onepoint`, `blockswap` ou `blockcopy`, ver seção 3.3). Nos três casos, o MOEA/D gera dois filhos: `crossover(A, B)` e `crossover(B, A)`, e o filho é numerado com ids únicos em todo o parque.

**Crossover de um ponto (`onepoint`, padrão)** — para cada zona:
1. Sorteia um ponto de corte `c` sobre o vetor de fundações.
2. Copia do pai A todas as turbinas das posições `0..c`.
3. Percorre o pai B **do fim para o começo**, adicionando as turbinas cujas posições ainda estão livres no filho, até atingir `turbines_per_zone[z]`.

O resultado sempre tem o número correto de turbinas por zona. Como o `availablePositions.txt` lista as fundações em varredura do grid (por coluna na instância A, por linha nas instâncias 240 e 10), o corte `0..c` separa uma **faixa contígua** do parque: o filho é, grosso modo, um lado de A colado ao outro lado de B.

**Block Swap (`blockswap`)** — adaptado de Tran et al. (GECCO 2013), que o usam com coordenadas contínuas e distância mínima entre turbinas. Em vez de uma faixa, troca um **retângulo** do parque. O primeiro pai é a base e o segundo é o doador. Para cada zona:
1. Sorteia um retângulo R centrado numa turbina do doador, com largura e altura entre 10% e 50% das dimensões da zona (retângulo envolvente das fundações).
2. Separa as turbinas da base fora de R (`O_A`), as da base dentro de R (`I_A`) e as do doador dentro de R (`I_B`).
3. Aceita R só se `|I_B| ≥ |I_A|`; senão, sorteia outro R, até 10 tentativas. Se nenhuma servir, a zona fica igual à da base.
4. O filho recebe `O_A` mais `I_B`. As fundações de dentro e de fora de R são disjuntas, então não há conflito.
5. O excesso (`|I_B| − |I_A|` turbinas) é retirado ao acaso entre as de `O_A`, para manter o bloco do doador intacto e o número exato de turbinas.

Toda turbina do filho vem de um dos pais (nenhuma posição aleatória nova) e fica na mesma fundação, com o mesmo custo. Se o filho sair com o mesmo layout da base, ela é devolvida sem nova avaliação. Diferente do artigo original, não há distância mínima, porque o problema não tem essa restrição.

**BlockCopy (`blockcopy`)** — adaptado do BlockCopy crossover de Mayo e Zheng (IEEE, 2015), que é mono-objetivo, contínuo e usa blocos fixos de 1 km × 1 km. Em vez de trocar a mesma região, **copia um padrão de turbinas do doador para outro lugar** da base, mantendo a arrumação relativa entre elas. O primeiro pai é a base e o segundo é o doador. Para cada zona:
1. Sorteia a origem: um retângulo no doador, centrado numa turbina do doador, com largura e altura entre 10% e 50% das dimensões da zona (como no Block Swap).
2. Sorteia o destino: uma fundação qualquer da zona, diferente do centro da origem. O deslocamento é Δ = destino − centro da origem, em x e em y.
3. Remove da base as turbinas dentro do retângulo de destino (mesmo tamanho, centrado no destino).
4. Cola as turbinas do doador que estão na origem, cada uma deslocada por Δ. As fundações de cada zona formam um grid regular, então a posição deslocada cai exatamente sobre uma fundação. Se ela não existir (fora da zona) ou estiver ocupada, a turbina é descartada.
5. Acerta a quantidade como Mayo e Zheng: o excesso é removido ao acaso entre as turbinas que vieram da base (o padrão colado fica intacto); o déficit é preenchido em fundações livres sorteadas na zona.

Δ ≠ 0 garante que o BlockCopy nunca repete o Block Swap. Diferente do Block Swap, o filho pode ter turbinas em fundações que não estavam em nenhum dos pais (as deslocadas e as do preenchimento) e, como o custo depende da fundação, o padrão colado preserva a configuração de wake, mas não o custo. Se o filho sair com o mesmo layout da base, ela é devolvida sem nova avaliação.

A mutação usada é escolhida pelo 6º argumento da linha de comando (`random` ou `worst`, ver seção 3.3). As duas movem **uma** turbina de uma zona sorteada para uma fundação livre da mesma zona, o que preserva a cardinalidade. A diferença está em qual turbina é movida.

**Mutação aleatória (`random`, padrão)**: sorteia uma posição ocupada e uma livre, e move a turbina para a posição livre.

**Worst Turbine Mutation (`worst`)** — adaptada da *relocation of the worst turbines* do AGA de Ju e Liu (Applied Energy, 2019), que é mono-objetivo e move a turbina de menor potência. Na iteração do subproblema i do MOEA/D:
1. Sorteia uma zona.
2. Com probabilidade λᵢ₁ (o peso do custo no subproblema), o foco é o **custo**; caso contrário, é a **potência**. Não há hiperparâmetro novo.
3. Escolhe a pior turbina da zona: a de fundação mais cara (custo) ou a de menor potência gerada, com o wake (potência). Se a pior em potência não sofrer wake, o foco passa para o custo.
4. Move essa turbina para uma fundação livre sorteada na mesma zona, mantendo o id.

A potência gerada por cada turbina móvel é guardada em `Turbine::generated_power` durante o `calculate_power`, então achar a pior turbina não custa avaliação extra. Detalhes em [docs/worst_turbine_mutation](docs/worst_turbine_mutation/operador-worst-turbine-mutation.md).

`mutationElisama.cpp` e `mutationFeature.cpp` são variantes antigas (usadas na extração de metafeatures) e não são compiladas.

---

## 6. Hiperparâmetros utilizados

| Parâmetro | Valor | Onde |
|---|---|---|
| Nº de subproblemas / tamanho da população (N) | 100 | [globals.h](src/moead/headers/globals.h) `SIZE_OF_POPULATION` |
| Tamanho da vizinhança (T) | 10 | [moead.cpp:37](src/moead/src/metaheuristics/moead/moead.cpp#L37) |
| Probabilidade de crossover | 0.5 | [moead.cpp:35](src/moead/src/metaheuristics/moead/moead.cpp#L35) |
| Probabilidade de mutação | 0.5 (um sorteio para os dois filhos) | [moead.cpp:36](src/moead/src/metaheuristics/moead/moead.cpp#L36) |
| Operador de crossover | `onepoint` (padrão), `blockswap` ou `blockcopy` | 5º argumento, `crossover_type` em [main.cpp](src/moead/main.cpp#L15) |
| Block Swap: tamanho do retângulo | 10% a 50% da largura e da altura da zona | `BLOCK_SWAP_MIN_FRACTION`, `BLOCK_SWAP_MAX_FRACTION` em [block_swap_crossover.cpp](src/moead/src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp) |
| Block Swap: tentativas por zona | 10 | `BLOCK_SWAP_MAX_ATTEMPTS` |
| BlockCopy: tamanho do retângulo | 10% a 50% da largura e da altura da zona | `BLOCK_COPY_MIN_FRACTION`, `BLOCK_COPY_MAX_FRACTION` em [block_copy_crossover.cpp](src/moead/src/global_modules/genetic_operators/crossover/block_copy_crossover.cpp) |
| Operador de mutação | `random` (padrão) ou `worst` | 6º argumento, `mutation_type` em [main.cpp](src/moead/main.cpp#L16) |
| Worst Turbine Mutation: turbinas movidas e foco | 1 turbina; foco em custo com probabilidade λᵢ₁ do subproblema | [worst_turbine_mutation.cpp](src/moead/src/global_modules/genetic_operators/mutation/worst_turbine_mutation.cpp) |
| Vetores de peso | Uniformes, `λᵢ = (i/(N−1), 1 − i/(N−1))` | `build_weight_vector` |
| Função de agregação | Tchebycheff, sem normalização | `calculate_gte` |
| Ponto de referência z* | Máximo de cada objetivo no arquivo Pareto, atualizado a cada subproblema | Step 2.3 |
| Seleção de pais | Apenas da vizinhança B(i) (sem cruzamento global) | Step 2.1 |
| Substituição | Todos os vizinhos em que o filho é melhor ou igual (sem limite `nr`) | Step 2.4 |
| Orçamento | 1.000.000 avaliações | [main.cpp:11](src/moead/main.cpp#L11) |
| Tamanho máximo do arquivo Pareto | 300 | [BoundedParetoSet.h](src/moead/headers/adaptative_grid/BoundedParetoSet.h) `MAXARCSIZE` |
| Profundidade do grid adaptativo | 5 bits por objetivo (32×32 células) | [param.h](src/moead/headers/adaptative_grid/param.h) |
| Vento / ângulo | 10 m/s / 30° (padrão, usado pelo `moead.sh`) | [instance_info.cpp](src/moead/instance_info.cpp#L20) |
| Execuções por instância | 20 | `scripts/main.sh` |
| Semente | `time(nullptr)` / `random_device` (não reprodutível) | `instance_info.cpp`, `moead.cpp` |

---

## 7. Problemas conhecidos

Pontos encontrados na leitura do código, do mais grave para o menos grave.

### 7.1 Status

| # | Problema | Categoria | Status |
|---|---|---|---|
| 1 | Vento fixo em 10 m/s e 30° | Invalidam a comparação | ❌ Aberto |
| 2 | Tchebycheff sem normalização | Qualidade dos resultados | ❌ Aberto |
| 3 | Sementes não reprodutíveis | Qualidade dos resultados | ❌ Aberto |
| 4 | Laço infinito da mutação em zona sem fundação livre | Menor | ❌ Aberto |
| 5 | Interpolação da curva da turbina fixa usa `power` | Menor | ❌ Aberto |
| 6 | Orçamento ultrapassado em até ~400 avaliações | Menor | ❌ Aberto |
| 7 | Código morto (`tch_vector`, `updateEP`, `nsga2.cpp`) | Menor | ❌ Aberto |
| 8 | `cout << argc` no topo do `log.txt` | Menor | ❌ Aberto |
| 9 | `__pycache__/*.pyc` versionado | Menor | ❌ Aberto |
| C1 | `makeall.sh` falhava (binário `nsga2`) e os `.sh` estavam sem permissão de execução | Impediam os experimentos | ✅ Corrigido (`02a50f2fc`) |
| C2 | `scripts/logs/` não existia e nenhum lote rodava | Impediam os experimentos | ✅ Corrigido (`02a50f2fc`) |
| C3 | Rodar o lote sobrescrevia os resultados de referência em `dataset/moead/` | Qualidade dos resultados | ✅ Corrigido (`de4bb2b85`) |
| C4 | Typo na dominância (`s2.fitness.second > s2.fitness.second`) | Qualidade dos resultados | ✅ Corrigido (`4a774c913`) |
| C5 | Ids de turbina repetidos entre zonas superestimavam a potência | Erravam os objetivos | ✅ Corrigido (`ae245d2c0`) |

### 7.2 Problemas abertos

**Invalidam a comparação dos experimentos**

1. **Vento fixo em 10 m/s e 30°.** O `moead.sh` passa só instância e pasta, e o sorteio pela rosa dos ventos está comentado em `instance_info.cpp`. Todas as instâncias e execuções usam o mesmo cenário, diferente dos resultados de referência em `dataset/reference/moead/`, que têm ventos sorteados.

**Afetam a qualidade ou a validade dos resultados**

2. **Tchebycheff sem normalização**: o custo (~10⁷) e a potência (~10⁶) estão em escalas diferentes. O termo do custo domina o `max` na maior parte dos pesos, e os subproblemas não se espalham uniformemente pela frente.
3. **Sementes não reprodutíveis** (`srand(time(nullptr))`, `random_device`): não dá para repetir uma execução específica.

**Menores**

4. A mutação entra em laço infinito se uma zona não tiver posição livre (todas as fundações ocupadas).
5. Na leitura da curva da turbina fixa, a interpolação usa `power` (da turbina móvel, ainda 0) em vez do valor anterior. Só afeta ventos que não estão exatamente na tabela.
6. O orçamento pode ser ultrapassado em até ~400 avaliações (critério checado só no início da geração).
7. Código morto: `tch_vector` e `updateEP` são calculados ou declarados, mas não usados; `nsga2.cpp` é compilado, mas não chamado.
8. `instance_info.cpp` imprime `argc` (`cout << argc`) no início de todo `log.txt`.
9. `src/instance_generation/__pycache__/*.pyc` está versionado.

### 7.3 Problemas corrigidos

- **C1. Build** (`02a50f2fc`): o `Makefile` gerava `nsga2`, mas o `makeall.sh` movia `moead`, e o `mv` falhava. Os `.sh` estavam versionados com modo 100644, então o `nohup ./main.sh` falhava com *Permission denied*. Agora `OUTPUT = moead` e os scripts têm bit de execução.
- **C2. Pasta de logs** (`02a50f2fc`): o `batch.sh` redirecionava para `./logs/batch_*.txt` sem a pasta existir, e nenhum lote rodava. Agora ele faz `mkdir -p ./logs`.
- **C3. Sobrescrita dos resultados de referência** (`de4bb2b85`): os resultados de A, 172, 240 e 456 ficavam em `dataset/moead/`, a mesma pasta em que o lote grava, e o `create_folders.sh` não criava nada porque a pasta já existia. Agora eles ficam em `dataset/reference/moead/`.
- **C4. Typo na dominância** (`4a774c913`): `dominatesP` comparava `s2.fitness.second > s2.fitness.second`, sempre falso, e o arquivo Pareto guardava soluções dominadas de mesmo custo (284 nos arquivos finais de referência). Agora compara `s1.fitness.second > s2.fitness.second`, como a função `dominates` de `dominates.cpp`.
- **C5. Ids de turbina repetidos entre zonas** (`ae245d2c0`): o crossover reiniciava o contador de ids em cada zona (`int id_new = fixd.size();` dentro de `crossoverAux`), e turbinas de zonas diferentes recebiam os mesmos ids. Como o `calculate_interference` devolve 0 para ids iguais, o wake entre elas era ignorado e a potência saía superestimada em toda instância com mais de uma zona. O mesmo layout chegava a ter potências diferentes (no teste com a instância A, 386.025,6 no pai e 386.071,5 no filho idêntico). Nos resultados de referência, os empates de custo aparecem quase só nas instâncias de 2 zonas (A: 15, 172: 189; 240 e 456, de 1 zona: 1 cada). Agora `id_new` é declarado em `crossover` e passado por referência a `crossoverAux`, com um contador único para todas as zonas. No mesmo teste, 0 de 200 filhos têm ids repetidos e as cópias do layout do pai têm a mesma potência. **Os resultados de referência em `dataset/reference/moead/` foram gerados antes desta correção.**

### 7.4 Desvios do MOEA/D original

Não são bugs, mas valem registro:

- B(i) exclui o próprio i; no artigo original, i faz parte da própria vizinhança.
- Não há limite `nr` de substituições: fiel ao MOEA/D de 2007, mas reduz a diversidade (um bom filho pode virar 10 cópias).

---

## 8. Referência de comandos do supercomputador

**Acesso**

```bash
ssh -p <porta> <usuario>@<host>
```

**Copiar arquivos**

```bash
scp -r -P <porta> "C:/Users/metafeatures_extraction/*" <usuario>@<host>:/home/gpu/metafeatures_extraction/
scp -r -P <porta> <usuario>@<host>:/home/gpu/meta_learning/scripts "C:/Users/last_run"
```

**Gerenciar arquivos**

```bash
cp -r meta_learning_v3/instances metafeatures_extraction/
mv source_directory target_directory
```

**Rodar em background e encerrar**

```bash
nohup ./script.sh &> output.log 2>&1 &
pkill -f "./main.sh"
```

---

## Referências

- Q. Zhang, H. Li. *MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition*. IEEE TEC, 2007.
- X. Ju, F. Liu. *Wind farm layout optimization using self-informed genetic algorithm with information guided exploitation*. Applied Energy 248, 2019. https://doi.org/10.1016/j.apenergy.2019.04.084
- R. Tran, J. Wu, C. Denison, T. Ackling, M. Wagner, F. Neumann. *Fast and effective multi-objective optimisation of wind turbine placement*. GECCO 2013. https://doi.org/10.1145/2463372.2463541
- M. Mayo, C. Zheng. *BlockCopy-based operators for evolving efficient wind farm layouts*. IEEE, 2015.
- W. Li, E. Özcan, R. John. *Multi-objective evolutionary algorithms and hyper-heuristics for wind farm layout optimisation*. Renewable Energy 105, 2017. https://doi.org/10.1016/j.renene.2016.12.022
- Y. Chen, H. Li, B. He, P. Wang, K. Jin. *Multi-objective genetic algorithm based innovative wind farm layout optimization method*.
- S. Rodrigues, P. Bauer, P. A. N. Bosman. *Multi-objective optimization of wind farm layouts – Complexity, constraint handling and scalability*.
