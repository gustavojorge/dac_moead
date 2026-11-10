# dac_moead — MOEA/D para o Wind Farm Layout Optimization Problem (WFLOP) multiobjetivo

Este repositório contém a implementação em C++ do **MOEA/D** (Multi-Objective Evolutionary Algorithm based on Decomposition), além do NSGA-II e de buscas locais Pareto (PLS e Anytime-PLS), aplicados ao problema de **posicionamento de turbinas em parques eólicos offshore** com dois objetivos: **minimizar o custo de fundação** e **maximizar a potência gerada**.

Também contém os scripts shell para compilar e executar os experimentos em lote (várias instâncias × 20 execuções).

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
├── makeall.sh                     # Compila src/comolsd e move o executável para scripts/exe/main
├── scripts/
│   ├── batch.sh                   # Dispara lotes de instâncias em background (nohup)
│   ├── main.sh                    # Para cada instância do lote, roda 20 execuções
│   ├── exe/                       # Executáveis compilados (main)
│   └── smaller_scripts/
│       ├── create_folders.sh      # Cria dataset/<...>/<instância>/<run>/
│       └── run.sh                 # Chama ./exe/main para (instância, run)
├── dataset/
│   └── moead/<inst>/<run>/        # Resultados de referência do MOEA/D (log.txt guarda vento/ângulo)
└── src/comolsd/
    ├── main.cpp                   # Ponto de entrada: escolhe algoritmo e busca local
    ├── instance_info.cpp          # Leitura da instância, curvas das turbinas, turbinas fixas
    ├── Makefile
    ├── headers/                   # .h espelhando a árvore de src/
    └── src/
        ├── metaheuristics/
        │   ├── moead/             # moead.cpp, moead_hybrid.cpp e modules/
        │   ├── nsga2/
        │   └── comolsd/
        ├── global_modules/
        │   ├── generate_initial_population/   # Solução aleatória + funções objetivo + wake effect
        │   ├── genetic_operators/             # crossover/ e mutation/
        │   ├── local_search/                  # PLS, Anytime-PLS, vizinhança, factory
        │   ├── dominates.cpp, isEqual.cpp
        ├── adaptive_grid/         # Arquivo Pareto limitado com grid adaptativo
        └── instance_generation/   # Gerador de instâncias sintéticas (Python)
```

Todo o código está em `src/comolsd` (o nome é histórico): um único executável `main` contém MOEA/D, NSGA-II e as buscas locais.

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

Implementação principal: [moead_hybrid.cpp](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp) (versão chamada pelo `main`). A versão [moead.cpp](src/comolsd/src/metaheuristics/moead/moead.cpp) é idêntica sem a busca local, mas **não é chamada** pelo `main.cpp` atual.

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
   Com probabilidade ls_prob = avaliações / 1.000.000: executa a busca local sobre o arquivo Pareto
```

### 2.3 Módulos principais

| Módulo | Arquivo | Função |
|---|---|---|
| Vetores de peso | [generate_weight_vectors.cpp](src/comolsd/src/metaheuristics/moead/modules/generate_weight_vectors.cpp) | Gera N pares (w1, w2) igualmente espaçados com w1 + w2 = 1 |
| Vizinhança | [generate_neighborhood.cpp](src/comolsd/src/metaheuristics/moead/modules/generate_neighborhood.cpp) | Para cada λᵢ, ordena os demais por distância euclidiana e pega os T primeiros |
| Ponto de referência | [get_best_z_point.cpp](src/comolsd/src/metaheuristics/moead/modules/get_best_z_point.cpp) | z* = máximo de cada objetivo na população |
| Tchebycheff | [tchebycheff.cpp](src/comolsd/src/metaheuristics/moead/modules/tchebycheff.cpp) | `calculate_gte`: max λk·\|fk − z*k\| (sem normalização dos objetivos) |
| População externa (EP) | [updateEP.cpp](src/comolsd/src/metaheuristics/moead/modules/updateEP.cpp) | Implementação clássica da EP — **não usada**; o papel de EP é do arquivo Pareto global |
| Arquivo Pareto global | [BoundedParetoSet.cpp](src/comolsd/src/adaptive_grid/BoundedParetoSet.cpp), [ParetoSet.cpp](src/comolsd/src/adaptive_grid/ParetoSet.cpp) | Guarda até 300 soluções não dominadas; quando cheio, remove uma solução da célula mais populosa de um grid adaptativo (32×32 células) |
| Busca local | [local_search_factory.cpp](src/comolsd/src/global_modules/local_search/local_search_factory.cpp), [pareto_ls.cpp](src/comolsd/src/global_modules/local_search/pareto_ls.cpp), [anytime_pls.cpp](src/comolsd/src/global_modules/local_search/anytime_pls.cpp) | `pls`, `apls` ou `none` |

### 2.4 Arquivo Pareto global (`pareto`)

Variável global `BoundedParetoSet * pareto` declarada em [main.cpp](src/comolsd/main.cpp). **Toda avaliação** (geração aleatória, crossover, mutação, vizinhos da busca local) chama `pareto->adicionarSol(...)`. É esse arquivo que é salvo como resultado do algoritmo — não a população do MOEA/D.

### 2.5 Busca local (versão híbrida)

- **PLS** (`pls`): copia o arquivo Pareto para um arquivo local; enquanto houver solução não explorada, escolhe uma aleatoriamente, gera 163 vizinhos (mutações distintas) e tenta inserir todos.
- **Anytime-PLS** (`apls`): fase 1 escolhe a solução não explorada de maior contribuição de hipervolume (`ohiv`) e para no primeiro vizinho que a domine (ou no primeiro inserido); fase 2 explora a vizinhança completa.
- **none**: não faz busca local (ver [problemas conhecidos](#7-problemas-conhecidos)).

Os vizinhos são gerados com a mesma `mutation()` do algoritmo genético, então cada vizinho conta como avaliação e é inserido no arquivo Pareto global. A população do MOEA/D **não** é alterada pela busca local.

---

## 3. Como executar

### 3.1 Compilação

```bash
./makeall.sh
```

Roda `make rebuild` em `src/comolsd` (g++ com `-std=c++17 -O2`), move o executável para `scripts/exe/main` e dá permissão de execução aos scripts.

### 3.2 Estrutura de instâncias esperada

As instâncias ficam no repositório [MO_WFLOP-archive/instances](https://github.com/gustavojorge/MO_WFLOP-archive/tree/main/instances) e não são versionadas aqui.

O executável lê as instâncias de **`../instances`, relativo ao diretório de trabalho** (`pathFolders = ".."` em [instance_info.cpp](src/comolsd/instance_info.cpp)). Como os scripts executam `./exe/main` a partir de `scripts/`, a pasta deve ficar na **raiz do repositório**.

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

As instâncias sintéticas são geradas por [script.py](src/comolsd/src/instance_generation/script.py) (detalhes na seção 4.1). O script cria `../instances/<N>/`; mova para `instances/site/<N>/`.

### 3.3 Linha de comando do executável

```bash
./exe/main <instância> <pasta_saida/> <ângulo> <vento> <algoritmo> <busca_local>
```

| Arg | Exemplo | Descrição |
|---|---|---|
| `instância` | `A`, `240` | Nome da pasta em `instances/site/` |
| `pasta_saida/` | `../dataset/x/A/1/` | Prefixo de saída — **precisa da barra final** |
| `ângulo` | `330` | Direção do vento em graus (só é lido se houver ≥ 4 argumentos; padrão 30) |
| `vento` | `13` | Velocidade do vento em m/s (padrão 10) |
| `algoritmo` | `moead_hybrid` | `moead_hybrid` ou `nsga2_hybrid` (padrão) |
| `busca_local` | `apls` | `pls`, `apls` (padrão) ou `none` |

### 3.4 Execução em lote

```bash
cd scripts
mkdir -p logs/moead_hybrid
ALGO=moead_hybrid LS=apls ./batch.sh
```

- `batch.sh` cria as pastas de saída e dispara 4 lotes em background com `nohup`, cada um chamando `main.sh` com uma lista de instâncias.
- `main.sh` roda **20 execuções** de cada instância, chamando `smaller_scripts/run.sh <instância> <run>`.
- `run.sh` lê o **ângulo e o vento** de `dataset/${ANGLE_WIND_SOURCE:-moead}/<inst>/<run>/log.txt`. Assim, todos os algoritmos usam o mesmo cenário de vento por (instância, run), sorteado originalmente pelo MOEA/D. Se o log não existir, ângulo e vento ficam **0.0** (potência zero).

### 3.5 Saídas

Em `<pasta_saida>/`:

| Arquivo | Conteúdo |
|---|---|
| `log.txt` | Parâmetros, nº de turbinas, vento/ângulo e tempo (`time`) |
| `infoRun.txt` | Por geração: avaliações acumuladas, tamanho do arquivo Pareto e se a busca local rodou |
| `<inst>_<algoritmo>_<N>.txt` | Arquivo Pareto a cada 100.000 avaliações (`N = 100000, …, 1000000`). Cada linha: `custo potência` (custo em valor absoluto) |

---

## 4. Modelagem do parque eólico

Código: [generate_rSolution.cpp](src/comolsd/src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp), [instance_info.cpp](src/comolsd/instance_info.cpp), [headers/wind_farm_entities/](src/comolsd/headers/wind_farm_entities/).

### 4.1 Instâncias: zonas e turbinas fixas

Uma instância descreve um sítio offshore com **até 3 zonas**. As coordenadas estão em metros, numa faixa de ~100–150 km. Há dois tipos de zona:

- **Zonas disponíveis** (1, 2, …): polígonos onde o algoritmo pode instalar turbinas novas.
- **Zona fixa** (o número seguinte, ex.: 3): um parque já existente, fora dos polígonos. Suas turbinas não se movem, mas entram no cálculo do wake e da potência.

![Instância 10](https://raw.githubusercontent.com/gustavojorge/MO_WFLOP-archive/main/instances/sites/10/plot.png)

*Instância 10: duas zonas disponíveis (contorno magenta, cor = profundidade em m) e o parque fixo (pontos laranja).*

**Como uma zona é construída** (gerador [script.py](src/comolsd/src/instance_generation/script.py) e [gen_layout.py](src/comolsd/src/instance_generation/gen_layout.py)):

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

As instâncias A–J têm um `info.json` com área, número de buracos e opções de quantidade de turbinas para densidades de potência de 4, 6 e 8 MW/km² (ex.: zona 1 de A: 13, 20 ou 26 turbinas de 15 MW). As quantidades usadas estão fixas em `turbines_instace` ([instance_info.cpp](src/comolsd/instance_info.cpp#L26)).

**Como o código usa esses dados**

- `availablePositions.txt` → `foundations[z]` (`Foundation {x, y, cost}`). A posição da fundação no arquivo, dentro da zona, é o índice `p` usado em `grid[z][p]`.
- `fixed_wf.txt` → `fixd` (turbinas NREL 10 MW).
- `geometry.txt`, `info.json` e `plot.png` **não** são lidos: os polígonos já estão representados pelas fundações que sobraram dentro deles.
- O grid de 160 m é menor que o diâmetro do rotor (240 m), então duas fundações vizinhas podem receber turbinas. Não há distância mínima, e a proximidade só é penalizada pelo wake.

### 4.2 Turbinas e vento (dados do wake effect)

- **Turbinas móveis**: as que o algoritmo posiciona. Modelo NREL 15 MW: diâmetro 240 m, altura do cubo 150 m, curva `wtg/NREL-15-240.txt`.
- **Turbinas fixas** (`fixd`): modelo NREL 10 MW: diâmetro 179 m, altura 119 m, curva `wtg/NREL-10-179.txt`.
- **Curvas** (`vento potência Ct`): para a velocidade do vento livre `U`, o código lê a potência (`P_curva`, em MW) e o coeficiente de empuxo `Ct` de cada modelo. Se `U` não estiver na tabela, usa a média entre as linhas vizinhas. Esses valores são os mesmos para todas as turbinas do mesmo modelo, pois são lidos com o vento livre, antes do wake.
- **Vento**: um único cenário por execução, com velocidade `wind` (m/s) e direção `angle` (graus; o vetor do vento é `(cos θ, sin θ)`). Os dois vêm da linha de comando. Os arquivos `wind/RVO_*.txt` (rosas dos ventos de sítios offshore holandeses, ex.: HKN = Hollandse Kust Noord) têm a distribuição `ângulo velocidade probabilidade`, mas a leitura deles está comentada em `instance_info.cpp`. Nos experimentos, cada (instância, run) reaproveita o cenário salvo em `dataset/moead/<inst>/<run>/log.txt` (ex.: A/run 1 → vento 13 m/s, ângulo 330°).

### 4.3 Solução

```cpp
struct Solution {
    vector<vector<int>> grid;          // grid[z][p] = id da turbina na fundação p da zona z, 0 = vazia
    vector<vector<Turbine>> turbines;  // turbines[z] = turbinas da zona z (com x, y, índice da fundação)
    pair<double,double> fitness;       // (−custo, potência)
};
```

Os ids das turbinas móveis começam em `fixd.size()`, para não colidir com os ids das fixas (0..F−1) nem com o 0 usado como "vazio".

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

Modelo do tipo Jensen/Park ([calculate_interference](src/comolsd/src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp#L46)). Para a turbina *a* gerando esteira sobre a turbina *b*:

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

- **Máximo**: `stop_criteria = 1.000.000` avaliações ([main.cpp](src/comolsd/main.cpp#L12)). O contador global é `countRevalue`.
- **O que conta como 1 avaliação**: cada chamada de `calculate_cost` + `calculate_power` sobre uma nova solução, isto é, cada solução da população inicial, cada filho de crossover, cada mutação e cada vizinho gerado na busca local.
- **Como é gasto no MOEA/D**:
  - 100 avaliações na população inicial;
  - por subproblema: 0, 2 ou 4 avaliações (crossover com prob. 0.5 gera 2 filhos; mutação com prob. 0.5 muta os 2) → em média **~200 avaliações por geração**, ou seja, cerca de **5.000 gerações** sem busca local;
  - busca local: até 163 avaliações por solução explorada (+ re-sorteios quando um vizinho repetido é gerado). No Anytime-PLS, a fase 1 devolve 162 avaliações ao contador por solução explorada (conta só 1).
- **Parada**: o laço do MOEA/D só verifica o critério no início de cada geração, então pode ultrapassar 1.000.000 em até ~400 avaliações. As buscas locais verificam a cada solução explorada.
- **Snapshots**: a cada 100.000 avaliações o arquivo Pareto é salvo em disco.

---

## 5. Operadores genéticos

Arquivos: [crossover.cpp](src/comolsd/src/global_modules/genetic_operators/crossover/crossover.cpp) e [mutation.cpp](src/comolsd/src/global_modules/genetic_operators/mutation/mutation.cpp). Ambos avaliam a solução gerada, inserem no arquivo Pareto global e incrementam o contador de avaliações.

**Crossover (um ponto por zona, com reparo)** — para cada zona:
1. Sorteia um ponto de corte `c` sobre o vetor de fundações.
2. Copia do pai A todas as turbinas das posições `0..c`.
3. Percorre o pai B **do fim para o começo**, adicionando as turbinas cujas posições ainda estão livres no filho, até atingir `turbines_per_zone[z]`.

O resultado sempre tem o número correto de turbinas por zona. O MOEA/D gera dois filhos: `crossover(A, B)` e `crossover(B, A)`.

**Mutação (realocação de uma turbina)**: sorteia **uma zona** e duas posições, uma ocupada e uma livre, e move a turbina para a posição livre. A cardinalidade é preservada.

`mutationElisama.cpp` e `mutationFeature.cpp` são variantes antigas (usadas na extração de metafeatures) e não são compiladas.

---

## 6. Hiperparâmetros utilizados

| Parâmetro | Valor | Onde |
|---|---|---|
| Nº de subproblemas / tamanho da população (N) | 100 | [globals.h](src/comolsd/headers/globals.h) `SIZE_OF_POPULATION` |
| Tamanho da vizinhança (T) | 10 | [moead_hybrid.cpp:39](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp#L39) |
| Probabilidade de crossover | 0.5 | [moead_hybrid.cpp:37](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp#L37) |
| Probabilidade de mutação | 0.5 (um sorteio para os dois filhos) | [moead_hybrid.cpp:38](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp#L38) |
| Vetores de peso | Uniformes, `λᵢ = (i/(N−1), 1 − i/(N−1))` | `build_weight_vector` |
| Função de agregação | Tchebycheff, sem normalização | `calculate_gte` |
| Ponto de referência z* | Máximo de cada objetivo no arquivo Pareto, atualizado a cada subproblema | Step 2.3 |
| Seleção de pais | Apenas da vizinhança B(i) (sem cruzamento global) | Step 2.1 |
| Substituição | Todos os vizinhos em que o filho é melhor ou igual (sem limite `nr`) | Step 2.4 |
| Orçamento | 1.000.000 avaliações | [main.cpp:12](src/comolsd/main.cpp#L12) |
| Tamanho máximo do arquivo Pareto | 300 | [BoundedParetoSet.h](src/comolsd/headers/adaptative_grid/BoundedParetoSet.h) `MAXARCSIZE` |
| Profundidade do grid adaptativo | 5 bits por objetivo (32×32 células) | [param.h](src/comolsd/headers/adaptative_grid/param.h) |
| Tamanho da vizinhança da busca local | 163 vizinhos | `pareto_ls.cpp`, `anytime_pls.cpp` |
| Probabilidade de busca local por geração | `avaliações / 1.000.000` (cresce de 0 a 1) | [moead_hybrid.cpp:160](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp#L160) |
| Execuções por instância | 20 | `scripts/main.sh` |
| Semente | `time(nullptr)` / `random_device` (não reprodutível) | `instance_info.cpp`, `moead_hybrid.cpp` |

---

## 7. Problemas conhecidos

Pontos encontrados na leitura do código, importantes para futuras alterações:

1. **`LS=none` corrompe o arquivo Pareto.** A busca local `none` devolve os próprios ponteiros do arquivo global, e o `moead_hybrid` dá `delete` neles ([moead_hybrid.cpp:150](src/comolsd/src/metaheuristics/moead/moead_hybrid.cpp#L150)). Para rodar o MOEA/D "puro", corrija o `none` para devolver cópias ou chame `moead()`.
2. **Dominância com typo**: em `dominatesP` ([ParetoSet.cpp:29](src/comolsd/src/adaptive_grid/ParetoSet.cpp#L29)) está `s2.fitness.second > s2.fitness.second` (sempre falso). Uma solução que melhora só a potência, com custo igual, não é considerada dominante.
3. **Tchebycheff sem normalização**: o custo (~10⁷) e a potência (~10⁶) estão em escalas diferentes, o que distorce a distribuição efetiva dos pesos.
4. **Caminhos nos scripts**: `"$ALGO_$LS"` em `create_folders.sh` e `run.sh` é lido pelo bash como a variável `$ALGO_` (vazia), então a saída vai para `dataset/<LS>/` e não para `dataset/<ALGO>_<LS>/`. O correto é `"${ALGO}_${LS}"`. A pasta `scripts/logs/<ALGO>/` também precisa existir antes do `batch.sh`.
5. **Vento zero**: se `dataset/moead/<inst>/<run>/log.txt` não existir, `run.sh` passa ângulo e vento 0. O dataset versionado só tem as instâncias A, 172, 240 e 456.
6. **Arquivo de layout nunca é gerado**: `printAllSolutionsLayout` só é chamado dentro da população inicial, quando o orçamento nunca está esgotado.
7. **Anytime-PLS**: `if(p->adicionarSol(...));{ break; }` ([anytime_pls.cpp:33](src/comolsd/src/global_modules/local_search/anytime_pls.cpp#L33)) tem um `;` extra e sempre sai do laço.
8. `ls_prob = 1 / stop_criteria` é divisão inteira (= 0) na primeira geração.
9. `tch_vector` e `updateEP` são calculados ou declarados, mas não usados.
10. A mutação entra em laço infinito se uma zona não tiver posição livre (todas as fundações ocupadas).
11. Na leitura da curva da turbina fixa, a interpolação usa `power` (da turbina móvel, ainda 0) em vez do valor anterior — só afeta ventos que não estão exatamente na tabela.

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
- Y. Chen, H. Li, B. He, P. Wang, K. Jin. *Multi-objective genetic algorithm based innovative wind farm layout optimization method*.
- S. Rodrigues, P. Bauer, P. A. N. Bosman. *Multi-objective optimization of wind farm layouts – Complexity, constraint handling and scalability*.
