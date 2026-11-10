# CLAUDE.md

Guia para trabalhar neste repositório. Detalhes completos do algoritmo, da modelagem e dos hiperparâmetros estão no [README.md](README.md).

## O que é

MOEA/D (e NSGA-II) em C++17 para o WFLOP biobjetivo: posicionar turbinas eólicas offshore em fundações discretas, **maximizando** `fitness.first = −custo` e `fitness.second = potência`. Todo o código fica em `src/comolsd/` (nome histórico) e gera um único executável `main`.

## Comandos

```bash
./makeall.sh                                   # make rebuild em src/comolsd → scripts/exe/main
make -C src/comolsd rebuild                    # só compila (gera src/comolsd/main)

cd scripts
./exe/main <inst> <saida/> <angulo> <vento> <moead_hybrid|nsga2_hybrid> <pls|apls|none>
mkdir -p logs/moead_hybrid && ALGO=moead_hybrid LS=apls ./batch.sh
```

Não há testes automatizados. Para validar uma mudança, compile e rode uma instância pequena (ex.: `A`) com orçamento reduzido.

## Entrada

O executável lê `../instances/` relativo ao **diretório de trabalho** (rodando de `scripts/`, é `<repo>/instances/`, que não é versionado). As instâncias estão em https://github.com/gustavojorge/MO_WFLOP-archive/tree/main/instances. Lá a pasta se chama `sites/`, mas o código lê `site/`; use `ln -s sites instances/site`.

- `instances/site/<inst>/availablePositions.txt` — `x y prof custo zona`: fundações candidatas num grid de ~160 m dentro dos polígonos das zonas (já sem buracos/estruturas), zona 1..3, ordenado por zona
- `instances/site/<inst>/fixed_wf.txt` — turbinas fixas (parque existente, "zona" seguinte às disponíveis)
- `geometry.txt`, `info.json`, `plot.png` e `instances/wind/RVO_*.txt` existem, mas não são lidos pelo código
- `instances/site/<inst>/turbines_per_zone.txt` — só para instâncias numéricas (A–J estão em `turbines_instace` em `instance_info.cpp`)
- `instances/wtg/NREL-15-240.txt` (móvel) e `NREL-10-179.txt` (fixa) — `vento potência Ct`

Os arquivos precisam terminar com quebra de linha: o loader faz `pop_back()` da última leitura duplicada.

## Arquitetura

- `main.cpp` — globais (`pareto`, `countRevalue`, `stop_criteria = 1000000`, `algorithm`, `root_folder`) e o dispatch entre `moead_hybrid` e `nsga2_hybrid`. `moead()` (sem busca local) existe, mas não é chamado.
- `instance_info.cpp` — globais do problema: `foundations[z]`, `fixd`, `turbines_per_zone`, `num_zones`, `wind`, `angle`, `power`, `thrust_coef`.
- `src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp` — solução aleatória, `calculate_cost`, `calculate_power` e `calculate_interference` (wake de Jensen com soma quadrática).
- `src/metaheuristics/moead/` — `moead_hybrid.cpp` (laço principal) e `modules/` (pesos, vizinhança, z*, Tchebycheff).
- `src/global_modules/genetic_operators/` — `crossover.cpp` e `mutation.cpp` (os únicos compilados).
- `src/global_modules/local_search/` — PLS, Anytime-PLS e `getNeighborhood` (vizinhos = chamadas de `mutation`).
- `src/adaptive_grid/` — `BoundedParetoSet` (arquivo global, máx. 300, grid adaptativo) e `ParetoSetLS` (arquivo local da busca local).
- Headers ficam em `headers/` espelhando `src/`. Um `.cpp` novo precisa entrar em `SRC_FILES` no `Makefile`.

## Invariantes que precisam ser mantidas

- **Contagem de avaliações**: quem avalia uma solução (`calculate_cost` + `calculate_power`) também deve chamar `pareto->adicionarSol`, incrementar `countRevalue` e salvar o snapshot quando `countRevalue % 100000 == 0`. Hoje isso está duplicado em `population.cpp`, `crossover.cpp` e `mutation.cpp`.
- **Representação**: `grid[z][p]` guarda o id da turbina ou `0` (vazio); `turbines[z][k].index` aponta para a fundação `p`. Os dois precisam ficar consistentes. Ids das móveis começam em `fixd.size()`.
- **Viabilidade por construção**: número exato de turbinas por zona, uma turbina por fundação, sem troca de zona. Novos operadores precisam preservar isso (não há penalização).
- **Sentido dos objetivos**: tudo é maximização; z* é o máximo e as saídas escrevem `abs(custo)`.
- **Memória**: soluções são `new`/`delete` manuais. `ParetoSet::adicionarSol` copia a solução recebida e o arquivo é dono das cópias; `ParetoSetLS::clear()` não libera memória (quem chama libera).

## Armadilhas conhecidas

- `LS=none` no `moead_hybrid` dá `delete` em ponteiros do arquivo Pareto global (use-after-free).
- `dominatesP` em `ParetoSet.cpp` tem typo (`s2.fitness.second > s2.fitness.second`).
- Tchebycheff sem normalização (custo ~10⁷ vs potência ~10⁶).
- Scripts usam `"$ALGO_$LS"` (vira `$ALGO_` vazio); a saída vai para `dataset/<LS>/`.
- `run.sh` lê vento/ângulo de `dataset/moead/<inst>/<run>/log.txt`; sem esse arquivo, usa 0.0 (potência zero).
- Sementes não reprodutíveis (`srand(time(nullptr))`, `random_device`).

A lista completa está na seção "Problemas conhecidos" do README.

## Estilo de código

- Não adicionar comentários ou docstrings em código que não os tem.
- Seguir o estilo existente: funções livres, globais `extern` declaradas nos headers, `using namespace std`.
