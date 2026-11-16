# CLAUDE.md

Guia para trabalhar neste repositório. Detalhes completos do algoritmo, da modelagem e dos hiperparâmetros estão no [README.md](README.md).

## O que é

MOEA/D em C++17 para o WFLOP biobjetivo: posicionar turbinas eólicas offshore em fundações discretas, **maximizando** `fitness.first = −custo` e `fitness.second = potência`. O código fica em `src/moead/` (vindo da branch `moead-runner` do MO_WFLOP-experiment-runner) e não tem busca local nem versões híbridas.

## Comandos

```bash
make -C src/moead rebuild && mkdir -p scripts/exe && mv src/moead/nsga2 scripts/exe/moead
chmod +x makeall.sh scripts/*.sh scripts/smaller_scripts/*.sh
bash makeall.sh                                # hoje falha: o Makefile gera "nsga2", o script move "moead"

cd scripts
./exe/moead <inst> <saida/> [<angulo> <vento>] # padrão: 30° e 10 m/s
mkdir -p logs && ./batch.sh                    # 52 processos nohup, 20 runs por instância
```

Não há testes automatizados. Para validar uma mudança, compile e rode uma instância pequena (ex.: `A`) com orçamento reduzido.

## Entrada

O executável lê `../instances/` relativo ao **diretório de trabalho** (rodando de `scripts/`, é `<repo>/instances/`, que não é versionado). As instâncias estão em https://github.com/gustavojorge/MO_WFLOP-archive/tree/main/instances. Lá a pasta se chama `sites/`, mas o código lê `site/`; use `ln -s sites instances/site`.

- `instances/site/<inst>/availablePositions.txt` — `x y prof custo zona`: fundações candidatas num grid de ~160 m dentro dos polígonos das zonas (já sem buracos/estruturas), zona 1..3, ordenado por zona
- `instances/site/<inst>/fixed_wf.txt` — turbinas fixas (parque existente, "zona" seguinte às disponíveis)
- `instances/site/<inst>/turbines_per_zone.txt` — só para instâncias numéricas (A–J estão em `turbines_instace` em `instance_info.cpp`)
- `instances/wtg/NREL-15-240.txt` (móvel) e `NREL-10-179.txt` (fixa) — `vento potência Ct`
- `geometry.txt`, `info.json`, `plot.png` e `instances/wind/RVO_*.txt` existem, mas não são lidos pelo código

Os arquivos precisam terminar com quebra de linha: o loader faz `pop_back()` da última leitura duplicada.

## Arquitetura

- `main.cpp` — globais (`pareto`, `countRevalue`, `stop_criteria = 1000000`, `algorithm = "moead"`, `instance`, `root_folder`); lê `argv[1]` e `argv[2]` e chama `moead()`.
- `instance_info.cpp` — globais do problema: `foundations[z]`, `fixd`, `turbines_per_zone`, `num_zones`, `wind`, `angle`, `power`, `thrust_coef`. Lê `argv[3]`/`argv[4]` (ângulo/vento) só se `argc >= 5`.
- `src/global_modules/generate_initial_population/generate_random_solution/generate_rSolution.cpp` — solução aleatória, `calculate_cost`, `calculate_power` e `calculate_interference` (wake de Jensen com soma quadrática).
- `src/metaheuristics/moead/` — `moead.cpp` (laço principal) e `modules/` (pesos, vizinhança, z*, Tchebycheff).
- `src/global_modules/genetic_operators/` — `crossover.cpp` e `mutation.cpp` (os únicos compilados).
- `src/adaptive_grid/` — `BoundedParetoSet` (arquivo global, máx. 300, grid adaptativo).
- `src/metaheuristics/nsga2/` é compilado, mas não é chamado; `metafeatures/` e `src/utils/` não são compilados.
- Headers ficam em `headers/` espelhando `src/`. Um `.cpp` novo precisa entrar em `SRC_FILES` no `Makefile`.
- Scripts: `batch.sh` → `main.sh` (20 runs) → `smaller_scripts/moead.sh`, que grava em `dataset/moead/<inst>/<run>/`.

## Invariantes que precisam ser mantidas

- **Contagem de avaliações**: quem avalia uma solução (`calculate_cost` + `calculate_power`) também deve chamar `pareto->adicionarSol`, incrementar `countRevalue`, salvar o snapshot quando `countRevalue % 100000 == 0` e o layout quando `countRevalue >= stop_criteria`. Hoje isso está duplicado em `population.cpp`, `crossover.cpp` e `mutation.cpp`.
- **Representação**: `grid[z][p]` guarda o id da turbina ou `0` (vazio); `turbines[z][k].index` aponta para a fundação `p`. Os dois precisam ficar consistentes. Ids das móveis começam em `fixd.size()`.
- **Viabilidade por construção**: número exato de turbinas por zona, uma turbina por fundação, sem troca de zona. Novos operadores precisam preservar isso (não há penalização).
- **Sentido dos objetivos**: tudo é maximização; z* é o máximo e as saídas escrevem `abs(custo)`.
- **Memória**: soluções são `new`/`delete` manuais. `ParetoSet::adicionarSol` copia a solução recebida e o arquivo é dono das cópias. `population` é `vector<Solution>` (valores): substituir um vizinho copia o filho.

## Armadilhas conhecidas

- `makeall.sh` quebra: `OUTPUT = nsga2` no `Makefile`, mas o script move `moead`. Os `.sh` estão versionados sem bit de execução (100644).
- Vento fixo em 10 m/s e 30°: `moead.sh` passa só 2 argumentos e o sorteio da rosa dos ventos está comentado. Os resultados versionados em `dataset/moead/` (A, 172, 240, 456) usaram vento sorteado e são sobrescritos ao rodar o lote.
- `scripts/logs/` precisa existir antes do `batch.sh`.
- Tchebycheff sem normalização (custo ~10⁷ vs potência ~10⁶).
- `dominatesP` em `ParetoSet.cpp` tem typo (`s2.fitness.second > s2.fitness.second`).
- Sementes não reprodutíveis (`srand(time(nullptr))`, `random_device`).

A lista completa está na seção "Problemas conhecidos" do README.

## Estilo de código

- Não adicionar comentários ou docstrings em código que não os tem.
- Seguir o estilo existente: funções livres, globais `extern` declaradas nos headers, `using namespace std`.
