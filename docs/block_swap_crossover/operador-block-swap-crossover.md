# Operador Block Swap Crossover no MoWFLOP

Este documento descreve o operador **Block Swap Crossover** como foi implementado neste repositório: de onde ele vem, o que mudou na adaptação, as decisões tomadas e como ele foi validado.

Implementação: [block_swap_crossover.cpp](../../src/moead/src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp).

---

## 1. Origem na literatura

### 1.1 Tran et al. (GECCO 2013)

O operador foi proposto por Tran, Wu, Denison, Ackling, Wagner e Neumann em *Fast and effective multi-objective optimisation of wind turbine placement* (GECCO 2013, [doi:10.1145/2463372.2463541](https://doi.org/10.1145/2463372.2463541)), na seção 3.2 (Algoritmo 3).

O cenário do artigo é diferente do nosso:

| | Tran et al. (2013) | Este repositório |
|---|---|---|
| Posições | Coordenadas contínuas num terreno retangular | Fundações discretas, num grid de ~160 m dentro de polígonos |
| Zonas | Uma área | Até 3 zonas, cada uma com número fixo de turbinas |
| Restrições | Distância mínima entre turbinas | Só as da representação (sem distância mínima) |
| Objetivos | Energia, área do fecho convexo e comprimento de cabo | Custo das fundações e potência |
| Vento | Rosa dos ventos (distribuição de Weibull) | Um único cenário (velocidade e direção) |

O operador do artigo, para gerar o filho 1 a partir de P1 (base) e P2 (doador):

1. Sorteia um retângulo A₁ dentro do terreno.
2. Separa as turbinas de P1 em I₁ (dentro de A₁) e O₁ (fora).
3. Ajusta as bordas de A₁ para manter a distância de segurança.
4. Aplica o mesmo retângulo em P2, obtendo I₂ e O₂.
5. Se |I₂| ≥ |I₁|, move |I₂| − |I₁| turbinas de O₁ para dentro do retângulo e coloca as turbinas do retângulo nas posições de I₂. Caso contrário, sorteia outro retângulo, até T tentativas. Se nenhum servir, devolve o pai.

O filho 2 é gerado da mesma forma, trocando os papéis dos pais. No artigo, o crossover é aplicado com probabilidade 0,3.

### 1.2 Li, Özcan e John (2017)

*Multi-objective evolutionary algorithms and hyper-heuristics for wind farm layout optimisation* (Renewable Energy 105, 2017, [doi:10.1016/j.renene.2016.12.022](https://doi.org/10.1016/j.renene.2016.12.022)) reutiliza os operadores de Tran et al. (Block Swap e Movement Mutation) com probabilidade de crossover 0,3. O foco do artigo são as hiper-heurísticas de seleção, e não o operador.

### 1.3 O que a literatura mostra e o que não mostra

- O Block Swap é um operador **específico do domínio**, usado em pelo menos dois trabalhos de WFLOP multiobjetivo.
- **Nenhum dos dois artigos compara o Block Swap com outro crossover.** Não há evidência isolada de que ele seja superior a um corte simples.
- No artigo original, boa parte da motivação é a **viabilidade**: com coordenadas contínuas e distância mínima, um crossover genérico gera filhos inviáveis. No nosso problema, a representação já garante a viabilidade.

---

## 2. O crossover atual (`onepoint`) já é espacial

O crossover padrão do repositório ([crossover.cpp](../../src/moead/src/global_modules/genetic_operators/crossover/crossover.cpp)), para cada zona:

1. sorteia um corte `c` sobre os índices `p` das fundações;
2. copia do pai A as turbinas das fundações `0..c`;
3. completa com as turbinas do pai B, percorrendo do fim para o começo.

O índice `p` é a posição da fundação no `availablePositions.txt`, e esse arquivo lista as fundações **em varredura do grid**:

- na instância A, coluna por coluna, com x crescente (verificado em todas as 3.194 linhas consecutivas);
- nas instâncias 240 e 10, linha por linha.

Por isso fundações com índices próximos são vizinhas no mapa, e o intervalo `0..c` é uma **faixa contígua** do parque. O filho do `onepoint` é, grosso modo, um lado de A colado ao outro lado de B.

```
 onepoint (faixa)            blockswap (retângulo)
 A A A A | B B B B           A A A A A A A A
 A A A A | B B B B           A A [B B B] A A
 A A A A | B B B B           A A [B B B] A A
 A A A A | B B B B           A A A A A A A A
```

O argumento comum, de que o crossover de um ponto "ignora a geometria", **não vale aqui**. O Block Swap é uma **generalização em duas dimensões** de um corte que já é espacial: permite trocar uma região no meio da zona, e não só uma faixa que encosta na borda. Por isso o ganho esperado é modesto, e o operador deve ser avaliado empiricamente.

---

## 3. Algoritmo implementado

A interface é a mesma do crossover atual: `Solution block_swap_crossover(Solution &base, Solution &doador)`. O MOEA/D continua gerando dois filhos, `crossover(A, B)` e `crossover(B, A)`, e cada chamada sorteia seu próprio retângulo.

Para cada zona `z`, com τ_z turbinas:

1. **Retângulo envolvente da zona:** calcula os mínimos e máximos de x e y das fundações da zona (largura W e altura H).
2. **Sorteio de R:** escolhe ao acaso uma turbina do doador nessa zona e centraliza R nela. A largura de R é `W · U(0,1; 0,5)` e a altura é `H · U(0,1; 0,5)`, sorteadas de forma independente.
3. **Classificação:**
   - `O_A`: turbinas da base fora de R;
   - `I_A`: turbinas da base dentro de R;
   - `I_B`: turbinas do doador dentro de R.
4. **Aceitação:** R só é aceito se `|I_B| ≥ |I_A|`. Caso contrário, volta ao passo 2, até 10 tentativas.
5. **Montagem:** o filho recebe `O_A` mais `I_B`. As fundações dentro e fora de R são disjuntas, então nenhuma fundação recebe duas turbinas.
6. **Reparo da cardinalidade:** o total é `|O_A| + |I_B| = τ_z − |I_A| + |I_B| ≥ τ_z`. As `|I_B| − |I_A|` turbinas excedentes são retiradas ao acaso de `O_A`, e o bloco vindo do doador fica intacto.
7. **Falha:** se as 10 tentativas falharem, ou se o doador não tiver turbinas na zona, a zona do filho fica igual à da base.

Depois de todas as zonas:

8. **Filho igual à base:** se o conjunto de fundações ocupadas for igual ao da base, a base é devolvida **sem nova avaliação**.
9. **Avaliação:** caso contrário, calcula custo e potência, insere no arquivo Pareto, incrementa `countRevalue` e salva os snapshots, como o crossover atual.

### Exemplo (τ_z = 6)

```
 base (A)              doador (B)            filho
 A . A . . .           . . . B . .           A . . . . .
 . . ┌─────┐ A         . . ┌─────┐ .         . . ┌─────┐ .
 A . │ A . │ .         . B │ B B │ .         A . │ B B │ .
 . . │ . . │ .         . . │ . B │ B         . . │ . B │ .
 . A └─────┘ .         B . └─────┘ .         . A └─────┘ .
```

Neste exemplo:
- A tem 1 turbina dentro de R e 5 fora; B tem 3 dentro. Como 3 ≥ 1, o retângulo é aceito.
- Juntando as 5 de fora com as 3 de dentro, o filho teria 8 turbinas. As 2 excedentes saem de `O_A`, sorteadas.
- O filho termina com 3 turbinas de A fora de R e 3 de B dentro, total 6.

---

## 4. Decisões de projeto

| Decisão | Escolha | Motivo |
|---|---|---|
| Regra de aceitação e reparo | A de Tran et al.: `\|I_B\| ≥ \|I_A\|` e excesso retirado de fora de R | Nenhuma posição aleatória nova entra no filho, e o bloco do doador chega inteiro. A alternativa de completar o déficit com fundações livres sorteadas seria uma mutação aleatória, ruim em custo. |
| Onde centralizar R | Numa turbina do doador sorteada ao acaso | As zonas são polígonos com buracos e a densidade é de 1 a 2%. Um R sorteado em qualquer ponto do retângulo envolvente cairia com frequência em áreas vazias. No artigo, R é sorteado no terreno inteiro. |
| Tamanho de R | 10% a 50% da largura e da altura da zona | Valor inicial, ainda não calibrado. Fica em `BLOCK_SWAP_MIN_FRACTION` e `BLOCK_SWAP_MAX_FRACTION`. |
| Tentativas por zona | 10 | Valor inicial (`BLOCK_SWAP_MAX_ATTEMPTS`). O artigo não fixa T. |
| Distância mínima σ | Não usada | O problema não tem essa restrição: o grid de 160 m é menor que o rotor de 240 m, e a proximidade só é penalizada pelo wake. |
| Ajuste das bordas de R | Não feito | No artigo, serve para manter a distância de segurança, que não existe aqui. |
| Zonas | Tratadas de forma independente, cada uma com seu R | Nenhuma turbina troca de zona, e cada zona termina com exatamente τ_z turbinas. |
| Ids das turbinas | Um único contador para o filho inteiro, a partir de `fixd.size()` | O wake ignora pares de turbinas com o mesmo id. O mesmo contador foi aplicado ao crossover atual no commit `ae245d2c0`. |
| Atributos das turbinas | Copiados do pai (diâmetro, altura, Ct, potência, zona, fundação) | Só o id é renumerado. |
| Filho igual à base | Devolvido sem avaliação | Não há solução nova, então não faz sentido gastar orçamento. Em 300 cruzamentos de teste, isso ocorreu de 3% a 4% das vezes. |
| Escolha do operador | 5º argumento da linha de comando, `onepoint` (padrão) ou `blockswap` | Os dois operadores rodam no mesmo binário. Um valor inválido encerra com erro. |
| Crossover atual | Mantido sem alteração | O Block Swap é uma alternativa, e não um substituto. |

---

## 5. Uso

```bash
cd scripts
./exe/moead A ../dataset/moead/A/1/ 30 10 blockswap   # Block Swap
./exe/moead A ../dataset/moead/A/1/                   # onepoint (padrão)
```

O operador é o 5º argumento, então é preciso passar também o ângulo e o vento. O `log.txt` registra o crossover usado (`Crossover: blockswap`).

Os scripts de lote (`batch.sh` → `main.sh` → `moead.sh`) passam só dois argumentos e, por isso, rodam sempre com o `onepoint`.

### Arquivos alterados

| Arquivo | Mudança |
|---|---|
| `src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp` | Operador (novo) |
| `headers/global_modules/genetic_operators/block_swap_crossover.h` | Header (novo) |
| `Makefile` | Novo arquivo em `SRC_FILES` |
| `main.cpp` | Global `crossover_type`, leitura de `argv[5]`, validação e impressão no log |
| `headers/global_modules/generate_initial_population/population.h` | `extern std::string crossover_type` |
| `src/metaheuristics/moead/moead.cpp` | Escolha da função de crossover pelo `crossover_type` |

---

## 6. Validação

O operador foi testado com o código real nas instâncias A (2 zonas, 34 turbinas) e 240 (1 zona, 77 turbinas). Em cada instância foram 300 cruzamentos: metade entre soluções aleatórias e metade entre uma solução e uma cópia dela com 10 mutações, o que imita pais vizinhos no MOEA/D.

Em todos os filhos, sem nenhuma falha:
- o número de turbinas por zona era exatamente τ_z;
- o `grid` batia com as turbinas (`grid[z][index] == id`), e as coordenadas batiam com as fundações;
- os ids eram únicos e maiores ou iguais a `fixd.size()`;
- os atributos físicos estavam corretos;
- toda turbina ocupava uma fundação ocupada na base ou no doador;
- o fitness conferia ao ser recalculado.

| Instância | Filhos iguais à base (sem avaliação) | Filhos avaliados | Turbinas vindas só do doador |
|---|---|---|---|
| A | 9/300 | 291/300 (1 avaliação cada) | 14,7% |
| 240 | 13/300 | 287/300 (1 avaliação cada) | 13,2% |

O MOEA/D completo também rodou até o fim com orçamento reduzido (50 mil avaliações) com os dois operadores. Uma única rodada curta não permite comparar o desempenho deles.

---

## 7. Próximos passos

- **Comparação experimental:** `onepoint` contra `blockswap`, com 20 execuções por instância, em instâncias de 1 e 2 zonas, comparando o hipervolume.
- **Lote:** para isso, o `moead.sh` precisa aceitar o tipo de crossover e gravar numa pasta separada.
- **Calibração:** o tamanho de R e o número de tentativas usam valores iniciais e podem ser ajustados.
- **Seleção adaptativa:** os dois operadores têm comportamentos diferentes (faixa até a borda contra ilha interna), o que é útil para seleção adaptativa de operadores mesmo que nenhum seja melhor sozinho.

---

## 8. Texto sugerido para a metodologia

> *"Além do crossover de um ponto, que, pela ordenação das fundações candidatas em varredura do grid, recombina faixas contíguas de cada zona, foi implementado o operador Block Swap Crossover, adaptado de Tran et al. (2013) e também utilizado por Li et al. (2017). Para cada zona, o operador sorteia um retângulo centrado numa turbina do segundo pai e substitui as turbinas do primeiro pai dentro dele pelas do segundo, desde que o segundo pai tenha ao menos tantas turbinas na região quanto o primeiro. O excedente é removido aleatoriamente entre as turbinas do primeiro pai fora do retângulo, o que preserva o número exato de turbinas por zona e não introduz posições novas. Diferentemente do trabalho original, não há ajuste de distância mínima, pois o problema considerado não impõe essa restrição."*

O operador deve ser apresentado como uma **estratégia de variação informada pelo domínio**, e não como superior ao crossover de um ponto: nenhum dos trabalhos de origem isola o efeito do operador, e a vantagem dele sobre um corte em faixa precisa ser demonstrada experimentalmente.
