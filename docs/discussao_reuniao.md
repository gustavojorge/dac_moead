# Discussão para a reunião

## 1. Filho idêntico a um pai deve gastar orçamento?

### Situação atual

O orçamento é de 1.000.000 avaliações. O contador (`countRevalue`) soma 1 a cada chamada de `calculate_cost` + `calculate_power`. Os dois crossovers tratam de forma diferente o filho que sai idêntico a um pai:

| Operador | Filho idêntico a um pai |
|---|---|
| `onepoint` ([crossover.cpp](../src/moead/src/global_modules/genetic_operators/crossover/crossover.cpp)) | Avalia e conta 1 |
| `blockswap` ([block_swap_crossover.cpp](../src/moead/src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp)) | Devolve a base, sem avaliar e sem contar |

### Proposta: não contar, nos dois operadores

1. **O orçamento mede esforço computacional, e o custo está em avaliar.** O custo e a potência de um filho idêntico a um pai já são conhecidos. Avaliar de novo não gera informação nova.
2. **O contador deve significar "avaliações realizadas".** Cobrar uma avaliação que não aconteceu descaracteriza o contador. Avaliar só para manter a contagem desperdiça CPU.
3. **A regra precisa ser a mesma nos dois operadores.** Hoje só o Block Swap a aplica, o que o favorece na comparação: com o mesmo orçamento, ele desperdiça menos avaliações.
4. **O número de filhos descartados vira uma métrica de cada operador,** em vez de ficar escondido no contador.

### Dados medidos

Pais evoluídos: soluções vizinhas do arquivo Pareto depois de 100 mil avaliações, nas instâncias A (2 zonas) e 240 (1 zona). Faixa entre dois arquivos, um gerado por cada operador. Foi uma rodada por configuração, então os valores são indicativos.

| Filho com custo e potência iguais aos de um pai | A | 240 |
|---|---|---|
| `onepoint` | 24–36% | 8–15% |
| `blockswap` | 46–48% | 21–29% |

- Os pares foram vizinhos do arquivo Pareto, e não os pais reais do MOEA/D. Em execuções completas, as taxas podem ser diferentes.
- A contagem foi verificada: em 4 execuções de 100 mil avaliações, `countRevalue` foi igual ao número de chamadas a `calculate_cost` e ao de chamadas a `calculate_power`.

### Alternativas

| Opção | Vantagem | Desvantagem |
|---|---|---|
| **A. Não contar nos dois (proposta)** | Critério único, sem desperdício | Altera o `crossover.cpp` atual |
| B. Manter como está | Nada a alterar | Critérios diferentes, favorece o Block Swap e precisa ser explicado na dissertação |
| C. Contar nos dois | Critério único, mesmo número de gerações | Desperdiça CPU (se avalia de novo) ou o contador deixa de significar avaliações (se cobra sem avaliar) |

### O que a proposta muda

- `crossover.cpp`: devolver o pai, sem avaliar, quando o filho tiver o mesmo layout dele. A regra já existe no Block Swap.
- O `onepoint` rende um pouco mais por orçamento.
- O invariante de contagem do [CLAUDE.md](../CLAUDE.md) passa a dizer que filho idêntico a um pai não é avaliado nem contado.
- Os resultados de referência em `dataset/reference/moead/` não são afetados, porque já não são comparáveis (vento sorteado e bug dos ids).

### Decisão

_A preencher depois da reunião._

---

## 2. O Block Swap cria hiperparâmetros a calibrar?

### Situação atual

O operador usa valores fixos que eu escolhi, sem calibração. Eles ficam em [block_swap_crossover.cpp](../src/moead/src/global_modules/genetic_operators/crossover/block_swap_crossover.cpp) e em [main.cpp](../src/moead/main.cpp):

| Parâmetro | Valor | Significado |
|---|---|---|
| `BLOCK_SWAP_MIN_FRACTION` | 0,1 | Menor fração da largura e da altura da zona que o retângulo R pode ter |
| `BLOCK_SWAP_MAX_FRACTION` | 0,5 | Maior fração |
| `BLOCK_SWAP_MAX_ATTEMPTS` | 10 | Tentativas de achar um R aceito, por zona |
| `crossover_type` | `onepoint` ou `blockswap` | Escolha do operador (categórico) |

### Ponto

**Sortear dentro de uma faixa não elimina o parâmetro.** A largura e a altura de R são sorteadas entre 10% e 50%, mas os limites 10% e 50% foram escolhidos. [10%, 50%], [30%, 70%] e [5%, 20%] dão operadores diferentes.

- O sorteio dá uma distribuição de tamanhos de bloco e tende a deixar o operador menos sensível ao valor exato dos limites. Não dispensa o ajuste (por exemplo, com o IRACE).
- As frações são de dimensões diferentes e não se somam: o maior R tem 50% da largura e 50% da altura, ou seja, 25% da área da caixa da zona.

### Opções

| Opção | Vantagem | Desvantagem |
|---|---|---|
| A. Manter a faixa [10%, 50%] e as 10 tentativas como valores iniciais, e calibrar depois | Nada a alterar agora, e a comparação básica não espera a calibração | 3 parâmetros numéricos a calibrar |
| B. Calibrar só a faixa e fixar as tentativas em 10 | Menos parâmetros | Não testamos se o número de tentativas importa |
| C. Sortear a fração em toda a faixa possível, sem limites escolhidos | Sem parâmetros de faixa | Aceita blocos muito pequenos ou muito grandes, e o desempenho pode piorar |

Dimensionar o bloco por k turbinas do doador, em vez de por área, foi considerado e deixado de lado por ora. Ele também teria um parâmetro de faixa (a de k).

### Sugestão

Começar pela opção A: o primeiro experimento compara `onepoint` e `blockswap` com os valores iniciais. Se o Block Swap se mostrar promissor, ou se a taxa alta de filhos devolvidos iguais à base (21% a 48% no meu teste) for um problema, calibrar os parâmetros em seguida.

### Decisão

_A preencher depois da reunião._
