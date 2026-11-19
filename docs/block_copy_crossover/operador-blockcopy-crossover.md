# Operador BlockCopy Crossover no MoWFLOP

Este documento descreve o operador **BlockCopy Crossover** como foi implementado neste repositório: de onde ele vem, o que mudou na adaptação, as decisões tomadas e como ele foi validado.

Implementação: [block_copy_crossover.cpp](../../src/moead/src/global_modules/genetic_operators/crossover/block_copy_crossover.cpp).

> **Base bibliográfica.** O operador é baseado no **BlockCopy Crossover de Mayo e Zheng [1]**, seção III (III.A e III.B), cujo PDF é [mayo_zheng.pdf](mayo_zheng.pdf). Nenhuma outra fonte foi usada para o operador. O que foi mantido do artigo e o que foi adaptado ao nosso problema está na seção 1.2, e a referência completa está na seção 9.

---

## 1. Origem na literatura

### 1.1 Mayo e Zheng, a base do operador [1]

O operador foi proposto por Michael Mayo e Chen Zheng (University of Waikato) em *BlockCopy-based operators for evolving efficient wind farm layouts* ([mayo_zheng.pdf](mayo_zheng.pdf)). O artigo propõe o BlockCopy em duas versões: como mutação (seção III.A) e como crossover (seção III.B). Nossa implementação segue o crossover.

A ideia central, nas palavras dos autores, é que "it is the relative configurations of turbines rather than their individual absolute positions on the layouts that count" (Resumo). Por isso, em vez de mexer em turbinas isoladas, o operador copia um **padrão** de turbinas de uma parte do parque para outra, preservando a arrumação relativa entre elas.

O cenário do artigo é diferente do nosso:

| | Mayo e Zheng [1] | Este repositório |
|---|---|---|
| Posições | Coordenadas contínuas num terreno retangular | Fundações discretas, num grid de ~160 m dentro de polígonos |
| Zonas | Uma área (com obstáculos em dois cenários) | Até 3 zonas, cada uma com número fixo de turbinas |
| Restrições | Distância mínima entre turbinas (8 raios do rotor) | Só as da representação (sem distância mínima) |
| Objetivos | Um: custo da energia (função de Kusiak e Song) | Dois: custo das fundações e potência |
| Algoritmo | Estratégia evolutiva, (1+1)-ES e (5,10)-ES | MOEA/D |
| Blocos | Grid fixo de blocos de ~1 km × 1 km | Retângulo sorteado em cada cruzamento |

O operador do artigo:

1. **Mutação (III.A):** divide o parque em blocos aproximadamente quadrados, sorteia um bloco de origem e um de destino, apaga as turbinas do destino e copia para ele as turbinas da origem. Só a posição absoluta muda; a configuração relativa é mantida ("essentially a copy/translate of the source block to a new position").
2. **Turbina que não cabe:** se uma turbina copiada colidir com um obstáculo ou ficar perto demais de outra, ela simplesmente não é colocada ("the copied turbine is simply not placed").
3. **Correção global:** como os blocos podem ter quantidades diferentes de turbinas, turbinas são colocadas ou removidas ao acaso até o total voltar ao valor fixo ("turbines are either randomly placed onto the layout, or randomly purged from it").
4. **Crossover (III.B):** funciona como a mutação, mas com dois layouts: o bloco é copiado/transladado de um pai para uma cópia do outro, e cada crossover gera **um** filho. Os autores limitam a troca a **um único bloco**, porque trocar metade dos blocos se mostrou destrutivo demais nos testes iniciais.

### 1.2 O que mantivemos do Mayo e Zheng e o que adaptamos

**Mantido do artigo:**
- **Copiar e transladar um padrão:** as turbinas de um bloco do doador são copiadas para outra posição da base, todas com o mesmo deslocamento. A arrumação relativa entre elas é preservada.
- **Limpar o destino antes de colar:** as turbinas da base dentro do bloco de destino são removidas.
- **Turbina que não cabe é descartada:** se a posição deslocada não for uma fundação da zona, ou estiver ocupada, a turbina não é colocada.
- **Correção global ao acaso:** o excesso é removido ao acaso e o déficit é preenchido em posições livres sorteadas.
- **Um único bloco por cruzamento** (em cada zona).

**Adaptado por nós:**
- **Bloco sorteado, e não um grid fixo de blocos de 1 km:** a origem é um retângulo centrado numa turbina do doador, com largura e altura entre 10% e 50% das dimensões da zona. É o mesmo critério do Block Swap, então não há parâmetro novo. Os valores são nossos e ainda não foram calibrados.
- **Destino:** uma fundação qualquer da zona, diferente do centro da origem. O artigo não diz como o destino é escolhido no crossover.
- **Deslocamento nunca zero:** com deslocamento zero, o BlockCopy repetiria o Block Swap.
- **Translação exata no grid de fundações:** como o destino é uma fundação, e as fundações de cada zona formam um grid regular, as turbinas deslocadas caem exatamente sobre fundações. No artigo, as posições são contínuas.
- **Excesso removido só entre as turbinas da base:** o padrão colado fica intacto. O artigo diz apenas que as turbinas são removidas ao acaso.
- **Sem distância mínima e sem obstáculos:** o nosso problema não tem essas restrições. A turbina só é descartada quando cai fora da zona ou numa fundação já ocupada.
- **Zonas:** o operador é aplicado zona por zona, cada uma com sua origem e seu destino, e nenhuma turbina troca de zona.
- **Dois filhos por cruzamento:** o MOEA/D chama `crossover(A, B)` e `crossover(B, A)`, como nos outros operadores. O artigo gera um filho.
- **Probabilidade de crossover:** a do MOEA/D (0,5).
- **Filho idêntico à base:** devolvido sem avaliar e sem contar, como no Block Swap. A decisão está em discussão em [discussao_reuniao.md](../discussao_reuniao.md).

### 1.3 O que o artigo mostra e o que não mostra

- **A evidência forte é da mutação.** Numa (1+1)-ES com 2.000 avaliações, a BlockCopy mutation superou o TDA e a perturbação aleatória nos quatro cenários, com 30 execuções e teste de Wilcoxon (seção V).
- **O crossover não teve o mesmo resultado.** Numa (5,10)-ES em que o BlockCopy crossover era o único operador, "the population-based approach fails generally to improve on the (1+1)-ES algorithm" (seção VI). Os autores atribuem isso ao orçamento de 2.000 avaliações, que favorece a busca de uma única solução. Aqui o orçamento é de 10⁶ avaliações.
- **O artigo é mono-objetivo.** Não há evidência do operador num problema biobjetivo.
- **No nosso problema, a translação preserva só parte do padrão.** O wake depende das posições relativas, que são preservadas. O custo depende da fundação (profundidade e custo próprios), então o padrão colado tem outro custo.

---

## 2. Diferença em relação aos outros crossovers do repositório

| | `onepoint` | `blockswap` (Tran et al., 2013) | `blockcopy` (Mayo e Zheng [1]) |
|---|---|---|---|
| O que recombina | Uma faixa de índices de fundação | Um retângulo, na **mesma posição** nos dois pais | Um retângulo do doador, **colado em outra posição** da base |
| Posição das turbinas do doador no filho | A mesma | A mesma | Deslocada por Δ |
| Turbinas em fundações fora dos dois pais | Não | Não | Sim (as deslocadas e as do preenchimento) |
| Déficit de turbinas | Completado com o outro pai | Evitado: o retângulo é rejeitado e sorteado de novo | Preenchido ao acaso |

```
 blockswap (mesma posição)         blockcopy (outra posição)
 doador         filho              doador         filho
 . [B B] .      . [B B] .          [B B] . .      . . . .
 . [. B] .      . [. B] .          [. B] . .      . . . .
 . . . .        . . . .            . . . .        . [B B] .
 . . . .        . . . .            . . . .        . [. B] .
```

---

## 3. Algoritmo implementado

A interface é a mesma dos outros crossovers: `Solution block_copy_crossover(Solution &base, Solution &doador)`. O MOEA/D gera dois filhos, `crossover(A, B)` e `crossover(B, A)`, e cada chamada sorteia sua própria origem e seu próprio destino.

Uma única vez, na primeira chamada:

0. **Mapa do grid:** para cada zona, calcula o passo do grid em x e em y (a menor distância entre coordenadas distintas) e monta o mapa `(coluna, linha) → fundação`.

Para cada zona `z`, com τ_z turbinas:

1. **Caso trivial:** se o doador ou a base não tiverem turbinas na zona, a zona do filho fica igual à da base.
2. **Origem:** calcula o retângulo envolvente das fundações da zona (largura W e altura H), escolhe ao acaso uma turbina do doador como centro e sorteia a largura `W · U(0,1; 0,5)` e a altura `H · U(0,1; 0,5)`, de forma independente.
3. **Destino:** sorteia uma fundação qualquer da zona, diferente da fundação do centro. O retângulo de destino tem o mesmo tamanho da origem e é centrado nessa fundação.
4. **Deslocamento:** `Δx = x_destino − x_centro` e `Δy = y_destino − y_centro`.
5. **Limpeza:** as turbinas da base dentro do destino são removidas; as de fora são mantidas.
6. **Colagem:** cada turbina do doador dentro da origem vai para `(x + Δx, y + Δy)`. A posição é procurada no mapa do grid. Se não houver fundação ali (fora da zona) ou ela já estiver ocupada, a turbina é descartada.
7. **Excesso:** enquanto o total passar de τ_z, remove ao acaso uma turbina que veio da base. As coladas não são removidas.
8. **Déficit:** enquanto faltar turbina, coloca uma turbina numa fundação livre sorteada da zona.

Depois de todas as zonas:

9. **Ids:** as turbinas do filho são numeradas com um único contador, a partir de `fixd.size()`.
10. **Filho igual à base:** se o conjunto de fundações ocupadas for igual ao da base, a base é devolvida **sem nova avaliação**.
11. **Avaliação:** caso contrário, calcula custo e potência, insere no arquivo Pareto, incrementa `countRevalue` e salva os snapshots, como os outros operadores.

### Exemplo (τ_z = 6)

```
 doador (B)               base (A)                 filho
 ┌─────┐ . . . .          A . . ┌─────┐ .          A . . ┌─────┐ .
 │ B B │ . . B .          . . A │ A . │ .          . . . │ B B │ .
 │ . B │ . . . B          A . . │ . . │ A          A . . │ . B │ A
 └─────┘ B . . .          . A . └─────┘ .          . . . └─────┘ .
   origem                     destino
```

Neste exemplo:
- A origem, no doador, tem 3 turbinas. O destino, na base, tem 1 turbina, que é removida; ficam 5 turbinas da base.
- As 3 turbinas da origem são coladas no destino, na mesma arrumação. O filho teria 5 + 3 = 8 turbinas.
- As 2 excedentes saem ao acaso entre as 5 da base. O filho termina com 3 turbinas da base e as 3 coladas, total 6.
- Se uma turbina colada caísse fora da zona, ela seria descartada e o filho ficaria com uma a menos. Essa turbina seria então colocada numa fundação livre sorteada.

---

## 4. Decisões de projeto

| Decisão | Escolha | Motivo |
|---|---|---|
| Tamanho e centro da origem | Centrada numa turbina do doador, com 10% a 50% da largura e da altura da zona | Mesmo critério do Block Swap: a origem nunca sai vazia e não há parâmetro novo. Fica em `BLOCK_COPY_MIN_FRACTION` e `BLOCK_COPY_MAX_FRACTION`. |
| Destino | Fundação qualquer da zona, sorteada | O artigo não especifica. Sortear uma fundação garante que o destino está dentro da zona e que o deslocamento é um número inteiro de passos do grid. |
| Deslocamento zero | Proibido (sorteia de novo) | Com Δ = 0 o operador seria um Block Swap. Proibir só o zero não cria parâmetro novo. |
| Translação | Exata, pelo mapa `(coluna, linha) → fundação` | Verificado nas instâncias A, 240 e 10: as fundações de cada zona formam um grid regular, alinhado aos eixos. Não é preciso aproximar para a fundação mais próxima. |
| Turbina que não cabe | Descartada | Como no artigo. Aqui só acontece quando a posição sai da zona ou está ocupada. |
| Excesso | Removido ao acaso entre as turbinas da base | O padrão colado fica intacto. O artigo não diz de onde remove. |
| Déficit | Preenchido em fundações livres sorteadas da zona | Como no artigo. É diferente do Block Swap, que rejeita o retângulo para não criar posições aleatórias. |
| Distância mínima σ | Não usada | O problema não tem essa restrição: o grid de 160 m é menor que o rotor de 240 m, e a proximidade só é penalizada pelo wake. |
| Zonas | Tratadas de forma independente | Nenhuma turbina troca de zona, e cada zona termina com exatamente τ_z turbinas. |
| Ids das turbinas | Um único contador para o filho inteiro, a partir de `fixd.size()` | O wake ignora pares de turbinas com o mesmo id. |
| Atributos das turbinas | Copiados do pai (diâmetro, altura, Ct, potência, zona) | Só a fundação, as coordenadas e o id mudam. |
| Filho igual à base | Devolvido sem avaliação | Mesma regra do Block Swap; em discussão em [discussao_reuniao.md](../discussao_reuniao.md). |
| Escolha do operador | 5º argumento da linha de comando, `blockcopy` | Os três crossovers rodam no mesmo binário. Um valor inválido encerra com erro. |

---

## 5. Uso

```bash
cd scripts
./exe/moead A ../dataset/moead/A/1/ 30 10 blockcopy   # BlockCopy
./exe/moead A ../dataset/moead/A/1/ 30 10 blockswap   # Block Swap
./exe/moead A ../dataset/moead/A/1/                   # onepoint (padrão)
```

O operador é o 5º argumento, então é preciso passar também o ângulo e o vento. O `log.txt` registra o crossover usado (`Crossover: blockcopy`).

Os scripts de lote (`batch.sh` → `main.sh` → `moead.sh`) passam só dois argumentos e, por isso, rodam sempre com o `onepoint`.

### Arquivos alterados

| Arquivo | Mudança |
|---|---|
| `src/global_modules/genetic_operators/crossover/block_copy_crossover.cpp` | Operador (novo) |
| `headers/global_modules/genetic_operators/block_copy_crossover.h` | Header (novo) |
| `Makefile` | Novo arquivo em `SRC_FILES` |
| `main.cpp` | Aceita `blockcopy` no 5º argumento |
| `headers/main.h` | Include do header |
| `src/metaheuristics/moead/moead.cpp` | Escolha da função de crossover pelo `crossover_type` |

---

## 6. Validação

O operador foi testado com o código real, numa cópia à parte do repositório, nas instâncias A (2 zonas, 34 turbinas) e 240 (1 zona, 77 turbinas). Em cada instância foram 300 cruzamentos: metade entre soluções aleatórias e metade entre uma solução e uma cópia dela com 10 mutações, o que imita pais vizinhos no MOEA/D.

Em todos os filhos, sem nenhuma falha:
- o número de turbinas por zona era exatamente τ_z;
- o `grid` batia com as turbinas (`grid[z][index] == id`), e as coordenadas batiam com as fundações;
- os ids eram únicos e maiores ou iguais a `fixd.size()`;
- o fitness conferia ao ser recalculado;
- a contagem de avaliações estava correta (1 por filho avaliado, 0 por filho igual à base).

**Translação:** o mapa do grid foi conferido em todas as fundações (cada fundação é encontrada na própria posição). Em 20 mil deslocamentos aleatórios por zona, toda posição encontrada no mapa coincidiu com a posição deslocada.

| Instância | Filhos iguais à base | Turbinas do filho fora das fundações da base | Turbinas da origem coladas | Por cruzamento: na origem / coladas / removidas da base / preenchidas ao acaso |
|---|---|---|---|---|
| A | 0/300 | 21% | 71% | 7,5 / 5,3 / 1,4 / 2,0 |
| 240 | 0/300 | 22% | 62% | 16,9 / 10,4 / 0,6 / 6,7 |

As turbinas da origem que não são coladas caíram fora da zona. Elas são repostas pelo preenchimento aleatório, que funciona como uma mutação aleatória e é o principal efeito colateral do operador.

O MOEA/D completo também rodou com o BlockCopy na instância A por 60 segundos (~507 mil avaliações), sem erros. Uma rodada curta não permite comparar o desempenho com os outros operadores.

---

## 7. Próximos passos

- **Comparação experimental:** `onepoint`, `blockswap` e `blockcopy`, com 20 execuções por instância, comparando o hipervolume.
- **Lote:** para isso, o `moead.sh` precisa aceitar o tipo de crossover e gravar numa pasta separada.
- **Calibração:** o tamanho do retângulo usa valores iniciais e pode ser ajustado.
- **BlockCopy mutation:** é a versão com evidência forte no artigo e poderia ser testada como mutação.

---

## 8. Texto sugerido para a metodologia

> *"Como terceiro operador de recombinação, foi implementado o BlockCopy Crossover, adaptado de Mayo e Zheng [1]. O operador parte da premissa de que a configuração relativa entre turbinas, e não a posição absoluta de cada uma, determina as interações de esteira. Para cada zona, sorteia-se um retângulo centrado numa turbina do segundo pai e uma fundação de destino na zona; as turbinas do primeiro pai dentro do retângulo de destino são removidas, e as do segundo pai dentro do retângulo de origem são transladadas para o destino, preservando sua configuração relativa. Como as fundações candidatas formam um grid regular, a translação mapeia exatamente fundações em fundações; turbinas que deixariam a zona são descartadas. O número de turbinas por zona é restaurado como no trabalho original: o excedente é removido aleatoriamente entre as turbinas herdadas do primeiro pai e o déficit é preenchido em fundações livres sorteadas. Diferentemente do trabalho original, os blocos não seguem uma grade fixa e não há restrição de distância mínima, que o problema considerado não impõe."*

O operador deve ser apresentado como uma **estratégia de variação informada pelo domínio**, e não como superior aos outros crossovers: no artigo de origem, a evidência favorável é da versão de mutação, e a vantagem do crossover precisa ser demonstrada experimentalmente.

---

## 9. Referências

**[1]** M. Mayo, C. Zheng. *BlockCopy-based operators for evolving efficient wind farm layouts*. IEEE, 2015. **Fonte do operador** (seção III). PDF: [mayo_zheng.pdf](mayo_zheng.pdf).

O PDF traz no rodapé `978-1-4799-8697-2/15 ©2015 IEEE`, mas não informa o nome do evento, as páginas nem o DOI. Esses dados precisam ser confirmados antes de citar o artigo.

```bibtex
@inproceedings{mayo2015blockcopy,
  author    = {Mayo, Michael and Zheng, Chen},
  title     = {BlockCopy-based operators for evolving efficient wind farm layouts},
  publisher = {IEEE},
  year      = {2015},
  isbn      = {978-1-4799-8697-2},
  note      = {Evento, páginas e DOI a confirmar}
}
```
