# EPET060-ER2-lectures

Este repositórios contém notas de aula para o curso de Engenharia de Reservatórios de Petróleo II do curso de Graduação em Engenharia de Petróleo da Universidade Federal de Alagoas (UFAL). Este repositório é utilizado para desenvolvimento e atualizações das aulas, trabalhos e projetos propostos ao longo do curso.

Todas as notas de aula e exemplos estão na linguagem [Julia](https://julialang.org) recentemente desenvolvida,  *opensource* (código aberto) e disponível para qualquer sistema operacional (Linux, OSX e Windows). A linguagem Julia foi desenvolvida explicitamente para cálculo numérico com foco em operações de vetores e matrizes, álgebra linear e execução paralela distribuída, com uma sintaxe semelhante ao [MATLAB](https://www.mathworks.com) porém projetada para alto desempenho, onde os programas são compilados automaticamente para código nativo eficientemente via LLVM. Além disso, é uma linguagem dinamicamente tipada (like a python) e tem bom suporte para uso interativo, podendo ser compilada separadamente.

Por se tratar de um curso que envolve problemas que apresentam um certo grau de dificuldade e elevada quantidade de cálculos, faz necessário que os discentes além de terem noções de estatística, compreensão sobre as principais propriedades dos fluidos (i.e. PVT, FVF, Rs, viscosidade, densidade, gas (id)real) e das rochas (Pc, porosidade, permeabilidade relativa, etc.) e dos tipos de reservatorios (BLKOIL, Dry — Wet Gas, retrogado, etc), é mandatório ter conhecimentos (mínimos) de MS Excel/Libre Calc, MATLAB, Julia, Python, etc. para resolverem alguns exercícios propostos ao longo do curso, bem como do projeto. Principalmente, referente aos projetos (feitos em grupos de alunos) é desejável que cada grupo escolha uma linguagem de programação (Python, MATLAB, Excel, ...) para o desenvolvimento, complemento e execução dos estudos/exercícios propostos.

Quando proposto projetos, estes devem ser apresentados como um breve relatório descrevendo o problema, abordagem computacional, código desenvolvido, resultados e conclusões. Preferencialmente (não obrigatoriamente), o texto pode ser escrito em LATEX para ser facilmente anexado a essas notas.

# Objetivos

Capacitar os alunos no uso e interpretação de cálculos aceitos e trivialmente realizados pela indústria, através da análise de dados de produção de poços - bem como serem capazes de realizar tais análises analiticamente (ou por auxílio de computador). Deduzir a relação de balanço de material para um líquido levemente compressível (óleo) na presença de outras fases (gás e água), bem como a relação de balanço de material para um gás seco. Analisar dados de produção (dados taxa-tempo ou taxa-pressão) para obter o volume do reservatório e estimativas das propriedades do reservatório para sistemas de reservatório de gás e líquido e, deste modo, fazer previsões de desempenho para tais sistemas. Conhecer a estrutura de produção secundária de reservatórios petrolíferos, suas distintas fases e comportamento para o entendimento da exploração e produção dos hidrocarbonetos.

# Ementa

1. Estimativa de reservas.
2. Balanço de materiais em reservatórios de gás e previsão de comportamento.
3. Balanço de materiais em reservatórios de óleo e previsão de comportamento.
4. Curvas de declínio de produção: previsão de recuperações futuras e ajuste de histórico.
5. Métodos de recuperação de petróleo e gás natural.

# Bibliografia

A maioria dos livros didáticos de engenharia de reservatórios. Existem vários livros que fazem uma boa referência sobre o assunto, que estão incluídos abaixo.

- DAKE, L. P. Fundamentals of Reservoir Engineering, Elsevier, New York, 1978.
- TOWLER, B. F., Fundamental Principles of Reservoir Engineering, SPE Textbook series, vol. 8, 2002.
- LAKE, L. W., Enhanced oil recovery. Prentice Hall, 1996.
- BLUNT, M., Reservoir Engineering, The Imperial College Lectures, 2015.
- ROSA, A. J. & CARVALHO, R. S. Engenharia de Reservatórios de Petróleo, Interciência, Rio de Janeiro, 2006.
- ROSA, A. J. & CARVALHO, R. S. Previsão de Comportamento de Reservatórios de Petróleo, Interciência, Rio de Janeiro, 2002.
- THOMAS, J. E. Fundamentos de engenharia de petróleo. 2. ed. Rio de Janeiro: Interciência, 2004.
- AMYX, J. W., BASS JR., D. M. & WHITING, R. L. Petroleum Reservoir Engineering, McGraw-Hill, New York, 1960.
- BEDRIKOVETSKY, P. G. Mathematical Theory of Oil & Gas Recovery, Kluwer Academic Publishers, London-Boston-Dordrecht, 1993.

