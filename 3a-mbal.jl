### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 78d10c30-d382-11ef-2a86-1d3c01e82daa
using Plots, PlutoUI, HypertextLiteral, LsqFit, MultiComponentFlash, PlutoTeachingTools,  Printf

# ╔═╡ 369a85ef-d489-4f2c-83aa-0348897aafda
begin
	# auxilary helpers
	struct TwoCols{w, L, R}
	leftcolwidth::w
    left::L
    right::R
	end
	
	function Base.show(io, mime::MIME"text/html", tc::TwoCols)
	    write(io, """<div style="display: flex;"><div style="flex: $(tc.leftcolwidth * 100)%;padding-right:1.5em;">""")
	    show(io, mime, tc.left)
	    write(io, """</div><div style="flex: $((1.0 - tc.leftcolwidth) * 100)%;padding-left=1.5em;">""")
	    show(io, mime, tc.right)
	    write(io, """</div></div>""")
	end
end

# ╔═╡ ce8cf761-58df-4769-af0b-b71efdf7d2f6
@htl"""
<style>
	.edit_or_run { position: absolute; }

	pluto-output code pre {
		font-size: 90%;
	}
	pluto-input .cm-editor .cm-content .cm-scroller {
		font-size: 90%;
	}
	pluto-tree {
		font-size: 90%;
	}
	pluto-output, pluto-output table>tbody td  {
	 	font-size: 100%;
	    line-height: 1.8;
	}
	.plutoui-toc {
	    font-size: 90%;
		border-radius: 50px;
  		background:rgba(209, 202, 202, 0.84);
	}
	.admonition-title {
	 	color: var(--pluto-output-h-color) !important;
	}
	pluto-output h1, pluto-output h2, pluto-output h3 {
	    border: none !important;
		line-height: 1.3 !important;
	}
	pluto-output h1 {
	    font-size: 200% !important;
	    /*border-bottom: 3px solid var(--pluto-output-color)!important;*/
	}
	pluto-output h2 {
	    font-size: 175% !important;
	    padding-top: 0.75em;
	    padding-bottom: 0.5em;
	}
	pluto-output h3 {
	    font-size: 130% !important;
	    padding-top: 0.3em;
	    padding-bottom: 0.2em;
	}
	img:not(picture *) {
	    width: 80%;
	    display: block;
	    margin-left: auto;
	    margin-right: auto;
	}
	title {
	    font-size: 200%;
	    display: block;
	    font-weight: bold;
	    text-align: left;
	    margin: 2em 0 0 0;
	}
	subtitle {
	    font-size: 140%;
	    display: block;
	    text-align: left;
	    margin: 0 0 1.5em 0;
	}
	author {
	    font-size: 120%;
	    display: block;
	    text-align: left;
	    margin: 0 0 1.5em 0;
	}
	email {
	    font-size: 100%;
	    display: block;
	    text-align: left;
	    margin: -1.8em 0 2em 0;
	}
	hr {
	color: var(--pluto-output-color);
	}
	semester {
	    font-size: 100%;
	    display: block;
	    text-align: left;
		padding-bottom: 0.5em;
	}
	blockquote {
		padding: 1.3em !important;
	}
	li::marker {
	    font-weight: bold;
	}
	li {
	    margin-top: 0.5em !important;
	    margin-bottom: 0.5em !important;
	}
	.small {
	    font-size: 80%;
	}
	.centered-image {
	    display: block;
	    margin-left: auto;
	    margin-right: auto;
	    margin-top: 1em;
	    margin-bottom: 1em;
	    width: 90%;
	}
	.hanged {
		padding-left: 3em;
		text-indent: -3em;
	}
	pluto-output details summary {
		font-weight: normal !important;
	}
	
</style>
"""

# ╔═╡ 785e1863-8001-4885-8829-9f98419d6cd9
TableOfContents(title="Índice")

# ╔═╡ 465c140f-4729-416c-b51e-901836d545e8
begin
	Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)
	@htl"""
	<button onclick="present()">Apresentar</button>
	<div style="margin-top:3em;margin-bottom:7em;">
	</div>
	<title>Engenharia de Reservatórios 2</title>
	<subtitle>MBAL para reservatórios de gás natural</subtitle>
	<author>Jonathan da Cunha Teixeira</author>
	<email><a href="mailto:jonathan.teixeira@ctec.ufal.br">jonathan.teixeira@ctec.ufal.br<a/></email>
	<semester>Engenharia de Petróleo<br>Universidade Federal de Alagoas</semester>
	<hr style="border-top:8px dashed;margin:2em 0em;"/>
	"""
end

# ╔═╡ 5c641485-6bea-498e-874b-c0fc0ed99743
md"""
# Balanço de Material (*MBAL*) de reservatórios de gás natural

**Objetivo**: Derivar as relações do balanço de material para gases compressíveis (gás natural) na presença de outras fases (água e condensado).

## Idéia geral

Partindo de equações simples (lineares) contabilizar todas as "massas" que entram/saem/geradas no reservatório, e desta forma, quantificar a reservas de HC, *medidos sob condições-padrão*, e realizar previsão do comportamento de reservatórios ao longo do tempo (na sua forma padrão o "tempo" é uma medida de depletação do reservatório). Além disso, quando existir, "determinar" o influxo de água proveniente de aqüíferos associados.

Esta análise é muitas vezes chamada de análise/modelagem concentrada, onde o reservatório é considerado como sendo um "tanque" em que adicionamos ou subtraimos uma certa quantidade de volume, e monitoramos a variação da pressão. Portanto, a pressão (média) do reservatório é **NECESSÁRIA**, porem é o elo fraco! Por que?...


![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal.png?raw=true)

Todos os volumes medidos/quantificados no MBAL estão sob condições padrão (*standart*) porém resultados/respostas da análise estão em condições de reservatório (RB/cf)

Equação de balanço de material: 

$Volume\ a\ P_i -\ Volume\ a\ P = Volume\ dos\ fluidos\ produzidos$
"""

# ╔═╡ 856145bf-8853-4d07-81dd-fcecbc45c422
md"""
# Reservatório de gás (natural) seco

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas.png?raw=true)

Inicialmente temos:

$$G = \frac{V_r\phi N_G S_{gi}}{B_{gi}}$$

Do Balanço de material (MBAL):

$$Gas\ produzido = Gas\ inicial(p_i)\ -\ Gas\ remanescente(p)$$

Portanto:

$$G_p = G - \frac{V_r\phi N_G S_{gi}}{B_{g}} = G\left(1-\frac{B_{gi}}{B_g}\right)$$

Podemos converter a equação em termos de pressão média (parâmetro medido), partindo da EOS dos gases reais:

$$pV = ZnRT$$

sabemos que FVF do gás natural é:

$$B_g = \frac{V_{r}}{V_{std}}$$

Logo:

$$B_g = \frac{Z T p_{std}}{T_{std} p} = 0.0283\frac{ZT}{p}$$

Finalmente temos as seguintes versões da equação de balanço:

Forma Gp x Bg:

$$Gp = G\left(1-\frac{B_{gi}}{B_g}\right)$$

Forma Gp x p/z

$$Gp = G\left(1-\frac{Z_{i}p}{p_i Z}\right)$$
"""

# ╔═╡ 7a2c0627-5a20-42f4-b99f-ef9a37b00aed
let
	press = [1920, 1850, 1802, 1720, 1638, 1475]
	z = [0.8542, 0.8672, 0.8802, 0.8932, 0.9072, 0.9230]
	Gp = [0.0, 1.36, 2.41, 3.50, 4.95, 6.84]
	p_z = press ./ z
	Bg = (0.0283 * (170 + 460)) .* z ./ press
	Bgi_Bg = Bg[1] ./ Bg
	pred = range(0,maximum(p_z))

	mbal(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	fit = curve_fit(mbal, p_z, Gp, p0)
	scatter(p_z, Gp, label="Dados de Produção", xlabel="p/z", ylabel="Gp, MMMscf")
	plot!(p_z,mbal(p_z,fit.param),label="MBAL", c=:red, lw=2)
	plot!(pred,mbal(pred,fit.param),label="Forecast", c=:red, lw=2, ls=:dash)
	quiver!([500],[mbal(pred[1],fit.param)],quiver=([-450],[0]), c=:black, lw=3)
	annotate!(650, mbal(pred[1],fit.param), "OGIP")
	xlims!(pred[1], pred[end]+50)
	
end

# ╔═╡ 6532b23d-0a83-4cac-82ce-166dfef2ffb1
details(
	md"""**Exercício 1.** Um reservatório de gás volumétrico de 1100 acres é caracterizado por uma temperatura de 170°F, espessura do reservatório de 50 pés, porosidade média de 0,15, saturação inicial de água de 0,39. O histórico de produção de 5 anos é representado na tabela abaixo.

	| Tempo (a) | Pressão Res. (psia) | Fator compress. (-) | Produção acumulada de gás (MMMscf) |
	| --------- | ------------------- | ------------------- | ---------------------------------- |
	|   0       |  1920               |     0.8542          |       0.00     |  
	|   1       |  1850               |     0.8672          |       1.36     |  
	|   2       |  1802               |     0.8802          |       2.41     |  
	|   3       |  1720               |     0.8932          |       3.50     |  
	|   4       |  1638               |     0.9072          |       4.95     |  
	|   5       |  1475               |     0.9230          |       6.84     |  
	
	Determinar o OGIP com base nos dados do histórico.
	""",
	md"""
	Primeiramente calculamos via método volumétrico, a estimativa de reserva:

	$$OGIP = 43560\times\frac{Ah\phi(1-S_{wi})}{B_{gi}}=\frac{43560\times 1100\times 50\times 0.15\times(1-0.39)}{0.00793}$$
	
	$$OGIP= 2.764\times 10^{10}\ scf\approx 27.64\ MMMscf$$

	Com os dados do histórico de produção, realizamos o ajuste conforme o balanço de material... que a partir do gráfico Gp x p/z encontramos OGIP através do valor onde intercepta $p\rightarrow 0$
	
	$$G = OGIP\approx 23.78\ MMMscf$$

	![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-lin.png?raw=true)
	
	"""
)

# ╔═╡ b07ac32c-9017-4103-8614-3d9203413adb
md"""
# Reservatório de gás seco: Anormalmente pressurizado

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/compaction-drive.png?raw=true)

Normalmente ocorre em reservatórios *não-consolidados* e/ou quando a geopressão calculada/estimada, $p_{geo}$ é **diferente**  da medida, $p_{measure}$.

## Fatores da ocorrência

* Intercalações de arenitos-shales (formações muito permeaáveis e pouco permeáveis): a osmose entre águas com diferentes salinidades, o folhelho de "vedação" atua como uma membrana semipermeável nesta troca iônica. Se a água dentro da "vedação" for mais salina do que a água ao redor, a osmose causará uma pressão anormalmente alta e vice-versa.

* Desequílibrio no processo de compactação/soterramento da bacia sedimentar;

* Processo de maturação térmica do querogênio (vulcanismo): mudança de temperatura; um aumento na temperatura de um grau Fahrenheit pode causar um aumento na pressão de 125 psi em um sistema de água doce selado.

* (Presença de aquíferos em) regiôes montanhosas: intenso tectonismos pod causar sobrecarga das rochas adjacentes sobre a rocha reservatório.

## Ocorrência & causas

* Compactação e subsidência
* Danos ao poço
* Seismicidade|Reativação de falhas geológicas
* Fechamento do espaço poroso
* Elevado aumento da produção seguido de diminuição|dano

![dog](https://aoghs.org/wp-content/uploads/2018/03/thums-Long-Beach-Subsidence-Long-Beach-Historical-Society-1959-%E2%80%93-Roger-Coar-768x958.jpg)

## Derivação segundo Ramagost

Do Balanço de material (MBAL) temos:

$$Gas\ inicial(p_i) =\ Gas\ remanescente(p)\ -\ Gas\ produzido$$

sendo $V\neq V_i$, pois a variação da porosidade agora é considerável, ou na **mesma ordem de grandeza da compressibilidade do gás**

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-abnormal.png?raw=true)

$$V = V_i - \left[\Delta V_w + \Delta V_{f}\right]$$

sendo: $\Delta V_w$ Variação do volume da água conata, e $\Delta V_{f}$ Variação do volume poroso (formação)

### Variação do Volume poroso

Da correlação entre $p\rightarrow V$ (definição da compressibilidade isotérmica):

$$c_f = \frac{1}{V_{p}}\cdot\left(\frac{\partial V_{p}}{\partial p}\right)_T\approx\frac{1}{V_{pi}}\cdot\frac{\Delta V_{pi}}{\Delta p}$$

sendo: $\Delta p_i = p - p_i$, e $V_{p_i} = \frac{GB_{gi}}{(1-S_{wi})}$ volume poroso. Rearranjando temos então que: 

$$\Delta V_{pi} = c_f(p-p_i) \frac{GB_{gi}}{(1-S_{wi})}\approx-\Delta V_f$$

### Variação da água conata

Da correlação entre $p\rightarrow V$ (definição da compressibilidade isotérmica):

$$c_w = -\frac{1}{V_w}\cdot\left(\frac{\partial V_w}{\partial p}\right)_T\approx -\frac{1}{V_{w_i}}\cdot\frac{\Delta V_w}{\Delta p}$$

sendo: $\Delta p_i = p - p_i$, e $V_{w_i} = S_{wi}\frac{GB_{gi}}{(1-S_{wi})}$ volume de água conata. Rearranjando temos então que: 

$$\Delta V_{wi} = c_w(p_i-p) \frac{S_{wi}GB_{gi}}{(1-S_{wi})}$$

### Equação MBAL

Voltando para MBAL:

$$\underbrace{Gas\ inicial(p_i)}_{GB_{gi}} =\underbrace{Gas\ remanescente(p)}_{GB_{g}-\Delta V_w-\Delta V_f}\ -\ \underbrace{Gas\ produzido}_{G_pB_g}$$

Portanto:

$$GB_{gi} = \left[GB_{g} - c_w(p_i-p) \frac{S_{wi}GB_{gi}}{(1-S_{wi})} - c_f(p_i-p) \frac{GB_{gi}}{(1-S_{wi})}\right] - G_pB_{g} $$

$$GB_{gi} = (G - G_p)B_{g} + \frac{(p_i-p)GB_{gi}}{(1-S_{wi})}\left[c_wS_{wi} + c_f\right]$$

$$G\left[1 - \left(c_wS_{wi} + c_f\right)\frac{(p_i-p)}{(1-S_{wi})}\right]\frac{B_{gi}}{B_{g}} = G - G_p$$

#### MBAL versão $\frac{p}{Z}$

Portanto,

$$\frac{p}{Z}(1-c_{ewf}Δp) = \frac{p_i}{Z_i} - \frac{p_i}{Z_i}\frac{G_p}{G}$$

sendo: $c_{ewf} = \frac{S_{wi}c_w + c_f}{1-S_{wi}}$
"""

# ╔═╡ a01e834a-c916-4599-bf89-08500e7790bb
details(
	md"""**Exercício 2.** Num campo hipotético foram dados duas evidências que poderíam indicar que o campo é anormalmente pressurizado, são elas:

	- A pressão inicial registrada foi de 11.444 psia (valor muito alto, comparado às pressões iniciais dos reservatórios da área).
	
	- Reportaram que a compressibilidade da formação é de 19,5 $\mu$ip[^1], (superior à compressibilidade da água).
	
	Para constatar tais evidência de sobrepressão será necessário plotar os gráficos de balanço de materiais e comparar com a relação $\frac{p}{Z}_{corr}$, considere a saturação inicial de água igual a 0,22. ($c_w=$ 3$\mu$ip)

	| Data |  pressão (psia) | Fator compressibilidade (z) | Produção de gás (Bscf) |  Bg (cf/scf) |
	| ----- | ---------------------- | ---------------------- | -------- | ------ |
	| 25 Jan 1966 | 11444.133631 	| 1.629631 |  0.000000  | 0.002087 |
	| 1 Fev 1967 	| 10674.004243 	| 1.444006 |  9.924243  | 0.006997 |
	| 1 Fev 1968 	| 10131.047006 	| 1.442243 |  28.667006 | 0.039825 |
	| 1 Jun 1969 	| 9253.028651 	| 1.358651 |  53.628651 | 0.041589 |
	| 1 Jun 1970 	| 8574.031685 	| 1.311685 |  77.701685 | 0.053474 |
	| 1 Jun 1971 	| 7905.955657 	| 1.185657 |  101.375657| 0.064116 |
	| 1 Jun 1972 	| 7379.951424 	| 1.148661 |  120.311424| 0.074527 |
	| 1 Set 1973 	| 6846.994661 	| 1.143424 |  145.004661| 0.080189 |
	| 1 Ago 1974 	| 6387.903678 	| 1.099656 |  160.533678| 0.092732 |
	| 1 Ago 1975 	| 5827.015656 	| 1.042507 |  182.355656| 0.194581 |
	| 10 Jun 1976 | 5408.871680 	| 1.025678 |  197.601680| 0.124326 |
	| 1 Jun 1977 	| 4999.968671 	| 1.001671 |  215.628671| 0.127106 |
	| 1 Ago 1978 	| 4500.037507 	| 0.985245 |  235.777507| 0.136303 |
	| 1 Ago 1979 	| 4169.997245 	| 0.928680 |  245.897245| 0.142072 |

	---
	[^1] ip significa Inverse pressure/psi
	""",
	md"""
	Após o ajuste dos dados de produção (ajuste de histórico) sobre a equação:

	$$Gp = G\left(1-\frac{Z_{i}p}{p_i Z}\right)$$

	Obtermos OGIP gás seco ~ 648 Bscf

	Já para a equação:

	$$\frac{p}{Z}(1-c_{ewf}Δp) = \frac{p_i}{Z_i} - \frac{p_i}{Z_i}\frac{G_p}{G}$$
	
	temos que OGIP para um reservatório anormalmente pressurizado ~ 506 Bscf

	Conforme figura a seguir.
	"""
)

# ╔═╡ a0689d18-6906-42b7-b474-052819995933
let
	p = [11444.133631, 10674.004243, 10131.047006, 9253.028651, 8574.031685, 7905.955657, 7379.951424, 6846.994661, 6387.903678, 5827.015656, 5408.871680, 4999.968671, 4500.037507, 4169.997245]
	z = [ 1.629631, 1.444006, 1.442243, 1.358651, 1.311685, 1.185657, 1.148661, 1.143424, 1.099656, 1.042507, 1.025678, 1.001671, 0.985245, 0.928680]
	Gp = [0.000000, 9.924243, 28.667006, 53.628651, 77.701685, 101.375657, 120.311424, 145.004661, 160.533678, 182.355656, 197.601680, 215.628671, 235.777507, 245.897245]
	#Bg = [0.002087, 0.006997, 0.039825, 0.041589, 0.053474, 0.064116, 0.074527, 0.080189, 0.092732, 0.194581, 0.124326, 0.127106, 0.136303, 0.142072]

	# dry-gas
	p_z = p ./ z
	pred = range(0,maximum(p_z))
	mbal(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	fit = curve_fit(mbal, p_z, Gp, p0)
	OGIP = ceil(fit.param[1])
	# abnormal pressured
	cw = 3e-6
	cf = 19.5e-6
	swi = 0.22
	cef = (swi * cw + cf)/(1-swi)
	Δp = p[1] .- p
	p0 = [0.5, 0.5]
	fit2 = curve_fit(mbal, p_z.*(1.0 .- cef.*Δp), Gp, p0)
	OGIP2 = ceil(fit2.param[1])
	pred2 = range(0,maximum(p_z.*(1.0 .- cef.*Δp)))
	
	plt = plot(size=(800, 400), layout=grid(1, 2), leg=false)
	# plot1
	scatter!(plt[1], p_z, Gp, label="Dados de Produção", xlabel="p/z", ylabel="Gp, MMMscf", title="Reservatório de gás seco")
	plot!(plt[1], p_z,mbal(p_z,fit.param),label="MBAL", c=:red, lw=2)
	plot!(plt[1], pred,mbal(pred,fit.param),label="Forecast", c=:red, lw=2, ls=:dash)
	annotate!(plt[1], 2350, mbal(pred[1],fit.param), "OGIP = $OGIP Bscf")
	xlims!(plt[1], pred[1], pred[end]+500)
	# plot2
	scatter!(plt[2], p_z.*(1.0 .- cef.*Δp), Gp, label="Dados de Produção", xlabel="p/z|\$_{corr}\$", ylabel="Gp, MMMscf", title="Reservatório anormalmente press.")
	plot!(plt[2], p_z.*(1.0 .- cef.*Δp),mbal(p_z.*(1.0 .- cef.*Δp),fit2.param),label="MBAL", c=:red, lw=2)
	plot!(plt[2], pred2,mbal(pred2,fit2.param),label="Forecast", c=:red, lw=2, ls=:dash)
	annotate!(plt[2], 2350, mbal(pred2[1],fit2.param), "OGIP = $OGIP2 Bscf")
	xlims!(plt[2], pred2[1], pred2[end]+500)
end

# ╔═╡ eb272cd0-2472-4b48-a35a-6868e39de7c5
md"""
# Reservatório de gás condensado

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-condensate.png?raw=true)

O balanço de material é similar ao reservatório de gás natural seco, exceto que para este tipo de reservatório devemos contabilizar os volumes de condensado (os $(GE)_w$ e $(GE)_c$) definido como Volume de gás+condesáveis total ($G_{p_T}$)
"""

# ╔═╡ f52cd2f2-ec88-4ec3-a33c-f0f96909e776
let
	ns = 200
	ubar = 1e5
	# Pressure range
	p0 = 1*ubar
	p1 = 180*ubar
	# Temperature range
	T0 = 273.15 + 1
	T1 = 263.15 + 350
	# Define mixture + eos
	# names = ["Methane", "CarbonDioxide", "n-Decane"]
	names = ["Methane", "n-Butane", "n-Decane"]
	props = MolecularProperty.(names)
	mixture = MultiComponentMixture(props)
	eos = GenericCubicEOS(mixture)
	# Constant mole fractions, vary p-T
	z = [0.5, 0.2, 0.3]
	p  = range(p0, p1, length = ns)
	T = range(T0, T1, length = ns)
	cond = (p = p0, T = T0, z = z)
	
	m = SSIFlash()
	S = flash_storage(eos, cond, method = m)
	K = initial_guess_K(eos, cond)
	data = zeros(ns, ns)
	for ip = 1:ns
	    for iT = 1:ns
	        c = (p = p[ip], T = T[iT], z = z)
	        data[ip, iT] = flash_2ph!(S, K, eos, c, NaN, method = m)
	    end
	end
	
	contour(p./ubar, T .- 273.15, data, levels = 10, fill=(true,cgrad(:jet)))
	ylabel!("Pressão [Bar]")
	xlabel!("Temperatura [°Celsius]")
end

# ╔═╡ 5319239b-f361-4375-b49d-596b36a89641
md"""
## Gás equivalente ao condensado

$$(GE)_c = \frac{nRT_{std}}{p_{std}} = 133000\frac{\gamma_o}{MM_o}[scf/stb]$$

sendo: $\gamma_o| d_o$ densidade do condensado, $MM_o$ massa molecular do condensado, determinado por:

$$MM_o = \frac{5954}{API° - 8.811}=\frac{42.43 d_o}{1.008 - d_o}$$

## Gás produzido

$$G_{p_T} = G_p + (GE)_c + (GE)_w$$

onde: $(GE)_w$ vapor de água condensado $\approx 7390 scf/STB$

## Considerações duranto o gerenciamento da produção

- A perda de líquido no reservatório durante a produção é considerada imóvel (abaixo da saturação residual do condensado). Portanto perdido para a produção.

- Recirculação de gás (ciclo de injeção de gás) é uma possibilidade.

Para evitar a condensação no reservatório, mantemos a condição monofásica! Portanto, o gerenciamento deste tipo de reservatório é:

1. O projeto opera com a reinjeção + importação de gás natural seco até o *breakthrough*[^2] do gás seco (importado) e venda do condensado.
2. Após o *breakthrough*, a importação do gás e a reinjeção são cerceadas.
3. A operação continua como um reservatório de gás seco

Por que não  depletar normalmente o reservatório? Se os componentes que vaporizam primeiro são os **HC pesados dos quais não necessitamos**.

!!! info "Por que não  depletar normalmente o reservatório?"
	- A **composição do fluido do reservatório muda**, fazendo com que a mistura fique mais rica (HC de cadeia longa).
	- O diagrama de fase fica "inchado", i.e. **move-se para a direita**.

---

[^2] Há a possibilidade do *breakthrough* precoce de gás, devido à heterogeneidade do reservatório ou porque a viscosidade do gás seco ser menor que do gás úmido.
"""

# ╔═╡ 2133ca27-a8ea-4e93-88b0-412d91fe5efb
let
	ns = 200
	ubar = 1e5
	# Pressure range
	p0 = 1*ubar
	p1 = 180*ubar
	# Temperature range
	T0 = 273.15 + 1
	T1 = 263.15 + 350
	# Define mixture + eos
	# names = ["Methane", "CarbonDioxide", "n-Decane"]
	names = ["Methane", "n-Butane", "n-Decane"]
	props = MolecularProperty.(names)
	mixture = MultiComponentMixture(props)
	eos = GenericCubicEOS(mixture)
	# Constant mole fractions, vary p-T
	z = [0.5, 0.3, 0.2]
	z2 = [0.6, 0.3, 0.1]
	p  = range(p0, p1, length = ns)
	T = range(T0, T1, length = ns)
	cond = (p = p0, T = T0, z = z)
	cond2 = (p = p0, T = T0, z = z2)
	
	m = SSIFlash()
	S = flash_storage(eos, cond, method = m)
	K = initial_guess_K(eos, cond)
	data = zeros(ns, ns)
	for ip = 1:ns
	    for iT = 1:ns
	        c = (p = p[ip], T = T[iT], z = z)
	        data[ip, iT] = flash_2ph!(S, K, eos, c, NaN, method = m)
	    end
	end

	S2 = flash_storage(eos, cond, method = m)
	K2 = initial_guess_K(eos, cond)
	data2 = zeros(ns, ns)
	for ip = 1:ns
	    for iT = 1:ns
	        c = (p = p[ip], T = T[iT], z = z2)
	        data2[ip, iT] = flash_2ph!(S2, K2, eos, c, NaN, method = m)
	    end
	end
	
	plt = plot(size=(800, 400), layout=grid(1, 2), leg=false)
	contour!(plt[1], p./ubar, T .- 273.15, data, levels = 10, fill=(true,cgrad(:jet)))
	ylabel!(plt[1], "Pressão [Bar]")
	xlabel!(plt[1], "Temperatura [°Celsius]")
	title!(plt[1], "[C1,C4,C10] ~ $z")

	contour!(plt[2], p./ubar, T .- 273.15, data2, levels = 10, fill=(true,cgrad(:jet)))
	ylabel!(plt[2], "Pressão [Bar]")
	xlabel!(plt[2], "Temperatura [°Celsius]")
	title!(plt[2], "[C1,C4,C10] ~ $z2")
end

# ╔═╡ dfb527ab-d7d8-4e5f-91fb-4c8055884a64
md"""
# Reservatório de gás com influxo

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/gas-influx.png?raw=true)

Do Balanço de material (MBAL):

$$Gas\ produzido = Gas\ inicial(p_i)\ -\ Gas\ remanescente(p)$$

Portanto:

$$G_p = G - \left[\frac{V_r\phi N_G S_{gi}}{B_{g}} - \frac{W_e}{B_g}\right] = G\left(1-\frac{B_{gi}}{B_g}\right) + \frac{W_e}{B_g}$$

$$G_pB_g = G\left(B_g-B_{gi}\right) + W_e$$

* O mecanismo de influxo de água adiciona um "suporte para a pressão", porém pode aprisionar o gás
* Os gráficos $\frac{p}{Z}\times G_p$ é não-linear
* **NÃO DEVE** extrapolar os dados para encontrar G, mesmo se a primeira vista os dados apresentarem comportamento linear (superestimação)

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-lin-influx.png?raw=true)

## Diagnosticando o influxo de água

![h:350 center](https://ars.els-cdn.com/content/image/3-s2.0-B978012813649200013X-f13-09-9780128136492.jpg)

Para identificar se a produção do reservatório está sob influência de um aquífero, fazemos o uso do método de Havlena-Odeh, através de:

$$\frac{F}{E_g} = \frac{G_pB_g + W_pB_w}{B_g - B_{gi}} = G + \frac{W_e}{B_g - B_{gi}}$$

Considerando o modelo de VE (aquífero linear ou radial),

$$W_e = C\sum_{i=1}^{n-1}\Delta p_i W_D(t_{Dn} - t_{Dj})$$

Inserindo na MBAL e rearranjando a equação em termos de energia produzida x energia do reservatório:

$$\frac{G_pB_g}{B_g - B_{gi}} = G + C\frac{\sum_{i}\Delta p_i W_{Di}}{B_g - B_{gi}}$$

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/Havlena-Odeh.png?raw=true)

!!! info "Dica"
	Quando o método Havlena-Odeh **NÂO** apresentar valores de OGIP positivos, utilizar o plot $G_p\times G+\frac{W_e}{B_g - B_{gi}}$

	![h:350 center](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-lin-influx-Gp.png?raw=true)

## Fator de recuperação máximo
![h:350 center](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/gas-influx.png?raw=true)

Máxima recuperação: $S_{gi}\rightarrow S_{gr}$

Influxo teórico máximo:

$$W_e = V_r\phi N_G(S_{gi}-S_{gr})\rightarrow W_e^{max} = GB_{gi}\left(1 - \frac{S_{gr}}{S_{gi}}\right)$$

Substituindo na MBAL:

$$G_p = G\left[1 - \frac{B_{gi}}{B_g}\times\frac{S_{gr}}{S_{gi}}\right]$$
"""

# ╔═╡ 400d72eb-469e-4019-b980-d3082a3d65ba
details(
	md"""**Atividade de Revisão 2.** Calcule o influxo de água para o sistema reservatório-aquífero, usando os valores do histórico de pressão da Tabela abaixo. Suponha um aquífero infinito. Encontre o influxo de água em cada passo de tempo, usando van Everdingen & Hurst e compare o resultado com os modelos de Carter-Tracy, Fetkovich, Schilthuis.

	| tempo (dias) | pressão rreservatório (psia) |
	| ----------- | ------------------------ |
	| 0 | 3793 |
	| 91,5 | 3788 |
	| 183,0 | 3774 |
	| 274,0 | 3748 |
	| 366,0 | 3709 |
	| 457,5 | 3680 |
	| 549,0 | 3643 |
	
	As propriedades estimadas do aquífero são:
	
	$\phi=$ 0,209; $\mu$ =0,25 cp; $\theta$ =180°; $k$ =275 md; $c_t=6\times 10^{-6}$ psia; $h$ = 19,2 pés; $r_e$ =5807 pés;
	""",
	md"""
	Reutilizar e incrementar as funções contidas [neste notebook](https://colab.research.google.com/drive/1zX8PzyHAiAPZ7L6w1rotHmN2Rf5C1kvk?usp=sharing) ou [nesta planilha](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/template-report-projeto.xlsm?raw=true), para resolução da atividade.
	"""
)

# ╔═╡ dab79865-074e-4df6-8d63-001979fb881b
md"""
## Determinação dos parâmetros do aquífero através do MBAL

O tamanho e a constante de influxo do aquífero ($U$ ou $C$ a depender da referência) assim como a reserva *in place* podem ser determinadas (estimadas) através do balanço de material, para isso, devemos seguir os passos:

1. Definindo a forma da equação de balanço mais adequada:

$$G_pB_g = G\left(B_g-B_{gi}\right) + W_e$$

ou

$$\frac{p}{Z} = \frac{p_i}{Z_i} \frac{\left(1 - \frac{G_p}{G}\right)}{\left(1 - \frac{W_e}{GB_{gi}}\right)}$$

2. Escolhendo $W_e$, o modelo de aquífero (VEH, Carter-Tracy, Fetkovich, Schilthuis, etc.)

3. Estimar os parâmetros do aquífero e reserva por ajuste de histórico.

Considere o reservatório de gás natural cuja a temperatura da formação é de 222°F e a saturação de água conata de 0.26. Além disso, o mecanismo de produção da reserva é puramente expansão do gás com influxo de água, apresenta o seguinte histórico.

| Pressão Res. (psia) | Fator compress. (-) | Produção acumulada de gás (MMMscf) |
| ------------------- | ------------------- | ---------------------------------- |
|  4930.0             |     0.597712        |         0.0     |  
|  4785.0             |     0.738801        |        72.0     |  
|  4640.0             |     0.829404        |       112.0     |  
|  4495.0             |     0.824446        |       129.0     |  
|  4350.0             |     0.899273        |       155.0     |  

Seguindo os passos temos que:

1. Definindo a forma da equação de balanço mais adequada:

$$G_pB_g = G\left(B_g-B_{gi}\right) + W_e$$

ou

$$\frac{p}{Z} = \frac{p_i}{Z_i} \frac{\left(1 - \frac{G_p}{G}\right)}{\left(1 - \frac{W_e}{GB_{gi}}\right)}$$

Visualmente percebemos que a forma mais "simples" é em função de $B_g$, desta forma, criamos a columa $B_g$

| Pressão Res. (psia) | Fator compress. (-) | Produção acumulada de gás (MMMscf) | Bg (res cf/scf)   |
| ------------------- | ------------------- | ---------------------------------- | --------------|
|  4930.0 | 0.597712 |   0.0 | 0.00234 |
|  4785.0 | 0.738801 |  72.0 | 0.00298 |
|  4640.0 | 0.829404 | 112.0 | 0.00345 |
|  4495.0 | 0.824446 | 129.0 | 0.00354 |
|  4350.0 | 0.899273 | 155.0 | 0.00399 |

2. Escolhendo $W_e$, o modelo de aquífero (VEH, Carter-Tracy, Fetkovich, Schilthuis, etc.)

Para estimarmos a força do aqúifero, utilizaremos o modelo de aquífero mais simples (*Pot Model*):

$$W_e = c_tW_i\Delta p_i$$

Portanto,

$$G_pB_g = G\left(B_g-B_{gi}\right) + c_tW_i\Delta p_i$$

$$\frac{G_pB_g}{B_g - B_{gi}} = G + Wc_t\frac{\Delta p_i}{B_g - B_{gi}}$$

Adicionando mais uma coluna $\frac{\Delta p_i}{B_g - B_{gi}}$

| Pressão Res. (psia) | Fator compress. (-) | Produção acumulada de gás (MMMscf) | Bg (res cf/scf)   | $\frac{\Delta p}{B_g-B_{gi}}$ (psia.scf/res cf) |
| ------------------- | ------------------- | ---------------------------------- | --------------| -------- |
|  4930.0             |     0.597712        |         0.0     |    0.00234   |       |
|  4785.0             |     0.738801        |        72.0     |    0.00298   |   226562.50 |
|  4640.0             |     0.829404        |       112.0     |    0.00345   |   261261.26 |
|  4495.0             |     0.824446        |       129.0     |    0.00354   |   362499.99 |
|  4350.0             |     0.899273        |       155.0     |    0.00399   |   351515.15 |


"""

# ╔═╡ 0c1776cd-337b-4a59-a0da-f024abb5973d
let
	x  = [1e-6, 226562.50, 261261.26, 362499.99, 351515.15]
	Bg = [0.00234, 0.00298, 0.00345, 0.00354, 0.00399]
	Gp = [0.0, 72.0, 112.0, 129.0, 155.0]
	y = Gp .* Bg ./(Bg .- Bg[1])
	pred = range(minimum(x[2:end]),maximum(x[2:end]))
	# MBAL
	mbal(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	fit = curve_fit(mbal, x[2:end], y[2:end], p0)
	G = ceil(fit.param[1])
	ctW = round(fit.param[2], digits=6)
	# dominance
	aquifer = ctW*diff(x)
	
	
	# plotting
	scatter(x, y, label="Dados de Produção", xlabel="Δp/(Bg - Bgi)", ylabel="Gp, MMMscf", title="Reservatório de gás seco com influxo")
	plot!(pred,mbal(pred,fit.param),label="MBAL", c=:red, lw=2)
	plot!(range(0,minimum(x[2:end])),mbal(range(0,minimum(x[2:end])),fit.param),label="Forecast", c=:red, lw=2, ls=:dash)
	annotate!(0.85e5, 350, "OGIP = $G MMscf")
	annotate!(0.85e5, 360, "       cₜWᵢ = $ctW MMscf")
	xlims!(0, ceil(maximum(x[2:end])+5500))
end

# ╔═╡ 7b6b3529-3cfe-4909-9182-00fcc2a3e7cc
md"""
Considerando que a compressibilidade total ($c_t$) é $\approx 3\times 10^{-6}\ \text{psi}^{-1}$ podemos determinação dos parâmetros do aquífero através de:

$$W_i\approx \frac{0.000322\times 10^{6}}{3\times 10^{-6}}\approx 1.07334\times 10^8\ \text{scf}$$

A partir da geometria do aquífero (linear ou radial) podemos determinar o raio interno e externo.

Uma importante análise em problemas de balanço de material, é determinar qual mecanismo de produção é domínante (neste exemplo seria expansão do gás e influxo de água). Para isso, reescrevemos a MBAL da seguinte forma:

$$\underbrace{\frac{G(B_g - B_{gi})}{G_pB_g}}_{\text{expansão do GN}} + \underbrace{\frac{c_tW\Delta p_i}{G_pB_g}}_{\text{influxo de água}} = 1$$

A partir desta formula podemos identificar qual mecanismos de produção tem maior dominância (termo de maior valor) ao longo do ciclo de vido da reserva.
"""

# ╔═╡ a657ba78-dcf0-43ca-a83b-a8f3d36cd9da
let
	pres = [4930.0, 4785.0, 4640.0, 4495.0, 4350.0]
	Bg = [0.00234, 0.00298, 0.00345, 0.00354, 0.00399]
	Gp = [0.0, 72.0, 112.0, 129.0, 155.0]
	
	y = Gp .* Bg ./(Bg .- Bg[1])
	x  = (pres[1] .- pres) ./ (Bg .- Bg[1] .+ 1e-16)
	
	# MBAL
	mbal(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	fit = curve_fit(mbal, x[2:end], y[2:end], p0)
	G = ceil(fit.param[1])
	ctW = round(fit.param[2], digits=6)

	#dominace
	expan  = G .* (Bg .- Bg[1]) ./ (Gp .* Bg .+ 1e-16)
	aquifer = ctW .* (pres[1] .- pres) ./ (Gp .* Bg .+ 1e-16)
	
	
	# plotting
	ndat = length(expan)
	plot(pres, expan.+aquifer, fillrange=aquifer, fillalpha=0.25, c = :gray, label ="Expansão gás", legend = :topright, dpi = 100, xlabel="pressão [psia]", ylabel="Fração da força motriz")
	plot!(pres, aquifer, fillrange=zeros(ndat), fillalpha=0.25, c = :blue, label = "Influxo de aquífero", legend = :topright, dpi = 100)
	# ylims!(0,1)
end

# ╔═╡ 07d540f5-b791-42c4-8138-aa98554fcf04
md"""
## Predição do comportamento

Quando determinamos a reserva (OGIP) e as características do aquífero associado, a etapa subsequente é realizar a predição do comportamento. Este procedimento para este tipo de reservatório é iterativo, para os demais casos, utilizamos a forma linear do ajuste do histórico (análise sobre as linha tracejada).

No processo iterativo, conhecemos a reserva (OGIP) e caracteristicas do aquífero associado, portanto o interesse na predição é calculara a produção acumulada de gás natural, água e a pressão média do reservatório em um instante a *posteriori*. Deste modo, utilizamos a forma da MBAL em termos de pressão e fator de compressibilidade:

$$\frac{p}{Z} = \frac{p_i}{Z_i} \frac{\left(1 - \frac{G_p}{G}\right)}{\left(1 - \frac{W_e}{GB_{gi}}\right)}$$

A depender do modelo de influxo ($W_e$) utilizado o algoritmo de predição sobre alterações, por isso, apresentaremos o algoritmo para o modelo de VEH e para os modelos de Carter-Tracy e Fetkovich (são similares).

#### van Everdingen & Hurst

Neste algoritmo utilizamos o modelo de van Everdingen & Hurst (1949)[^3] para cálculo do influxo de água no instante $t_n$

$$W_e(t_{Dn})=U\sum_{i=1}^n\Delta p_i W_{D}(t_{Dn}-t_{Di-1})$$

O algoritmo para a previsão de comportamento utilizando este modelo consiste em:

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-prediction-influx-veh.png?raw=true)

#### Fetkovich

Para este algoritmo, o influxo acumulado[^3] é computado através de:

$$W_e = \frac{W_{ei}}{p_i}(p_i - p)\left[1 - e^{-\frac{jp_i}{W_{ei}}t}\right]$$

Desta forma, o influxo acumulado durante um intervalo de tempo $\Delta t_n$ é:

$$\Delta W_{en} = \frac{W_{ei}}{p_i}\left(\bar{p}_{a, n-1}-\bar{p}_n\right)\left[1 - e^{-\frac{jp_i}{W_{ei}}\Delta t_n}\right]$$

sendo $W_{en}$ o influxo acumulado até o instante $t_n$, $\bar{p}_{a, n-1}=p_i\left(1 - \frac{W_{e,n-1}}{W_{ei}}\right)$ pressão média do aquífero no instante $t_{n-1}$ e $\bar{p}_n=\frac{p_{n-1} + p_n}{2}$ é a média das pressões no contato[^4] durante o intervalo de tempo $\Delta t_n$, e $p$ é a pressão no contato.

Considerando intervalos constantes de $\Delta t_n$, de forma compacta temos:

$$\Delta W_{en}(p_n) = \alpha_n\left(\bar{p}_{a, n-1}-0.5[p_{n-1} + p_n]\right)$$

sendo, $$\alpha_n = \frac{W_{ei}}{p_i}\left[1 - e^{-\frac{jp_i}{W_{ei}}\Delta t_n}\right]$$

Inserindo na equação de balanço de material:

$$\frac{p_n}{Z_n} = \frac{p_i}{Z_i}\frac{\left(1 - \frac{G_{pn}}{G}\right)}{\left(1 - \frac{W_{e,n-1}+\Delta W_{en}(p_n)}{GB_{gi}}\right)}$$

Devido ao caráter implícito e não-linear da equação ($W_{en}(p_n)$ e $Z_n\rightarrow Z(p_n)$), o processo de solução deve ser iterativo, este consite em:

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-gas-prediction-influx-fet.png?raw=true)


!!! warn "Implementar algoritmos"
	Fazer a implementação dos algoritmos apresentados (VEH e Fetkovich) e adaptar esta último algoritmo para o modelo de Carter-Tracy

----

[^3] Visto na disciplina Engenharia de Reservatórios 1

[^4] por simplicidade admite-se que a pressão no contato seja igual à pressão (média) do reservatório.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
MultiComponentFlash = "35e5bd01-9722-4017-9deb-64a5d32478ff"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
HypertextLiteral = "~0.9.5"
LsqFit = "~0.15.0"
MultiComponentFlash = "~1.1.16"
Plots = "~1.40.9"
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.7"
manifest_format = "2.0"
project_hash = "fa863f2d0ab31335c395a344f68d79ecb661b412"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "0ba8f4c1f06707985ffb4804fdad1bf97b233897"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.41"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "f36e5e8fdffcb5646ea5da81495a5a7566005127"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7901a6117656e29fa2c74a58adb682f380922c47"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.116"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e51db81749b0777b2147fbe7b783ee79045b8e99"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "84e3a47db33be7248daa6274b287507dd6ff84e8"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.2"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "21fac3c77d7b5a9fc03b0ec503aa1a6392c34d2b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.15.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "786e968a8d2fb167f2e4880baba62e0e26bd8e4e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "846f7026a9decf3679419122b49f8a1fdb48d2d5"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.16+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "424c8f76017e39fdfcdbb5935a8e6742244959e8"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "b90934c8cb33920a8dc66736471dc3961b42ec9f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "c67b33b085f6e2faf8bf79a61962e7339a81129c"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.15"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "b1c2585431c382e3fe5805874bda6aea90a95de9"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.25"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "71b48d857e86bf7a1838c4736545699974ce79a2"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.9"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "10da5154188682e5c0726823c2b5125957ec3778"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.38"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "df37206100d39f79b3376afb6b9cee4970041c61"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.51.1+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89211ea35d9df5831fca5d33552c02bd33878419"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e888ad02ce716b319e6bdb985d2ef300e7089889"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "688d6d9e098109051ae33d126fcfc88c4ce4a021"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.1.0"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "40acc20cfb253cf061c1a2a2ea28de85235eeee1"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultiComponentFlash]]
deps = ["ForwardDiff", "LinearAlgebra", "Roots", "StaticArrays"]
git-tree-sha1 = "bf038b993aeaad2f99fee2690cf96a6a389a5e03"
uuid = "35e5bd01-9722-4017-9deb-64a5d32478ff"
version = "1.1.16"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "030ea22804ef91648f29b7ad3fc15fa49d0e6e71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+3"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ed6834e95bd326c52d5675b4181386dfbe885afb"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.55.5+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "b5e7c125687aa818de948ef8a0b9dc59f043342c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.1"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "f233e0a3de30a6eed170b8e1be0440f732fdf456"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.4"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "47091a0340a675c738b1304b58161f3b0839d454"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.10"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5db3e9d307d32baba7067b13fc7b5aa6edd4a19a"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.36.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "beef98d5aad604d9e7d60b2ece5181f7888e2fd6"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.4+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e9216fdcd8514b7072b43653874fd688e4c6c003"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.12+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "807c226eaf3651e7b2c468f687ac788291f9a89b"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.3+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89799ae67c17caa5b3b5a19b8469eeee474377db"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.5+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+3"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "6fcc21d5aea1a0b7cce6cab3e62246abd1949b86"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.0+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "984b313b049c89739075b8e2a94407076de17449"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.2+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a1a7eaf6c3b5b05cb903e35e8372049b107ac729"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.5+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "b6f664b7b2f6a39689d822a6300b14df4668f0f4"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.4+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a490c6212a0e90d2d55111ac956f7c4fa9c277a6"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+1"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c57201109a9e4c0585b208bb408bc41d205ac4e9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.2+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "dbc53e4cf7701c6c7047c51e17d6e64df55dca94"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+1"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "ab2221d309eda71020cdda67a973aa582aa85d69"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+1"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6dba04dbfb72ae3ebe5418ba33d087ba8aa8cb00"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.1+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "622cf78670d067c738667aaa96c553430b65e269"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0ba42241cb6809f1a278d0bcb976e0483c3f1f2d"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7b5bbf1efbafb5eca466700949625e07533aff2"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.45+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "63406453ed9b33a0df95d570816d5366c92b7809"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+2"
"""

# ╔═╡ Cell order:
# ╟─78d10c30-d382-11ef-2a86-1d3c01e82daa
# ╟─369a85ef-d489-4f2c-83aa-0348897aafda
# ╟─ce8cf761-58df-4769-af0b-b71efdf7d2f6
# ╟─785e1863-8001-4885-8829-9f98419d6cd9
# ╟─465c140f-4729-416c-b51e-901836d545e8
# ╟─5c641485-6bea-498e-874b-c0fc0ed99743
# ╟─856145bf-8853-4d07-81dd-fcecbc45c422
# ╟─7a2c0627-5a20-42f4-b99f-ef9a37b00aed
# ╟─6532b23d-0a83-4cac-82ce-166dfef2ffb1
# ╟─b07ac32c-9017-4103-8614-3d9203413adb
# ╟─a01e834a-c916-4599-bf89-08500e7790bb
# ╟─a0689d18-6906-42b7-b474-052819995933
# ╟─eb272cd0-2472-4b48-a35a-6868e39de7c5
# ╟─f52cd2f2-ec88-4ec3-a33c-f0f96909e776
# ╟─5319239b-f361-4375-b49d-596b36a89641
# ╟─2133ca27-a8ea-4e93-88b0-412d91fe5efb
# ╟─dfb527ab-d7d8-4e5f-91fb-4c8055884a64
# ╟─400d72eb-469e-4019-b980-d3082a3d65ba
# ╟─dab79865-074e-4df6-8d63-001979fb881b
# ╟─0c1776cd-337b-4a59-a0da-f024abb5973d
# ╟─7b6b3529-3cfe-4909-9182-00fcc2a3e7cc
# ╟─a657ba78-dcf0-43ca-a83b-a8f3d36cd9da
# ╟─07d540f5-b791-42c4-8138-aa98554fcf04
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
