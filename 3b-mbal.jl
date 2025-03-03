### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 78d10c30-d382-11ef-2a86-1d3c01e82daa
using Plots, PlutoUI, HypertextLiteral, LsqFit, MultiComponentFlash, Interpolations, NLsolve, Printf

# ╔═╡ 286c1e61-8d3c-4570-8611-985f4e88fe56
begin
	# ----------------------------< auxilary function >------------------------
	# specific gravity of oil
	γoil(a) = 141.5 / (131.5 + a)
	
	# solution gas-oil ratio of saturated oil (Rs scf/stb)
	Rs(a,yg,T,p, FO1) = yg * (10 ^ (0.0125 * a - 0.00091 * T) * p / (18 * FO1)) ^ (1 / 0.83)

	# bubble point pressure (Pb psia) from Frick.
	Pbubble(a,yg,T,gor,FO1) = 18 * FO1 * (gor / yg) ^ 0.83 * 10 ^ (0.00091 * T - 0.0125 * a)

	# oil formation volume factor (Bo rbl/stb)
	function Bo(yg, yo, rs, gor, Co, p, Pb, T, FO2, FO3)
		# yg, Gas specific gravity
		# yo, oil specific gravity
		# rs, solution gas-oil ratio of saturated oil, scf/stb
		# gor, Gas oil ratio in scf/stb
		# Co, Oil compressibility, 1/psi
		# P, Pressure, psia
		# Pb, bubble point pressure, psia
		# T, in °F
		# FO2, Tuning factor
		# FO3, Tuning factor
	    if (p <= Pb) 	# saturated oil
	        F = rs * (yg / yo) ^ 0.5 + 1.25 * T
	        FVF = 0.972 * FO2 + 0.000147 * FO3 * F ^ 1.175
		else 			# undersaturated oil
	        F = gor * (yg / yo) ^ 0.5 + 1.25 * T
	        Bob = 0.972 * FO2 + 0.000147 * FO3 * F ^ 1.175
	        FVF = Bob * exp(Co * (Pb - p))
		end
		FVF
	end

	# oil compressibility (Co 1/psi) from Vazquez and Beggs
	function Co(yg, yo, Rsb, API, p, T, FO4)
		#yg, Gas specific gravity
		#yo, oil specific gravity
		#Rsb, gas solubility at the bubble-point pressure in scf/stb
		#API, API of the oil
		#P, Pressure, psia
		#T, in °F
		#FO4, Tuning factor
		#Pb Bubble point pressure, psia
		ygS = yg * (1 + 0.00005912 * yo * 60 * log10(14.7 / 114.7))
		Co = FO4 * (-1433 + 5 * Rsb + 17.2 * T - (1180 * ygS) + 12.61 * API) / (p * 100000)
		return Co
	end

	
	# oil density [saturated and undersaturated] in kg/m³ from Beggs and Robinson.
	function ρₒ(yo, yg, Rs, GOR, Bo, p, Pb)
		#yg, Gas specific gravity
		#yo, oil specific gravity
		#Rs, solution gas-oil ratio of saturated oil, m^3/m^3
		#GOR, Gas oil ratio in m^3/m^3 (undersaturated oil)
		#Bo, oil formation volume factor, m^3/m^3
		#p pressure, psia
		#Pb bubblepoint pressure, psia
    	if (p < Pb)
        	ρ = (yo * 1000 + Rs * 1.223 * yg) / Bo         # [kg/m3]
    	else
        	ρ = (yo * 1000 + GOR * 1.223 * SGg) / Bo
    	end
	end

	function zFactor(pᵣ, Tᵣ)
		T = 1 / Tᵣ
		y = 0.001
		α = 0.06125 * T * exp(-1.2 * (1- T) ^ 2)

		f = 1
		while abs(f) > 1e-9
			f = -α * pᵣ + (y + y ^ 2 + y ^ 3 - y ^ 4) / (1 - y) ^ 3 - (14.76 * T - 9.76 * T ^ 2 + 4.58 * T ^ 3) * y ^ 2 + (90.7 * T - 242.2 * T ^ 2 + 42.4 * T ^ 3) * y ^ (2.18 + 2.82 * T)
			
			∂f = (1 + 4 * y + 4 * y ^ 2 - 4 * y ^ 3 + y ^ 4) / (1 - y) ^ 4 - (29.52 * T - 19.52 * T ^ 2 + 9.16 * T ^ 3) * y + (2.18 + 2.28 * T) * (90.7 * T - 242.2 * T ^ 2 + 42.4 * T ^ 3) * y ^ (1.18 + 2.82 * T)
			
			y = y - f / ∂f
		end
		return α * pᵣ / y
	end

	Ppc(x) = 677 + 15 * x - 37.5 * x ^ 2
	Tpc(x) = 168 + 325 * x - 12.5 * x ^ 2

	Bg(z, p, T) = z * (14.7 / p) * (T + 460) / (60 + 460) # in scf/stb 
end

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

# ╔═╡ c1cb7823-f336-4c66-b92e-de78bae8aad0
begin
	Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)
	@htl"""
	<button onclick="present()">Apresentar</button>
	<div style="margin-top:3em;margin-bottom:7em;">
	</div>
	<title>Engenharia de Reservatórios 2</title>
	<subtitle>MBAL para reservatórios de óleo</subtitle>
	<author>Jonathan da Cunha Teixeira</author>
	<email><a href="mailto:jonathan.teixeira@ctec.ufal.br">jonathan.teixeira@ctec.ufal.br<a/></email>
	<semester>Engenharia de Petróleo<br>Universidade Federal de Alagoas</semester>
	<hr style="border-top:8px dashed;margin:2em 0em;"/>
	"""
end

# ╔═╡ 5c641485-6bea-498e-874b-c0fc0ed99743
md"""
# Balanço de Material (*MBAL*) de reservatórios de petróleo

## Objetivo

Derivar as relações do balanço de material para líquidos levemente compressíveis (óleo) na presença de outras fases (gás e água).

# Reservatório de petróleo

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/reservatohrio.png?raw=true)

Partindo de um reservatório inicialmente preenchido por:

* Gás (*capa de gás*)
* Óleo
* Água (*aquífero*)

Portanto temos que, o volume de óleo ($N$) é:

$$N = V\phi N_G\frac{s_{oi}}{B_{oi}}$$

e o volume de gás na capa ($G$):

$$G = V_{gr}\phi N_G\frac{s_{gi}}{B_{gi}}$$


estas equações contabilizam as "massas" de HC, *medidos sob condições-padrão*, e por simplificação definimos a variável auxiliar $m$ como sendo a razão entre entre o HC gasosos e os líquidos:

$$m = \frac{V_{gr}\phi N_G s_{gi}}{V\phi N_Gs_{oi}}$$

sendo: $V_{gr}\phi N_G s_{gi} = G B_{gi}$ e $V\phi N_Gs_{oi} = NB_{oi}$, são a quantidade inicial de óleo e gás (em RB), daí:

$$m = \frac{GB_{gi}}{NB_{oi}}$$

A análise aquí consistirá em quantificar a adição e/ou subtração uma certa quantidade de volume por meio do processo de depletação que envolve:

1. Expansão do óleo + gás em solução, $E_o$
1. Expansão da capa de gás, $E_{gc}$
1. Expansão da água/aquífero, $E_w$
1. Expansão da água conata e compressão da rocha, $E_{wr}$

Portanto, de forma geral o balanço de materiais em um reservatório de petróleo é:

$$\text{Expansão do volume Reservatório }$$
$$||$$
$$\text{Volume dos fluidos produzidos }(N_p, G_p, W_p)\text{ e injetados}$$

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil.png?raw=true)

"""

# ╔═╡ 856145bf-8853-4d07-81dd-fcecbc45c422
md"""
# Equação geral MBAL

## Expansão do óleo + gás em solução

Na fase líquido teremos duas situações:

1. Abaixo da pressão de bolha:
2. Acima da pressão de bolha:

| Status |  Fase | Propriedade | Volume  |
| -------| ------| ------------| -------| 
| Inicialmente | óleo | $B_{oi}$ | $NB_{oi}$|
| Depletando...| óleo | $B_{o}$ | $NB_{o}$|
| Inicialmente | gás em solução | $R_{si}$ | $NR_{si}$|
| Depletando...| gás em solução | $R_{s}$ | $NR_{s}$|

Portanto,

$$NE_o = N \left[(Bo - B_{oi}) + (R_{si} - R_s)B_g\right]$$

## Expansão da capa de gás

Da definição de $m$, reescrevemos:

$$G = m\frac{B_{oi}}{B_{gi}}N$$

| Status |  Fase | Propriedade | Volume  |
| -------| ------| ------------| ------- |
| Inicialmente | gás | $B_{gi}$ | $GB_{gi}$|
| Depletando...| gás | $B_{g}$ | $GB_{g}$|

Portanto,

$$NE_g = G (B_g - B_{gi}) = mN B_{oi}\left(\frac{B_g}{B_{gi}} -1\right)$$

## Expansão da água/aquífero

Apenas é a variação do volume do influxo (quantidade de água que invade a zona oleosa)

$$W_e = \text{modelo de influxo}$$

Modelos de influxo de aquífero:

* Schilthuis
* van Everdingen e Hurst
* Fetkovitch
* Carter-Tracy

## Expansão da água conata e volume poroso

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil.png?raw=true)

$$V_{hc} = [V_{oi} + V_{gi}]\phi\left(1-s_{wi}\right) = V_{pi}\left(1-s_{wi}\right)$$

Como ocorre a expansão do volume poroso e da água conata teremos:

$$V_{hc} = V_{pi} - V_{wi}$$

Expansão = Variação de volume

$$\frac{\partial V_{hc}}{\partial p} = \frac{\partial V_{pi}}{\partial p} -\frac{\partial V_{wi}}{\partial p}$$

### Compressibilidade da água conata e volume poroso

Da definição da compressibilidade,

$$\frac{\partial V_{hc}}{\partial p} = c\phi V_{pi} + c_w V_{wi} = V_{pi}\cdot\left(c_\phi + c_w s_{wi}\right)$$

Daí a varação de volume devido a variação de pressão será:

$$\frac{\Delta V_{hc}}{\Delta p_i} = \frac{V_{hci}}{1-s_{wi}}\left(c_\phi + c_w s_{wi}\right)$$

Portanto a variação da água conata e do meio poroso será:

$$\Delta V_{hc} = (1+m)NB_{oi}\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p_i = (1+m)N E_{wr}$$

Finalmente,

$$E_{wr} = B_{oi}\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p_i$$

## Fluidos produzidos e injetados

$$N_p Bo + G_p B_g + W_p B_w - \underbrace{W_IB_{wI}-G_IB_{gI}}_{Injetados}$$

*OBS.:* Uma parte do gás que é produzido advém do óleo, **gás dissolvido**!!!

$$N_p Bo + \underbrace{G_p B_g}_{\underbrace{N_p R_p B_g}_{Gás\ produzido} - \underbrace{N_p R_s B_g}_{Gás\ dissolvido}} + W_p B_w$$

sendo: $R_p=\frac{G_p}{N_p}$ razão gás-óleo produzido

Portanto,

$$F = N_pB_o + N_p(R_p - R_s)B_g + W_p B_w - W_IB_{wI}-G_IB_{gI}$$

#### Equação geral

$$\underbrace{NE_o}_{Expansão\ do\ óleo} + \underbrace{NE_g}_{Expansão\ do\ gás} + \underbrace{W_e}_{Expansão\ do\\ aquífero} + \underbrace{(1+m)NE_{wr}}_{Expansão\ da\ água\\ conata\ e\ rocha} = F$$

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil.png?raw=true)

$N_p\left[B_o + (R_p - R_s)B_g\right] + W_p B_w - W_IB_{w}-G_IB_{g}=$

$N \left[(Bo - B_{oi}) + (R_{si} - R_s)B_g\right]$

$+$

$mN B_{oi}\left(\frac{B_g}{B_{gi}} -1\right)$

$+$

$(1+m)NB_{oi}\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p$

$+$

$W_e$

Quando realizando analise MBAL em reservatorios de petróleo, devemos ser cuidadosos em relação ao estado atual e passado dos fluidos presentes, isto é, se o fluido está subsaturado ou saturado:

1. Acima do ponto de bolha ($p_b$)
2. Abaixo do ponto de bolha ($p_b$)

Na situação 1. é conveniente considera a compressibilidade da água conata e volume poroso, enquanto que 2. isto não se faz necessário. 
"""

# ╔═╡ b4af5d4e-01d4-491a-9e9e-40eabbd02c22
md"""
# Reservatório com gás em solução

Quando produzindo neste mecanismo de produção, as caracteristicas marcantes (dominantes) neste reservatório são:

* Sem influxo de água (aquífero)
* Sem capa de gás

Além disso, ao longo do ciclo de vida da produção, este reservatório pode estar:

1. Acima do ponto de bolha ($p_b$)
2. Abaixo do ponto de bolha ($p_b$)

"""

# ╔═╡ 7a2c0627-5a20-42f4-b99f-ef9a37b00aed
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

# ╔═╡ f22d739d-a5b4-4543-8a75-7612cf6b44bb
md"""
### 1. Acima do ponto de bolha (óleo subsaturado)

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil-bubble.png?raw=true)

A partir das caracteristicas teremos:

$$NE_o + \cancel{m}NE_g + \cancel{W_e} + (1+\cancel{m})NE_{wr} = N_pB_o + N_p\cancel{(R_p - R_s)}B_g + \cancel{W_p B_w}$$

Portanto,

$$N\underbrace{(E_o + E_{wr})}_{\text{Variação efetiva da zona de óleo}} = N_pB_o$$

Avaliando a expansão da fase óleo:

$$E_o =(Bo - B_{oi}) + \underbrace{\cancel{(R_{si} - R_s)}}_{\text{óleo subsaturado}}B_g$$

Somando as expansões da fase óleo, água conata e rocha (multiplicando e dividindo por $B_{oi}$ e depois por $s_{oi}$) e rearranjando temos:

$$E_{wr} + E_o =B_{oi}\left(\frac{c_\phi + c_w s_{wi} + c_o s_{oi}}{1-s_{wi}}\right)\Delta p = B_{oi} c_{eo}\Delta p$$

Finalmente,

$$N_pB_o = NB_{oi}c_{eo}\Delta p = N(E_o + E_{wr})$$

Em termos de declínio da pressão (previsão do comportamento):

$$p = p_i - \frac{B_o}{NB_{oi}c_{eo}} N_p$$

ou,

$$\underbrace{N_pB_o}_{F} = NB_{oi}c_{eo}\Delta p = N(E_o + E_{wr})$$

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil-bubble-linfit.png?raw=true)
"""

# ╔═╡ e9946302-c337-402d-984e-ab14e3932469
let
	API = 45
	γg = 1.026
	T = 150
	pb = 1775
	pᵢ = 3000
	cᵩ = 1e-6
	cₐ = 3e-6
	Bₒᵢ = 1.4987
	Rₛᵢ = 845
	co = 2.665e-5
	N = 2492514
	# If experimental data, fill bubble point values (at the end) and tune using this parameters
	TF01 = 1.0 # Tuning factor 1 
	TF02 = 1.0 # Tuning factor 2
	TF03 = 1.0 # Tuning factor 3
	TF04 = 1.0 # Tuning factor 4
	
	# ----------------------------< auxilary function >------------------------
	# specific gravity of oil
	γoil(a) = 141.5 / (131.5 + a)
	
	# solution gas-oil ratio of saturated oil (Rs scf/stb)
	Rs(a,yg,T,p, FO1) = yg * (10 ^ (0.0125 * a - 0.00091 * T) * p / (18.2 * FO1)) ^ (1 / 0.83)

	# bubble point pressure (Pb psia) from Frick.
	Pbubble(a,yg,T,gor,FO1) = 18 * FO1 * (gor / yg) ^ 0.83 * 10 ^ (0.00091 * T - 0.0125 * a)

	# oil formation volume factor (Bo rbl/stb)
	function Bo(yg, yo, rs, gor, Co, p, Pb, T, FO2, FO3)
		# yg, Gas specific gravity
		# yo, oil specific gravity
		# rs, solution gas-oil ratio of saturated oil, scf/stb
		# gor, Gas oil ratio in scf/stb
		# Co, Oil compressibility, 1/psi
		# P, Pressure, psia
		# Pb, bubble point pressure, psia
		# T, in °F
		# FO2, Tuning factor
		# FO3, Tuning factor
	    if (p <= Pb) 	# saturated oil
	        F = rs * (yg / yo) ^ 0.5 + 1.25 * T
	        FVF = 0.972 * FO2 + 0.000147 * FO3 * F ^ 1.175
		else 			# undersaturated oil
	        F = gor * (yg / yo) ^ 0.5 + 1.25 * T
	        Bob = 0.972 * FO2 + 0.000147 * FO3 * F ^ 1.175
	        FVF = Bob * exp(Co * (Pb - p))
		end
		FVF
	end

	
	# oil density [saturated and undersaturated] in kg/m³ from Beggs and Robinson.
	function ρOil(yo, yg, Rs, GOR, Bo, p, Pb)
		#yg, Gas specific gravity
		#yo, oil specific gravity
		#Rs, solution gas-oil ratio of saturated oil, m^3/m^3
		#GOR, Gas oil ratio in m^3/m^3 (undersaturated oil)
		#Bo, oil formation volume factor, m^3/m^3
		#p pressure, psia
		#Pb bubblepoint pressure, psia
    	if (p < Pb)
        	ρ = (yo * 1000 + Rs * 1.223 * yg) / Bo         # [kg/m3]
    	else
        	ρ = (yo * 1000 + GOR * 1.223 * SGg) / Bo
    	end
	end

	# ------------------------< plotting >---------------------------
	p = range(2100,pᵢ,5)
	n = length(p)
	# solution gas-oil ratio plot
	solub = zeros(n)
	oilFVF = zeros(n)
	GOR = zeros(n)
	yo = γoil(API)
	Rp = 0.
	for i=1:n
		if p[i] > pb
			solub[i] = Rs(API,γg,T,pb,TF01)
		else
			solub[i] = Rs(API,γg,T,p[i],TF01)
		end
		oilFVF[i] = Bo(γg, yo, solub[i], solub[i], co, p[i], pb, T, TF02, TF04)
		GOR[i] = Rp * solub[i] * (Rp * 0.01 + 1.)
		Rp = 1
	end
	# cum. prod
	Np = N .* (oilFVF .- oilFVF[end]) ./ oilFVF
	# MBAL
	F = Np .* oilFVF + 2000 .* randn(n)
	# fix F[end] to zero (MANDATORY!)
	F[end] = .0
	Eo = (oilFVF .- oilFVF[end])	
	
	# fitting
	mbal(x, p) = x .* p[1]
	p0 = [0.5]
	fit = curve_fit(mbal, Eo, F, p0)
	scatter(Eo, F, label="Dados de Produção", xlabel="Eₒ = Bₒ - Bₒᵢ [rb/stb]", ylabel="F = Nₚ x Bₒ [stb]")
	plot!(Eo,mbal(Eo,fit.param),label="MBAL", c=:red, lw=2)
	
	
	# # prediction
	pred = range(pb,1800,5)
	solub .= 0.
	oilFVF .= 0.
	for i=1:n
		solub[i] = Rs(API,γg,T,pb,TF01)
		oilFVF[i] = Bo(γg, yo, solub[i], solub[i], co, pred[i], pb, T, TF02, TF04)
	end
	# MBAL
	Eopred = range(Eo[1],.05,5)
	plot!(Eopred,mbal(Eopred,fit.param),label="Forecast", c=:red, lw=2, ls=:dash)

	# M = zeros(n, 6) # Np, GOR, Bo, Rs, Eo, F
	# pr = [3000, 2874, 2600, 2376, 2100]
	# Rp = 0
	# rn = randn(n)
	# for i=1:n
	# 	M[i,4] = Rₛᵢ
	# 	M[i,3] = Bo(γg, yo, solub[i], solub[i], co, pr[i], pb, T, TF02, TF04)
	# 	M[i,2] = Rp * M[i,4] * (Rp * 1.01)
	# 	M[i,1] = Rp * (N * (M[i,3] - M[1,3]) + 2000 * abs(rn[i])) / M[i,3]
	# 	M[i,5] = M[i,3] - M[1,3]
	# 	M[i,6] = M[i,1] * M[i,3]
	# 	Rp = 1
	# end
	# M
	# # scatter(M[:,5], M[:,6], label="Dados de Produção", xlabel="Eₒ = Bₒ - Bₒᵢ [rb/stb]", ylabel="F = Nₚ x Bₒ [stb]")
end

# ╔═╡ 6532b23d-0a83-4cac-82ce-166dfef2ffb1
details(
	md"""**Exercício 1.** Um reservatório de petróleo apresenta um histórico de produção apresentado na tabela abaixo. Este reservatório está subsaturado e estima-se que não há aquífero associado e compressibilidade da água de formação e rocha são depreziveis. Estime o *Original Oil in place* (OOIP, as vezes chamado de *stock tank oil initially in place* STOIIP) e a reserva atual.

	| Pressão Res. (psia) | Produção acumulada (stb) | Razão gás-óleo (scf/stb) |
	| ------------------- | ------------------- | ----------------------------- |
	|  2940               |     8473.2          |       845     |  
	|  2600               |     26828.3         |       853     |  
	|  2376               |     41796.2         |       862     |  
	|  2100               |     59921.4         |       871     |  
	
	Dados adicionais (PVT do óleo): °API = 45, γg = 1.026, T = 150 °F, Pb = 1775 psia, Pᵢ = 3000 psia, cᵩ = 1e-6 psia$^{-1}$, cₐ = 3e-6 psia$^{-1}$, $s_{wi}=0.22$, e
	cₒ = 2.665e-5 psia$^{-1}$, Rₛᵢ = 845 scf/stb
	""",
let
	md"""
	Aplicando as considerações:

	1. Óleo está subsaturado
	1. Não há aquífero associado
	1. Compressibilidade da água de formação e rocha são depreziveis

	Desta forma, a equação de balanço fica:

	$$\underbrace{N_p\left[B_o + (R_p - Rs)B_g\right]}_{F} = N\underbrace{\left[(B_o - B_{oi})\right]}_{E_o}$$

	Portando precisamos realizar alguns cálculos PVTs (determinar $B_o, Z,$ e $B_g$) para poder realizar o ajuste de histórico. Portanto construimos a tabela:

	| P (psia) | Np (stb) | Rp (scf/stb) | Bo (rb/stb) | Z  | Bg (ft³/scf) |
	| -------- | -------- | ------------ | ----------- | --- | ------------ |
	|  2940    | 8473.2   |      845     | 1.534   | 0.885684 | 0.00519488 |
	|  2600    | 26828.3  |      853     | 1.54796 | 0.868748 | 0.00576188 |
	|  2376    | 41796.2  |      862     | 1.55723 | 0.854483 | 0.00620156 |
	|  2100    | 59921.4  |      871     | 1.56872 | 0.831828 | 0.00683059 |



	Com os dados do histórico de produção, realizamos o ajuste conforme o balanço de material através do gráfico F x Eo, onde encontramos OOIP através da inclinação da reta (forçando a **intercepção em Zero** ver figura a seguir)
	
	$$N = OOIP = STOIIP \approx 2.57813\ MMstb$$

	A reserva atual é de $\approx$ 2.5182 MMstb
	"""
end
)

# ╔═╡ e6ed5d77-eba5-422a-93a6-bdb71391d939
let
	API = 41
	pb = 1775 # psia
	pᵢ = 3000 # psia
	cᵩ = 1e-6 #1/psia
	cₐ = 3e-6 #1/psia
	swi = 0.22
	cₒ = 2.665e-5 #1/psia
	Rₛᵢ = 845 #scf/stb
	press = [3000, 2940, 2600, 2376, 2100]
	T = 150 # °F
	Nₚ = [0, 8473.2, 26828.3, 41796.2, 59921.4]
	Rₚ = [0, 845, 853, 862, 871]
	# param: API, T, pb
	p = [API, T, pb]
	fobj!(γ) = Rₛᵢ - Rs(p[1], γ[1], p[2], p[3], 1.0) 
	res = nlsolve(fobj!,[0.5])
	γg = res.zero[1] # should give out [1;5] as the solut
	
	# compute PVT properties	
	z = zFactor.(Ppc(γg) ./ press, Tpc(γg) / (T+460))
	bg = Bg.(z, press, T) 
	
	n = length(press)
	bo = zeros(n)
	for i=1:n
		bo[i] = Bo(γg, γoil(API), Rₛᵢ, Rₛᵢ, cₒ, press[i], pb, T, 1, 1)
	end
	bo

	Eo = bo .- bo[1]
	F = Nₚ .* (bo .+ (Rₚ .- Rₛᵢ).* bg ./ 5.615)
	mbal(x, p) = p[1] .* x
	p0 = [0.5]
	fit = curve_fit(mbal, Eo, F, p0)
	scatter(Eo, F, label="Dados de Produção", xlabel="Eₒ, rb/stb", ylabel="F, rb")
	plot!(Eo,mbal(Eo,fit.param),label="MBAL", c=:red, lw=2)
	
end

# ╔═╡ d81ddba8-30c4-458e-a1ff-d2e30d03f402
md"""

!!! warn "Atividade de Fixação"
	Repetir o procedimento do exercício anterior utilizando o gráfico Np x p. Além disso, comparar com o MBAL considerando a compressibilidade da água de formação e rocha, qual o erro quando deprezamos esta parcela? Podemos fazer isso?
"""

# ╔═╡ b07ac32c-9017-4103-8614-3d9203413adb
md"""
### 2. Abaixo do ponto de bolha (óleo saturado)

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil-abovebubble.png?raw=true)

Nesta fase, termos a ocorrência de:

1. gás liberado do óleo (gás em solução);
1. desprezar a compressibilidade da formação e da água residual;

A partir das suposições anteriores, teremos:

$$NE_o + \cancel{m}NE_g + \cancel{W_e} + \cancel{(1+m)NE_{wr}} = F = N_pB_o + N_p(R_p - R_s)B_g +$$
$$+\cancel{W_p B_w}$$

Portanto,

$$NE_o = N_pB_o + N_p(R_p - R_s)B_g = N(Bo - B_{oi}) + N(R_{si} - R_s)B_g$$

$$N (B_o - B_{oi} + (R_{si} - R_s)B_g) = N_p\left[B_o + (R_p - R_s)B_g)\right]$$

Observando a expressão final, temos que:

$$\frac{N_p}{N} = \frac{B_o - B_{oi} + (R_{si} - R_s)B_g}{B_o + (R_p - R_s)B_g}$$

Linearizando a expressão, a fim de utiliza-la para estimar a reserva original ($N$) é:

$$\underbrace{N_pB_o + N_p(R_p - R_s)B_g}_{F} = N\underbrace{(Bo - B_{oi}) + (R_{si} - R_s)B_g}_{E_o}$$

**Obs.#1**: Esta forma é utilizada como análise de diagnóstico do comportamento (produção). Após $N$ ser estimado, procede-se para o procedimento de "simulação"/predição através de métodos iterativos.

**Obs.#2**: Se durante o ciclo de vida da produção do reservatório, caso inicialmente subsaturado, realiza-se a análise de MBAL considerando a equação para óleo subsaturado até o ponto de bolha. A partir deste ponto, a análise é realizada via processos iterativos.


"""

# ╔═╡ daf6ba60-2c2b-4097-8283-3d661d1b247b
md"""
### Predição do comportamento

Quando abaixo do ponto de bolha, existem os métodos de predição:

1. Tarner (1944), onde a variável de controle é o $f_R = \frac{N_p}{N}$
1. Muskat (1949), onde a variável de controle é o $S_o$
1. Tracy (1955), onde a variável de controle é o RGO
1. Schilthuis (19XX), onde a variável de controle é a MBAL.

##### Método de Tarner (1944)

Reescrevemos a equação de balanço de material (MBAL) da seguinte forma:

$$\frac{G_p}{N_b} = \left(\frac{B_o}{B_g} - R_s\right)\cdot\left(1 - \frac{N_p}{N_b}\right) - \left(\frac{B_{ob}}{B_g} - R_{sb}\right)$$

sendo: $G_p$ Produção acumulada de gás **a partir do ponto de bolha**, $N_b$ Volume remanescente de óleo **a partir do ponto de bolha** (obtido na análise do reservatório acima da pressão de bolha), $B_{ob}$ fator volume-formação do óleo **a partir do ponto de bolha**, $R_{sb}$ razão de solubilidade **a partir do ponto de bolha**.

O método basea-se em **comparar a variação da produção de gás por volume de óleo existente definido pela MBAL e pelos dados de GOR**:

$$\left|\left(\frac{\Delta G_{ps}}{N_b}\right)_{MBAL} - \left(\frac{\Delta G_{ps}}{N_b}\right)_{GOR}\right| < Tolerância$$

A variação dos parâmetros do reservatório, no periodo de produção **a partir do ponto de bolha** de $t_b \Rightarrow t_j$, serão:

$$\left(\frac{\Delta G_{ps}}{N_b}\right)_{EBM} = \frac{G_{ps}|_{j+1} - G_{ps}|_{j}}{N_b}$$

$$\frac{\color{red}{ G_{ps}|_j}}{N_b} = \left(\color{red}{\frac{B_{oj}}{B_{gj}}} - \color{red}{R_{sj}}\right)\cdot\left(1 - \frac{\color{red}{ N_{pj}}}{N_b}\right) - \left(\frac{B_{ob}}{\color{red}{B_{gj}}} - R_{sb}\right)$$

Os termos $\color{red}{\text{em vermelho}}$ variam ao longo da produção do reservatório

Assim a expressão da variação fica:

$$\left(\frac{\Delta G_{ps}}{N_b}\right)_{EBM} = B_{ob}\cdot\left(\frac{1}{B_{gj}} - \frac{1}{B_{gj+1}}\right) + \left(\frac{B_{oj+1}}{B_{gj+1}}-R_{sj+1}\right)\left(1 - \frac{N_{pj+1}}{N_b}\right) - \left(\frac{B_{oj}}{B_{gj}}-R_{sj}\right)\left(1 - \frac{N_{pj}}{N_b}\right)$$

Para prever o comportamento do reservatório, utilizamos as seguintes variáveis auxiliares: razão gás/óleo instantânea ($R_j$) no j-ésimo tempo, valor médio da razão gás/óleo ($\bar{R}$), dados por:

$$R_{j} = \left(\frac{k_g}{k_o}\right)_j\left(\frac{\mu_o}{\mu_g}\right)_j\left(\frac{B_o}{B_g}\right)_j + R_{sj}$$ 

$$\bar{R} = 0.5\left(R_{j} + R_{j+1}\right)$$ 

Portanto, a variação da produção de gás por volume de óleo existente através dos dados da GOR será:

$$\left(\frac{\Delta G_{ps}}{N_b}\right)_{GOR} = \bar{R}\cdot\left(\frac{N_{psj+1}}{N_b} - \frac{N_{psj}}{N_b}\right)$$

Calculando e checando com a tolerância:

$$\left|\left(\frac{\Delta G_{ps}}{N_b}\right)_{MBAL} - \left(\frac{\Delta G_{ps}}{N_b}\right)_{GOR}\right| < Tolerância$$

##### Algorítmo

**Passo 1**: Escolhe um valor de pressão $p_{j+1}<p_j$ (chutado/estimado);

**Passo 2**: Determinar as propriedades dos fluidos (óleo e gás) para a pressão $p_{j+1}$;

**Passo 3**: Estimar o valor de $\frac{N_{psj+1}}{N_b}$.

**Passo 4**: Calcular o incremento da produção de gás por volume de óleo através da EBM $\left(\frac{\Delta G_{ps}}{N_b}\right)_{EBM}$.

**Passo 5**: Com o valor do fator de recuração do óleo, determinar as saturações de líquidos:

$$S_{Lj+1} = \left(1 - \frac{N_{psj+1}}{N_b}\right)\left(\frac{B_{oj+1}}{B_{ob}}\right)(1-s_{wb}) + s_{wb}$$

**Passo 6**: Calcular as razões gás/óleo instantâneas ($R_{j+1}$)

**Passo 7**: Calcular o valor médio RGO ($\bar{R}$).

**Passo 8**: Determinar o incremento da produção de gás por volume de óleo através da RGO $\left(\frac{\Delta G_{ps}}{N_b}\right)_{RGO}$.

**Passo 9**: Comparar as variações obtidas nos passos 4 e 8.

**Passo 10**: Se o erro for maior que a tolerância permitida estimar um novo valor de fator de recuperação ($\frac{N_{psj+1}}{N_b}$) e repetir o processo a partir do passo 3, caso contrário, ir para o passo 1.

#### Método de Muskat (1949)

Muskat desenvolveu um método para a previsão do desempenho do reservatório em qualquer estágio de depletação de pressão, expressando a **equação de balanço de material** para um reservatório em depleção na **forma diferencial**.

Os volumes iniciais e em qualquer instante (depletação) são dados por:

* Óleo inical: $N = \frac{V_rS_{oi}}{B_{oi}}$

* Óleo: $N_r=N - N_p = \frac{V_pS_{o}}{B_{o}}$

* Gás dissolvido: $G_{diss}=\frac{V_pS_{o}}{B_{o}}R_{so}$

* Gás livre: $G_{livre}=\frac{V_pS_{g}}{B_{g}}=\frac{V_p(1-S_o-S_w)}{B_{g}}$

Quantidade de gás total restante no reservatório em pés cúbicos padrão (scf) é a soma dos gases livres e dissolvidos dada como:

$$G_r=\frac{V_pS_{o}}{B_{o}}R_{so}+\frac{V_p(1-S_o-S_w)}{B_{g}}$$

Diferenciando o volume de gás e óleo remanescente com respeito a presão:

$$\frac{\partial G_r}{\partial P}=V_p\left[\frac{S_o}{B_o}\frac{\partial R_{so}}{\partial P}-\frac{R_{so}S_o}{B_o^2}\frac{\partial B_{o}}{\partial P}+\frac{R_{so}}{B_o}\frac{\partial S_{o}}{\partial P} - \frac{1}{B_g}\frac{\partial S_{o}}{\partial P}-\frac{1 - S_o - S_w}{B_g^2}\frac{\partial B_{g}}{\partial P}\right]$$

$$\frac{\partial N_r}{\partial P}=V_p\left[\frac{1}{B_o}\frac{\partial S_o}{\partial P}-\frac{S_o}{B_o^2}\frac{\partial B_o}{\partial P}\right]$$

###### Razão gás-óleo instanânea

Definição exata: $R=\frac{\frac{\partial G_r}{\partial P}}{\frac{\partial N_r}{\partial P}}$

e a partir da equação de balanço de material (lei de Darcy):

$$R=\frac{B_ok_{rg}\mu_o}{B_gk_{ro}\mu_g}+R_{so}$$

###### Saturação de óleo em função da pressão

Igualando as razões de gás-óleo ($R$), e inserindo as formas diferenciais de $G_r$ e $N_r$, após alguns algebrismos os termos obtermos:

$$\frac{\partial S_o}{\partial P}=\frac{\frac{S_oB_g}{B_o}\frac{\partial R_{so}}{\partial P}-\frac{1-S_o-S_w}{B_g}\frac{\partial B_g}{\partial P}+\frac{k_{rg}\mu_o}{k_{ro}\mu_g}\frac{S_o}{B_o}\frac{\partial B_o}{\partial P}}{1+\frac{k_{rg}\mu_o}{k_{ro}\mu_g}}$$

Resolvemos a EDO atraves de qualquer método numérico!
"""

# ╔═╡ 5e3cd049-335d-4384-a903-87659e52587a
md"""
!!! info "Revisão: EDO solvers (resolvedores)"
	Uma equação diferencial é uma relação entre uma função, $f(x)$, sua variável independente, x, e qualquer número de suas derivadas. *Uma equação diferencial ordinária ou EDO* é uma equação diferencial onde a variável independente e suas derivadas estão em uma dimensão. Portanto, assumimos que uma EDO pode ser escrita como:

	$$\frac{dS(t)}{dt} = F(t, S(t))$$
	
	é uma EDO de primeira ordem explicitamente definida, ou seja, $F$ é uma função que retorna a derivada, ou mudança, de um estado dado um tempo e valor de estado. Além disso, seja $t$ uma grade numérica do intervalo $[t_0,t_f]$ com espaçamento $\Delta t$. Sem perda de generalidade, assumimos que $t_0 = 0$ e que $t_f = n\times \Delta t$ para algum inteiro positivo, $n$. A aproximação linear de $S(t)$ em torno de $t_j$ em $t_{j+1}$ é:
	
	$$S(t_{j+1}) = S(t_{j})+ (t_{j+1} - t_J)\frac{dS(t)}{dt}$$
	
	que também pode ser escrito,
	
	$$S(t_{j+1}) = S(t_j) + \Delta t\times F(t_j,S(t_j))$$
	
	Esta fórmula é chamada de Fórmula de Euler Explícita (direta). Ela nos permite calcular uma aproximação para o estado em S(tj+1) dado o estado em S(tj). Neste método, usamos apenas o item de primeira ordem na série de Taylor para aproximar linearmente a próxima solução. Partindo de um dado valor inicial de $S_0 = S(t_0)$, podemos usar esta fórmula para integrar os estados até $S(t_f )$; estes valores de $S(t)$ são então uma aproximação para a solução da equação diferencial. A fórmula de Euler explícita é o método mais simples e intuitivo para resolver problemas de valor inicial, isto é notório pela simplicidade do pseudo-código apresentado:

	![](https://raw.githubusercontent.com/pranabendra/articles/master/Euler-method/images/eqn_new_1.png)

	sendo xᵢ$\rightarrow t_i$ e yᵢ$\rightarrow S(t_i)$

	Em qualquer estado $(t_j,S(t_j))$, usamos $F$ naquele estado para “apontar” linearmente em direção ao próximo estado e então se move naquela direção uma distância de $\Delta t = h$, como mostrado aqui:

	![image.png](data:image/png;base64,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)

"""

# ╔═╡ 44ff5e94-08a4-4df0-9361-5f1cd924080a
md"""


##### Forma discreta simplificada

[Craft et al. 1991](https://www.amazon.com.br/Applied-Petroleum-Reservoir-Engineering-Craft/dp/0130398845) simplificou a equação acima, agrupando os termos dependentes de propriedades PVT como funções de:

$$X(P)=\frac{B_g}{B_o}\frac{\partial R_{so}}{\partial P}$$
$$Y(P)=\frac{\mu_o}{B_o\mu_g}\frac{\partial B_{o}}{\partial P}$$
$$Z(P)=\frac{1}{B_g}\frac{\partial B_{g}}{\partial P}$$

Portanto, a equação diferencial na forma discreta, mais simples, fica:

$$\frac{\Delta S_o}{\Delta P}=\frac{S_oX(P)+\frac{k_{rg}}{k_{ro}}S_oY(P)-(1-S_o-S_w)Z(P)}{1+\frac{k_{rg}\mu_o}{k_{ro}\mu_g}}$$

sendo: $\Delta S_o=S_o|_{i-1}-S_o|_{i}$ e $\Delta P=P|_{i-1}-P|_i$

**Obs.:** Como os valores de $B_g$ são menores que 1 e muito próximos!, para uma melhor precisão, o cálculo de $\frac{\partial B_g}{\partial P}$ é obtido quando $\frac{1}{B_g}$ é plotado com a pressão (P). Assim teremos:

$$\frac{\partial \frac{1}{B_g}}{\partial P}=-\frac{1}{B_g^2}\frac{\partial B_g}{\partial P}$$

$$\frac{\partial B_g}{\partial P}=-B_g^2\frac{\partial \frac{1}{B_g}}{\partial P}$$

$$Z(P)=\frac{1}{B_g}\left[-B_g^2\frac{\partial \frac{1}{B_g}}{\partial P}\right]=-B_g\frac{\partial \frac{1}{B_g}}{\partial P}$$

##### Algorítmo

Preparando os dados:

**Passo 1**: Encontrar as relações $\frac{\partial R_{so}}{\partial P},\frac{\partial B_{o}}{\partial P},$ e $\frac{\partial \frac{1}{B_g}}{\partial P}$ (i.e. determinar a inclinação das retas, coeficiente linear).

**Passo 2**: Calcular os termos dependentes das propriedades PVT's (X(P), Y(P), e Z(P)).

**Passo 3**: Calcular $\frac{\Delta S_o}{\Delta P}$ usando a saturação de óleo (normalmente é a saturação de óleo inicial ou do passo anterior, i-1) para a queda de pressão incremental $\Delta P_{i-1}$,

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{i-1}=\frac{(S_o)_{i-1}X(P_{i-1})+\frac{k_{rg}}{k_{ro}}(S_o)_{i-1}Y(P_{i-1})-(1-(S_o)_{i-1}-S_{wi})Z(P_{i-1})}{1+\left[\frac{k_{rg}\mu_o}{k_{ro}\mu_g}\right]_{i-1}}$$

**Passo 4**: Determinar a saturação de óleo na queda de pressão média $\Delta P$,

$$(S_o)_{i} = (S_o)_{i-1} - [P_{i-1} - P_{i}]\left(\frac{\Delta S_o}{\Delta P}\right)_{i-1}$$

**Passo 5**: Usando $(S_o)_i$ e $P_i$, recalcular $\left(\frac{\Delta S_o}{\Delta P}\right)_{i}$,

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{i}=\frac{(S_o)_{i}X(P_{i})+\frac{k_{rg}}{k_{ro}}(S_o)_{i}Y(P_{i})-(1-(S_o)_{i}-S_{wi})Z(P_{i})}{1+\left[\frac{k_{rg}\mu_o}{k_{ro}\mu_g}\right]_{i}}$$

**Passo 6**: Calcular a média para $\left(\dfrac{\Delta S_o}{\Delta P}\right)_{m}$ a partir dos valores previamente computados,

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{m}=\frac{\left(\frac{\Delta S_o}{\Delta P}\right)_{i-1}+\left(\frac{\Delta S_o}{\Delta P}\right)_{i}}{2}$$

**Passo 7**: Usando $\left(\dfrac{\Delta S_o}{\Delta P}\right)_{m}$, calcular a saturção de óleo,

$$(S_o)_{i} = (S_o)_{i-1} - [P_{i-1} - P_{i}]\left(\frac{\Delta S_o}{\Delta P}\right)_{m}$$

**Passo 8**: Calcular: $S_g$ e $N_p$,

$$(S_g)_i = 1 - (S_o)_{i} - S_{wi}$$

$$N_p = N\left[1 - \left(\frac{S_o}{1-S_{wi}}\right)\frac{B_{oi}}{B_o}\right]$$

**Passo 9**: Calcular a $R$ na pressão $P_i$,

$$R_{i+} = \left[\frac{B_ok_{rg}\mu_o}{B_gk_{ro}\mu_g}\right]_i + (R_{so})_i$$

**Passo 10**: Calcular a razão gás-óleo média,

$$R_m = \frac{R_{i+}+R_i}{2}$$

**Passo 11**: Calcular a produção acumulada de gás ($G_p$),

$$G_p= R_mN_p$$

**Passo 12**: Repetir os passos 3-11 para todas as quedas de pressão de interesse
"""

# ╔═╡ 87c10293-85e0-4b30-8af6-45588fb3d92a
details(
	md"""
	**Exercício 2:** Dado um reservatório de óleo saturado localizado no campo de óleo de Amassoma no estado de Bayelsa sem capa de gás; cuja pressão inicial é de 3620 psia e a temperatura do reservatório de 220°F. A saturação de água inicial (conata) é de 0.195 e da análise volumétrica, o OOIIP foi estimado em 45 MMSTB. Não há influxo aquífero. Os dados PVT são fornecidos na tabela a seguir.
	
	| Pressão(psia) | Bo(bbl/STB)| Rs(scf/STB) |Bg (bbl/scf) | μo (cp) | μg (cp) |
	| -----| ------ | ---- | -------- | ------ | ------ |
	| 3620 | 1.5235 |  858 | 0.001091 | 0.7564 | 0.0239 |
	| 3335 | 1.4879 |  796 | 0.001202 | 0.8355 | 0.0233 |
	| 3045 | 1.4533 |  734 | 0.001332 | 0.9223 | 0.0227 |
	| 2755 | 1.4187 |  672 | 0.001499 | 1.0199 | 0.0222 |
	| 2465 | 1.3841 |  610 | 0.0017   | 1.1253 | 0.0216 |
	| 2175 | 1.3496 |  549 | 0.001961 | 1.2431 | 0.0211 |
	| 1885 | 1.314  |  487 | 0.002296 | 1.3749 | 0.0205 |
	| 1595 | 1.2794 |  425 | 0.002762 | 1.5206 | 0.0199 |

	Neste campo, não há dados de permeabilidade relativa disponíveis. Assim, a correlação abaixo é usada para gerar a curva de permeabilidade relativa.
	
	$$\frac{k_{rg}}{k_{ro}} = 0.000149\cdot e^{12.57\times S_g}$$
	
	Calcular a produção acuulada de óleo e gás quando a pressão de 3335, 3045 e 2755 psia. Utilizando os métodos de Tarner e **Muskat**.
	""",
md"""
Vamos resolver utilizando o algoritmo de Muskat, seguindo o passo-a-passo detalhando anteriormente. Realizaremos o passo a passo para a pressão 3335 para as demais faremos uma tabelas mostrando os valores dos principais passos. Portanto, para $P=3335$ psia

**Passo 1**: Encontrar as relações $\frac{\partial R_{so}}{\partial P},\frac{\partial B_{o}}{\partial P},$ e $\frac{\partial \frac{1}{B_g}}{\partial P}$ (i.e. determinar a inclinação das retas, coeficiente linear).

$$\frac{\partial R_{so}}{\partial P} =0.2134828 \qquad\frac{\partial B_{o}}{\partial P}= 0.0001202\qquad\frac{\partial \frac{1}{B_g}}{\partial P}=0.2739908$$

**Passo 2**: Calcular os termos dependentes das propriedades PVT's (X(P), Y(P), e Z(P)) para $P_{i-1}$ e $P_{i}$, i.e $i-1=$ 3620 e $i=$ 3335.

$$X(3335)=\frac{Bg(3335)}{Bo(3335)}\frac{\partial R_{so}}{\partial P}\approx 0.0001724$$
$$Y(3335)=\frac{\mu_o(3335)}{Bo(3335)\mu_g(3335)}\frac{\partial B_{o}}{\partial P}\approx 0.0028673$$
$$Z(3335) = -Bg(3335)\frac{\partial \frac{1}{B_g}}{\partial P}\approx -0.000329$$

e

$$X(3620)=\frac{Bg(3620)}{Bo(3620)}\frac{\partial R_{so}}{\partial P}\approx 0.0001541$$
$$Y(3620)=\frac{\mu_o(3620)}{Bo(3620)\mu_g(3620)}\frac{\partial B_{o}}{\partial P}\approx 0.0023844$$
$$Z(3620) = -Bg(3620)\frac{\partial \frac{1}{B_g}}{\partial P}\approx -0.00030119$$

**Passo 3**: Calcular $\frac{\Delta S_o}{\Delta P}$ usando a saturação de óleo (normalmente é a saturação de óleo inicial ou do passo anterior, i-1) para a queda de pressão incremental $\Delta P_{i-1}$,

$$(S_o)_{i-1} = 1 - S_{wi} = 1 - 0.195\approx 0.805$$

Como não há gás livre inicialmente no reservatório $S_g = 0$

$$\frac{k_{rg}}{k_{ro}}= 0$$

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{i-1}=\frac{(S_o)_{i-1}X(P_{i-1})+\frac{k_{rg}}{k_{ro}}(S_o)_{i-1}Y(P_{i-1})-(1-(S_o)_{i-1}-S_{wi})Z(P_{i-1})}{1+\left[\frac{k_{rg}\mu_o}{k_{ro}\mu_g}\right]_{i-1}}$$

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{3620}=\frac{0.805\cdot 0.0001541 + 0 - 0}{1+ 0} = 0.00012406$$

**Passo 4**: Determinar a saturação de óleo na queda de pressão média $\Delta P$,

$$(S_o)_{i} = (S_o)_{i-1} - [P_{i-1} - P_{i}]\left(\frac{\Delta S_o}{\Delta P}\right)_{i-1}$$

$$(S_o)_{3335}\approx 0.769642$$

**Passo 5**: Usando $(S_o)_i$ e $P_i$, recalcular $\left(\frac{\Delta S_o}{\Delta P}\right)_{i}$,

$$(S_g)_{3335} = 1 - (S_o)_{3335} -  S_{wi}\approx 0.0353572$$

$$\frac{k_{rg}}{k_{ro}}= 0.000149\cdot e^{12.57\times S_g}\approx 0.00023238$$

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{3335}=\frac{0.769\cdot 0.00017 + 0.00023\cdot 0.769\cdot 0.00286 - 0.0353\cdot( -0.00033)}{1+ \frac{0.00023\cdot 0.8355}{0.0233}}$$

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{3335}\approx 0.0001437$$

**Passo 6**: Calcular a média para $\left(\dfrac{\Delta S_o}{\Delta P}\right)_{m}$ a partir dos valores previamente computados,

$$\left(\frac{\Delta S_o}{\Delta P}\right)_{m}=\frac{0.00012406+0.0001437}{2} \approx 0.0001338$$

**Passo 7**: Usando $\left(\dfrac{\Delta S_o}{\Delta P}\right)_{m}$, calcular a saturção de óleo,

$$(S_o)_{3335} = (S_o)_{3620} - [P_{3620} - P_{3335}]\left(\frac{\Delta S_o}{\Delta P}\right)_{m}$$

$$(S_o)_{3335}\approx 0.7668424$$

**Passo 8**: Calcular: $S_g$ e $N_p$,

$$(S_g)_i = 1 - (S_o)_{i} - S_{wi} = \approx 0.038157$$

$$N_p = N\left[1 - \left(\frac{S_o}{1-S_{wi}}\right)\frac{B_{oi}}{B_o}\right] = 45\times 10^6\left[1 - \left(\frac{0.7668424}{1-0.195}\right)\frac{1.5235}{1.4879}\right]\approx 1.107\times 10^6$$

$$(N_p)_{3335}\approx 1.107\times 10^6$$

**Passo 9**: Calcular a $R$ na pressão $P_i$,

$$\frac{k_{rg}}{k_{ro}}= 0.000149\cdot e^{12.57\times S_g}\approx 0.0002407$$

$$R_{i+} = \left[\frac{B_ok_{rg}\mu_o}{B_gk_{ro}\mu_g}\right]_{3335} + (R_{so})_{3335}$$

$$R_{i+} = \left[\frac{1.4879\cdot 0.0002407\cdot 0.8355}{0.001202\cdot 0.0233 }\right] + 796\approx 806.58$$

**Passo 10**: Calcular a razão gás-óleo média,

$$R_m = \frac{R_{i+}+Rs_i}{2}= \frac{806.58+796}{2}=801.29$$

**Passo 11**: Calcular a produção acumulada de gás ($G_p$),

$$G_p= R_mN_p\Rightarrow 801.29\times 1.107\times 10^6= 8.87332\times 10^8$$

**Passo 12**: Repetir os passos 3-11 para todas as quedas de pressão de interesse

Ao final teremos a seguinte tabela de valores

| Pressão(psia) | Sₒ (-) | Nₚ(MMstb) | Rₚ (scf/STB) | Gₚ (MMscf) |
| -----| ------ | ---- | -------- | ------ |
| 3620 | 0.805  | 0.0 | 858 | 0.0 |
| 3335 | 0.766842 | 1.10738 | 801.29 | 887.33 |
| 3045 | 0.721459 | 2.72192 | 743.564 | 2023.92 |
| 2755 | 0.668024 | 4.89851 | 690.585 | 3382.84 |

!!! info "Atividade de fixação"
	Terminar o exercício para as outras pressão (fazer à 🤚) e repetir com o método de Muskat, porém, agora resolvendo numericamente a EDO (método de Euler). Além disso, refazer o exercício, agora utilizando o método de Tarner.
"""
)

# ╔═╡ 8caa9071-4f3b-4d79-8164-6e3e4aee580b
let
	swi = 0.195
	N = 45e6
	p = [3620, 3335 , 3045, 2755, 2465, 2175, 1885, 1595]
	bo = [1.5235, 1.4879, 1.4533, 1.4187, 1.3841, 1.3496, 1.314, 1.2794]
	rs = [858, 796, 734, 672, 610, 549, 487, 425]
	bg = [0.001091, 0.001202,0.001332, 0.001499, 0.0017, 0.001961, 0.002296, 0.002762]
	μo = [0.7564, 0.8355, 0.9223, 1.0199, 1.1253, 1.2431, 1.3749, 1.5206]
	μg = [0.0239, 0.0233, 0.0227, 0.0222, 0.0216, 0.0211, 0.0205, 0.0199]
	kgo(x) = 0.000149 * exp(12.57*x)
	sg = zeros(length(p))
	so = zeros(length(p))
	Nₚ = zeros(length(p))
	Gₚ = zeros(length(p))
	Rₚ = zeros(length(p))
	
	so[1] = 1 - swi 
	
	# passo#1
	pvt(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	pBo = curve_fit(pvt, p, bo, p0)
	∂Bo = pBo.param[2]
	pRs = curve_fit(pvt, p, rs, p0)
	∂Rs = pRs.param[2]
	p1_Bg = curve_fit(pvt, p, 1 ./bg, p0)
	∂1_bg = p1_Bg.param[2]
	pμo = curve_fit(pvt, p, μo, p0)
	pμg = curve_fit(pvt, p, μg, p0)

	for i=2:length(p)
		p_ = p[i-1]
		pp = p[i]
		# passo#2
		x = 1 / pvt(pp, p1_Bg.param) / pvt(pp, pBo.param) * ∂Rs
		y = pvt(pp, pμo.param) / (pvt(pp, pBo.param) * pvt(pp, pμg.param)) * ∂Bo
		z = - 1/pvt(pp, p1_Bg.param) * ∂1_bg
		x_ = 1 / pvt(p_, p1_Bg.param) / pvt(p_, pBo.param) * ∂Rs
		y_ = pvt(p_, pμo.param) / (pvt(p_, pBo.param) * pvt(p_, pμg.param)) * ∂Bo
		z_ = - 1/pvt(p_, p1_Bg.param) * ∂1_bg
		# passo#3
		if i==2
			kg_ko = 0.
		else
			kg_ko = kgo(1 - so[i-1] - swi)
		end
		ΔSoi_1 = (so[i-1] * x_ + kg_ko * so[i-1] * y_ - (1 - so[i-1] - swi) * z_) / (1 + kg_ko * pvt(p_, pμo.param) / pvt(p_, pμg.param))
		# passo#4
		so[i] = so[i-1] - (p_ - pp) * ΔSoi_1
		# passo#5
		kg_ko = kgo(1 - so[i] - swi)
		ΔSoi = (so[i] * x + kg_ko * so[i] * y - (1 - so[i] - swi) * z) / (1 + kg_ko * pvt(pp, pμo.param) / pvt(pp, pμg.param))
		# passo#6
		ΔSo = 0.5 * (ΔSoi + ΔSoi_1)
		# passo#7
		so[i] = so[i-1] - (p_ - pp) * ΔSo
		# passo#8
		Nₚ[i] = N * (1 - so[i] / (1 - swi) * bo[1] / bo[i])
		# passo#9
		kg_ko = kgo(1 - so[i] - swi)
		R = kg_ko * bo[i] / bg[i] * pvt(pp, pμo.param) / pvt(pp, pμg.param) + rs[i]
		Rₚ[i] = 0.5 * (R + rs[i])
		Gₚ[i] =	Rₚ[i] * Nₚ[i]
	end
	plot(p, Nₚ./1e9, label="Nₚ", color=:red, lw=2, ylabel="Bscf| MMMstb")
	plot!(p, Gₚ./1e12, label="Gₚ", color=:blue, lw=2, xlabel="pressão, psia", ylabel="Bscf| MMMstb", legend=:right)
	plot!(twinx(),p, Rₚ, label="Rₚ", color=:black, lw=2, ylabel="scf/stb")
end

# ╔═╡ 7f033f77-1e6a-4954-bd43-ed09a10c8ade
md"""
# Reservatório com capa de gás

Para reservatórios produzindo sob efeito de uma capa, obtermos o seguinte MBAL:

$$\overbrace{N_p\left[B_o + (R_p - R_s)B_g\right] + W_p B_w}^{\text{Fluidos produzidos}}$$

$$||$$

$$\underbrace{N\left[(Bo - B_{oi}) + (R_{si} - R_s)B_g\right] + mN B_{oi}\left(\frac{B_g}{B_{gi}} -1\right)}_{\text{Expansão dos fluidos do reservatório}}$$

Na forma compacta:

$$\frac{F}{E_o} = N + mN\frac{E_g}{E_o}$$

Então basta determinar a expansão do óleo ($E_o = B_o - B_{oi}$) e do gás natural ($E_g = B_{oi}\left(\frac{B_g}{B_{gi}}-1\right)$) e ajustar ao dados de histórico de produção para estimar $N$ e o tamanho da capa de gás ($m$).

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil-gascap.png?raw=true)

Por exemplo, para encontrar o tamanho da capa de gás, $m$, o STOIIP, N, e o GIIP, G, a partir do histórico de produção podemos ter este tipo de comportamento no gráfico $\frac{E_g}{E_o} \times \frac{F}{E_o}$
"""

# ╔═╡ 3d91d2d1-b199-483a-86c5-bca3121f1f57
let
	p =[4200, 3850, 3708, 3590, 3410, 3300, 2985, 2752, 2500]
	Nₚ = [0, 8.92, 12.023, 13.213, 14.776, 17.268, 20.16, 26.704, 28.204]
	bo = [1.3696, 1.4698, 1.4542, 1.4423, 1.4304, 1.4185, 1.4056, 1.3827, 1.3603]
	bg = [0.00095, 0.00109, 0.00117, 0.00121, 0.00129, 0.00138, 0.0015, 0.00183, 0.00194]
	Rₚ = [0, 1249, 1261, 1370, 1469, 1505, 1547, 1645, 1666]
	Rₛ = [640, 568, 535, 513, 477, 446, 419, 403, 362]

	F = Nₚ .* (bo .+ (Rₚ .- Rₛ) .* bg)
	Eo = (bo .- bo[1]) .+ (Rₛ[1] .- Rₛ) .* bg .+ 1e-9 
	Eg = bo[1].*(bg./bg[1].-1)

	# fitting
	mbal(x, p) = p[1] .+ x .* p[2]
	p0 = [0.5, 0.5]
	fit = curve_fit(mbal, Eg[2:end]./Eo[2:end], F[2:end]./Eo[2:end], p0)
	a = round(fit.param[1], digits=3)
	b = round(fit.param[2], digits=3)
	scatter(Eg[2:end]./Eo[2:end], F[2:end]./Eo[2:end], label="Dados de Produção", xlabel="Eg/Eo", ylabel="F/Eo")
	plot!(Eg./Eo,mbal(Eg./Eo,fit.param),label="MBAL", c=:red, lw=2)
	title!("F = $a + $b Eg/Eo")
	xlims!(0,3)
	ylims!(0,250)
end

# ╔═╡ 9c8f574c-1930-456a-a345-8d7fdc719a5c
md"""
# Reservatório com influxo de água

A partir da MBAL geral, chegamos a:

$$F = N\left[(Bo - B_{oi}) + (R_{si} - R_s)B_g\right] + NB_{oi}\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p_i + W_e$$

Onde temos as seguintes contribuições:

* Fluidos produzidos:

$$F = N_p\left[B_o + (R_p - R_s)B_g\right] + W_p B_w$$

* Expansão dos fluidos do reservatório:

  $$N\left[(Bo - B_{oi}) + (R_{si} - R_s)B_g\right] + NB_{oi}\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p_i$$

* Influxo de água (aquífero associado)

$$W_e$$

Rearranjando a equação MBAL, temos:

$$\frac{F}{E_o + B_{oi}E_{wf}} = N + \frac{W_e}{E_o + B_{oi}E_{wf}}$$

sendo: $E_o = (Bo - B_{oi}) + (R_{si} - R_s)B_g$ e $E_{wf} =\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p_i$, e $W_e$ modelos de influxo de aquífero.

**Obs.**: Se OOIP é conhecido a partir de uma análise volumetríca, o volume do influxo do aquífero pode ser determinado rearranjando o MBAL, para:

$$W_e = B_w W_p + F - NE_o - NB_{oi}E_{fw}$$

## Método de Havlena-Odeh (1963)

![](https://github.com/johnteixeira-ctec-ufal/EPET060-ER2-lectures/blob/main/images/mbal-oil-Havlena-Odeh.png?raw=true)

A análise deste gráfico é bastante utilizando para determinar OOIP, ao mesmo tempo determinar, para o mecânismo de produção, a constante de influxo do aquífero do modelo de VEH:

$$\frac{F}{E_o + B_{oi}E_{wf}} = N + C\frac{\sum\Delta p W_{eD}}{E_o + B_{oi}E_{wf}}$$

Para o modelo de Fetkovich:

$$\frac{F}{E_o + B_{oi}E_{wf}} = N + J\frac{\Delta\bar{p}}{E_o + B_{oi}E_{wf}}$$

Modelo simplificado

$$\frac{F}{E_o + B_{oi}E_{wf}} = N + C\frac{\Delta p_i}{E_o + B_{oi}E_{wf}}$$
"""

# ╔═╡ 853a6608-8c24-4c8a-b959-96414e0fad5d
md"""
## Método de Sills (195)

Para avaliar/diagnosticar a presença de influxo de aquífero no mecanismo de produção, faz-se  necessario a análise do comportamento do gráfico $\frac{F}{E_o + B_{oi}E_{wf}}\times N_p$

![](https://media.springernature.com/lw685/springer-static/image/chp%3A10.1007%2F978-3-030-02393-5_6/MediaObjects/467770_1_En_6_Figc_HTML.png)
"""

# ╔═╡ b8a90732-8731-48f6-8770-43ff15a6700f
details(
	md"""**Exercício 3.** Considere um reservatório, saturado, produzindo sob mecanismos de influxo de água, apresenta os segintes dados PVTs e histórico de produção:

	| P (psia) | Np (stb) | Wp (stb) | Gp (Mscf) | Bo (rb/stb) | Rs (scf/stb) | Bg (rb/scf) | Bw (rb/stb) |
	| -------- | -------- | -------- | --------- | ----------- | ------------ | ------------- | ----------- |
	| 3093 | 0 | 0 | 0 | 1.3101 | 504 | 0.000950 | 1.0334 |
	| 3017 | 200671 | 0 | 98063 | 1.3113 | 504 | 0.000995 | 1.0336 |
	| 2695 | 1322730 | 7 | 814420 | 1.2986 | 470.9 | 0.001133 | 1.0345 |
	| 2640 | 1532250 | 10 | 894484 | 1.2942 | 461.2 | 0.001150 | 1.0346 |
	| 2461 | 2170810 | 29 | 1359270 | 1.2809 | 430.7 | 0.001239 | 1.0350 |
	| 2318 | 2579850 | 63 | 1826800 | 1.2700 | 406.2 | 0.001324 | 1.0353 |
	| 2071 | 3208410 | 825 | 2736410 | 1.2489 | 361.7 | 0.001505 | 1.0359 |
	| 1903 | 3592730 | 11138 | 3401290 | 1.2360 | 331.5 | 0.001663 | 1.0363 |
	| 1698 | 4011570 | 97446 | 4222680 | 1.2208 | 294.6 | 0.001912 | 1.0367 | 
	
	Caracterize a dominância (forte, fraco ou indiferente) do aquífero utilizando o método de Sills, dados adicionais:
	cᵩ = 2.28$\times\ 10^{-6}$ psia$^{-1}$, ϕ = 0.25, Swi = 0.208, cₐ = 3.2$\times\ 10^{-6}$ 
	""",
	md"""
	Para utilizarmos o método de Sills, devemos calculas os termos:
	
	1. Fluidos produzidos: $F = N_p\left[B_o + (R_p - R_s)B_g\right] + W_p B_w$
	1. Expansão dos fluidos do reservatório: $E_o = (Bo - B_{oi}) + (R_{si} - R_s)B_g$ e $E_{wf} =\left(\frac{c_\phi + c_w s_{wi}}{1-s_{wi}}\right)\Delta p$. Aqui precisamos calcular $E_{wf}$ pois o reservatório passou do estado subsaturado para o saturado

	Assim:

	|  $F$  | $E_o$ | $E_{wf}$ |
	| ----- | ----- | ------- |
	| 0 | 0.0 | 0.0 |
	| 162605  | 0.0012 | 0.000282659 |
	| 1012910 | 0.0260023 | 0.00148024 |
	| 1171400 | 0.03332 | 0.00168479 |
	| 1623880 | 0.0616187 | 0.00235053 |
	| 1891430 | 0.0893872 | 0.00288237 |
	| 2265430 | 0.152961 | 0.00380101 |
	| 2477200 | 0.212767 | 0.00442584 |
	| 2746800 | 0.311073 | 0.00518827 |

	Portanto, observando os gráficos a seguir pode-se deduzir que há um suporte de aquífero de fraca dominância.
	"""
)

# ╔═╡ ea8c958a-d285-4f2b-a512-cf2938d2d363
let
	cᵩ = 2.28e-6
	sw = 0.208
	cₐ = 3.2e-6
	p = [3093, 3017, 2695, 2640, 2461, 2318, 2071, 1903,  1698]
	Nₚ = [1e-12, 200671, 1322730, 1532250, 2170810, 2579850, 3208410, 3592730, 4011570]
	Wₚ = [0, 0, 7, 10, 29, 63, 825, 11138, 97446]
	Gₚ = [0, 98063, 814420, 894484, 1359270, 1826800, 2736410, 3401290, 4222680]
	bo = [1.3101, 1.3113, 1.2986, 1.2942, 1.2809, 1.2700, 1.2489, 1.2360, 1.2208]
	Rso =[504, 504, 470.9, 461.2, 430.7, 406.2, 361.7, 331.5, 294.6]
	bg = [0.000950, 0.000995, 0.001133, 0.001150, 0.001239, 0.001324, 0.001505, 0.001663, 0.001912]
	bw = [1.0334, 1.0336, 1.0345, 1.0346, 1.0350, 1.0353, 1.0359, 1.0363, 1.0367]
	# computations
	Eo = bo .- bo[1] .+ (Rso[1] .- Rso) .* bg
	Ewf = (cᵩ + cₐ * sw) / (1 - sw) .* (p[1] .- p)
	F = Nₚ .* (bo .+ (Gₚ./Nₚ .- Rso) .* bg) .+ Wₚ .* bw
	y = F ./ (Eo .+ bo[1] .* Ewf); y[1]=0
	plot(Nₚ, y, label="", lw=3, xlabel="Nₚ", ylabel="F / (Eₒ + BₒᵢEᵩₐ)", title="Método de Sills (1996)")
	scatter!(Nₚ, y, label="", xlabel="Nₚ", ylabel="F / (Eₒ + BₒᵢEᵩₐ)", title="Método de Sills (1996)")	
end

# ╔═╡ f4d19003-3654-4d4a-9e56-b659c730bc4d
md"""
!!! info "Diagnóstico preliminar"
	O método de Sills é uma importante ferramenta que deve ser utilizada sempre antes de iniciar uma análise de MBAL. Desta forma, realize esta análise para os demais exercícios.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
MultiComponentFlash = "35e5bd01-9722-4017-9deb-64a5d32478ff"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
HypertextLiteral = "~0.9.5"
Interpolations = "~0.15.1"
LsqFit = "~0.15.0"
MultiComponentFlash = "~1.1.16"
NLsolve = "~4.5.1"
Plots = "~1.40.9"
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.8"
manifest_format = "2.0"
project_hash = "82544714f2c2222e5c9d52dc403f872afca0caf5"

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

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

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

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

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

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

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

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "030ea22804ef91648f29b7ad3fc15fa49d0e6e71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "5e1897147d1ff8d98883cda2be2187dcf57d8f0c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.15.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

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

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

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

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

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
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

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
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

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
# ╟─286c1e61-8d3c-4570-8611-985f4e88fe56
# ╟─369a85ef-d489-4f2c-83aa-0348897aafda
# ╟─ce8cf761-58df-4769-af0b-b71efdf7d2f6
# ╟─785e1863-8001-4885-8829-9f98419d6cd9
# ╟─c1cb7823-f336-4c66-b92e-de78bae8aad0
# ╟─5c641485-6bea-498e-874b-c0fc0ed99743
# ╟─856145bf-8853-4d07-81dd-fcecbc45c422
# ╟─b4af5d4e-01d4-491a-9e9e-40eabbd02c22
# ╟─7a2c0627-5a20-42f4-b99f-ef9a37b00aed
# ╟─f22d739d-a5b4-4543-8a75-7612cf6b44bb
# ╟─e9946302-c337-402d-984e-ab14e3932469
# ╟─6532b23d-0a83-4cac-82ce-166dfef2ffb1
# ╟─e6ed5d77-eba5-422a-93a6-bdb71391d939
# ╟─d81ddba8-30c4-458e-a1ff-d2e30d03f402
# ╟─b07ac32c-9017-4103-8614-3d9203413adb
# ╟─daf6ba60-2c2b-4097-8283-3d661d1b247b
# ╟─5e3cd049-335d-4384-a903-87659e52587a
# ╟─44ff5e94-08a4-4df0-9361-5f1cd924080a
# ╟─87c10293-85e0-4b30-8af6-45588fb3d92a
# ╟─8caa9071-4f3b-4d79-8164-6e3e4aee580b
# ╟─7f033f77-1e6a-4954-bd43-ed09a10c8ade
# ╟─3d91d2d1-b199-483a-86c5-bca3121f1f57
# ╟─9c8f574c-1930-456a-a345-8d7fdc719a5c
# ╟─853a6608-8c24-4c8a-b959-96414e0fad5d
# ╟─b8a90732-8731-48f6-8770-43ff15a6700f
# ╟─ea8c958a-d285-4f2b-a512-cf2938d2d363
# ╟─f4d19003-3654-4d4a-9e56-b659c730bc4d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
